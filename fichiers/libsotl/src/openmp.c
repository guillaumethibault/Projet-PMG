#include "default_defines.h"
#include "global_definitions.h"
#include "device.h"
#include "openmp.h"
#include "sotl.h"
#include <omp.h>
#include <string.h>

#ifdef HAVE_LIBGL
#include "vbo.h"
#endif

#include <stdio.h>
#define NUM_THREADS 8
static int *atom_state = NULL;
static calc_t z_atoms_distance(sotl_atom_set_t* dev, unsigned pos_atom_1 , unsigned pos_atom_2);

struct {
  float r , g, b;
} colors[NUM_THREADS] = {{1.0 , 0.0 , 0.0},
			 {0.1 , 1.0 , 0.0},
			 {0.0 , 0.0 , 1.0},
			 {1.0 , 0.0 , 1.1},
			 {1.0 , 1.0 , 1.0},
			 {0.0 , 0.0 , 0.0},
			 {0.0 , 1.0 , 0.0},
			 {1.0 , 1.0 , 1.0}};

#ifdef HAVE_LIBGL

#define SHOCK_PERIOD 50

// Update OpenGL Vertex Buffer Object
//
static void omp_update_vbo (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  sotl_domain_t *domain = &dev->domain;
#pragma omp for schedule(runtime)
  for (unsigned n = 0; n < set->natoms; n++) {
    vbo_vertex[n*3 + 0] = set->pos.x[n];
    vbo_vertex[n*3 + 1] = set->pos.y[n];
    vbo_vertex[n*3 + 2] = set->pos.z[n];

    // Atom color depends on z coordinate
    {
      float ratio = (set->pos.z[n] - domain->min_ext[2]) / (domain->max_ext[2] - domain->min_ext[2]);
      int tid = atom_state[n];
      vbo_color[n*3 + 0] = (1.0 - ratio) * atom_color[0].R + ratio * 1.0 + colors[tid].r;
      vbo_color[n*3 + 1] = (1.0 - ratio) * atom_color[0].G + ratio * 0.0 + colors[tid].g;
      vbo_color[n*3 + 2] = (1.0 - ratio) * atom_color[0].B + ratio * 0.0 + colors[tid].b;
      //   atom_state[n]--;
    }
  }
}
#endif

// Update positions of atoms by adding (dx, dy, dz)
//
static void omp_move (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;

#pragma omp for schedule(runtime)
  for (unsigned n = 0; n < set->natoms; n++) {
    set->pos.x[n] += set->speed.dx[n];
    set->pos.y[n] += set->speed.dy[n];
    set->pos.z[n] += set->speed.dz[n];
  }
}

// Apply gravity force
//
static void omp_gravity (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  const calc_t g = 0.005;
  //TODO
  sotl_atom_speed_t* atom_speed = &set-> speed;
  unsigned nb_atoms = set-> natoms;

#pragma omp for schedule(runtime)
  for (unsigned i = 0 ; i < nb_atoms ; i++)
    atom_speed-> dy[i] -= g;
}

static void omp_bounce (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  sotl_domain_t *domain = &dev->domain;
  for (unsigned c = 0; c < 3; c++) {// c == 0 -> x, c == 1 -> y, c == 2 -> z
    calc_t *pos = set->pos.x + c * set->offset;
    calc_t *speed = set->speed.dx + c * set->offset;
 
#pragma omp for //simd safelen(8)
    for (unsigned n = 0; n < set->natoms; n++){
      if (pos[n] < domain->min_ext[c] || pos[n] > domain->max_ext[c])
	{
	  speed[n] *= -1;	  
	  //	  speed[n] /= (LENNARD_CUTOFF)/2;
	}
    }
  }
}

static calc_t squared_distance (sotl_atom_set_t *set, unsigned p1, unsigned p2)
{
  calc_t *pos1 = set->pos.x + p1,
    *pos2 = set->pos.x + p2;

  calc_t dx = pos2[0] - pos1[0],
    dy = pos2[set->offset] - pos1[set->offset],
    dz = pos2[set->offset*2] - pos1[set->offset*2];

  return dx * dx + dy * dy + dz * dz;
}

static calc_t lennard_jones (calc_t r2)
{
  calc_t rr2 = 1.0 / r2;
  calc_t r6;

  r6 = LENNARD_SIGMA * LENNARD_SIGMA * rr2;
  r6 = r6 * r6 * r6;

  return 24 * LENNARD_EPSILON * rr2 * (2.0f * r6 * r6 - r6);
}

static calc_t z_distance (sotl_atom_set_t *set, unsigned p1, unsigned p2)
{
  calc_t */*restrict*/ pos1 = set->pos.x + p1,
    * /*restrict*/ pos2 = set->pos.x + p2;

  calc_t  dz = pos2[set->offset*2] - pos1[set->offset*2];

  return  dz * dz;
}

// Calcul de la force avec le tri selon l'axe z
static void omp_force_z (sotl_device_t *dev)
{
  sotl_atom_set_t */*restrict*/ set = &dev->atom_set;

#pragma omp for 
  for (unsigned current = 0; current < set->natoms; current++) {
    calc_t force[3] = { 0.0, 0.0, 0.0 };

    //#pragma omp simd 
    for (unsigned other = current-1; other < set->natoms; other--)
      {

	if (z_distance(set, current, other) > LENNARD_SQUARED_CUTOFF)
	  break;
	calc_t sq_dist = squared_distance (set, current, other);

	if (sq_dist < LENNARD_SQUARED_CUTOFF) {
	  calc_t intensity = lennard_jones (sq_dist);
	  
	  calc_t * /*restrict*/ posx = set->pos.x ;

	  force[0] += intensity * (posx[current] - posx[other]);
	  force[1] += intensity * (posx[set->offset + current] -
				   posx[set->offset + other]);
	  force[2] += intensity * (posx[set->offset * 2 + current] -
				   posx[set->offset * 2 + other]);
	}

      }
    //#pragma omp simd 
    for (unsigned other = current + 1; other < set->natoms; other++)
      {	
	if (z_distance(set, current, other) > LENNARD_SQUARED_CUTOFF)
	  break;
	
	
	calc_t sq_dist = squared_distance (set, current, other);

	if (sq_dist < LENNARD_SQUARED_CUTOFF) {
	  calc_t intensity = lennard_jones (sq_dist);
	  
	  calc_t * /*restrict*/ posx = set->pos.x ;

	  force[0] += intensity * (posx[current] - posx[other]);
	  force[1] += intensity * (posx[set->offset + current] -
				   posx[set->offset + other]);
	  force[2] += intensity * (posx[set->offset * 2 + current] -
				   posx[set->offset * 2 + other]);
	}

      }   
    set->speed.dx[current] += force[0];
    set->speed.dx[set->offset + current] += force[1];
    set->speed.dx[set->offset * 2 + current] += force[2];
  }
}

// Calcul de la force avec le tri par boite
static void omp_force (sotl_device_t *dev)
{
  sotl_atom_set_t *set = &dev->atom_set;
  sotl_domain_t *dom = &dev->domain;

  int *sum_prefix = atom_set_prefix_sum(dom, set);

#pragma omp for schedule(runtime)
  for (unsigned current = 0; current < set->natoms; current++) {
    int tid = omp_get_thread_num();
    atom_state[current] = tid % NUM_THREADS;
    calc_t force[3] = { 0.0, 0.0, 0.0 };

    int center = atom_get_num_box(dom, set->pos.x[current], set->pos.y[current], set->pos.z[current], BOX_SIZE_INV);
    int boxes_to_test[27];
    int atoms_to_test[270];
    memset(boxes_to_test , '\0' , 27);

    // Mise en place des bornes pour ne pas récupérer de boite inexistante
    int minX = -1, maxX = 1, minY = -1, maxY = 1, minZ = -1, maxZ = 1;
    if(set->pos.x[current] <= dom->min_ext[0])
      minX=0;
    else if(set->pos.x[current] >= dom->max_ext[0])
      maxX=0;
    if(set->pos.y[current] <= dom->min_ext[1])
      minY=0;
    else if(set->pos.y[current] >= dom->max_ext[1])
      maxY=0;
    if(set->pos.z[current] <= dom->min_ext[2])
      minZ=0;
    else if(set->pos.z[current] >= dom->max_ext[2])
      maxZ=0;

    // Calcul des boites à tester
    int cpt = 0;
    
    for(int z=minZ; z<=maxZ; z++)
      for(int y=minY; y<=maxY; y++)
	for(int x=minX; x<=maxX; x++)
	  {
	    boxes_to_test[cpt] = center + dom->boxes[0] * dom->boxes[1] * z + dom->boxes[1] * y + x; 
	    cpt++;
	  }

#pragma omp parallel for schedule(runtime)
    for(int i = cpt; i < 27; i++)
      boxes_to_test[i] = -1;

    // Calcul des atomes à tester
    cpt = 0;
    
    for(int i=0; i<27; i++)
      {
      	if(boxes_to_test[i] == -1)
	  continue;

	for(int j=sum_prefix[boxes_to_test[i]]; j<sum_prefix[boxes_to_test[i]+1]; j++)
	  {
	    atoms_to_test[cpt] = j;
	    cpt++;
	  }
      }

#pragma omp parallel for schedule(runtime)
    for(int i = cpt; i < 270; i++)
      atoms_to_test[i] = -1;

    /*      for(int i=0; i<27; i++)
	    printf("box %d = %d\n", i, boxes_to_test[i]);

	    for(int i=0; i<120; i++)
	    printf("atom %d = %d\n", i, atoms_to_test[i]); */
  
    for(unsigned i = 0; i < 270; i++)
      {
	unsigned other = atoms_to_test[i];

	if(other == current)
	  continue;
	if((int)other == -1)
	  break;
	
	calc_t sq_dist = squared_distance (set, current, other);
	
	if (sq_dist < LENNARD_SQUARED_CUTOFF)
	  {
	    calc_t intensity = lennard_jones (sq_dist);
	    
	    calc_t * /*restrict*/ posx = set->pos.x ;
	    
	    force[0] += intensity * (posx[current] - posx[other]);
	    force[1] += intensity * (posx[set->offset + current] -
				     posx[set->offset + other]);
	    force[2] += intensity * (posx[set->offset * 2 + current] -
				     posx[set->offset * 2 + other]);	
	  }
      }
    set->speed.dx[current] += force[0];
    set->speed.dx[set->offset + current] += force[1];
    set->speed.dx[set->offset * 2 + current] += force[2];
  
  }

  free(sum_prefix);
}

// Main simulation function
//
void omp_one_step_move (sotl_device_t *dev)
{
  // Tri par boite
  box_set_sort(&dev->domain , &dev->atom_set);

  // Tri selon l'axe z
  // atom_set_sort(&dev-> atom_set);

#pragma omp parallel firstprivate(dev)
  { // Apply gravity force
    //
    if (gravity_enabled)
      omp_gravity (dev);

    // Compute interactions between atoms
    //
    if (force_enabled)
      // Force avec le tri par boite
      omp_force (dev);
    // Force avec me tri selon l'axe z
    // omp_force_z (dev);
    
    // Bounce on borders
    //
    if(borders_enabled)
      omp_bounce (dev);

    // Update positions
    //
    omp_move (dev);

#ifdef HAVE_LIBGL
    // Update OpenGL position
    //
    if (dev->display)
      omp_update_vbo (dev);
#endif
  }
}

void omp_init (sotl_device_t *dev)
{
#ifdef _SPHERE_MODE_
  sotl_log(ERROR, "Sequential implementation does currently not support SPHERE_MODE\n");
  exit (1);
#endif

  borders_enabled = 1;

  dev->compute = SOTL_COMPUTE_OMP; // dummy op to avoid warning
}

void omp_alloc_buffers (sotl_device_t *dev)
{
  atom_state = calloc(dev->atom_set.natoms, sizeof(int));
  printf("natoms: %d\n", dev->atom_set.natoms);
}

void omp_finalize (sotl_device_t *dev)
{
  free(atom_state);

  dev->compute = SOTL_COMPUTE_OMP; // dummy op to avoid warning
}
