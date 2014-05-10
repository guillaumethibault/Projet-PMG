#define _XOPEN_SOURCE 600

#include <errno.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "atom.h"
#include "sotl.h"
#include "default_defines.h"
#include "ocl.h"
#include "seq.h"
#include "ocl_kernels.h"
#include "window.h"
#include "profiling.h"

#ifdef HAVE_LIBGL
#include "vbo.h"
#include "shaders.h"
#endif

static sotl_atom_set_t global_atom_set;

static void sift_down(sotl_atom_set_t *set, const int start, const int count);

#define SWAP(a, b)				\
  do {						\
    calc_t t = a; a = b; b = t;			\
  } while(0)

#define SWAP_INT(array, i , j)			\
  do {						\
    int t = array[i]; array[i] = array[j];	\
    array[j] = t;				\
  } while(0)

#define SWAP_ATOM(a, b)				\
  do {						\
    SWAP(set->pos.z[a], set->pos.z[b]);		\
    SWAP(set->pos.y[a], set->pos.y[b]);		\
    SWAP(set->pos.x[a], set->pos.x[b]);		\
    SWAP(set->speed.dz[a], set->speed.dz[b]);	\
    SWAP(set->speed.dy[a], set->speed.dy[b]);	\
    SWAP(set->speed.dx[a], set->speed.dx[b]);	\
  } while(0)

static void heap_sort(sotl_atom_set_t *set, const unsigned count)
{
  int start, end;

  /* heapify */
  start = (count - 2) / 2;
  for (; start >= 0; start--) {
    sift_down(set, start, count);
  }

  for (end = count - 1; end > 0; end--) {
    SWAP_ATOM(end, 0);
    sift_down(set, 0, end);
  }
}

static void sift_down(sotl_atom_set_t *set, const int start, const int end)
{
  int root = start;

  while (root * 2 + 1 < end) {
    int child = 2 * root + 1;
    if ((child + 1 < end) && (set->pos.z[child] < set->pos.z[child + 1])) {
      child += 1;
    }
    if (set->pos.z[root] >= set->pos.z[child])
      return;

    SWAP_ATOM(child, root);
    root = child;
  }
}

int atom_set_init(sotl_atom_set_t *set, const unsigned long natoms,
                  const unsigned long maxatoms)
{
  if (maxatoms < natoms)
    return SOTL_INVALID_VALUE;

  set->natoms = natoms;
  set->current = 0;
  set->offset = ROUND(maxatoms);

  if (!sotl_have_multi()) {
    /* No need to have ghosts for a single device. */
    /* XXX: This is absolutely wrong when the torus mode is enabled. :-) */
    set->offset_ghosts = 0;
  } else {
    set->offset_ghosts = ROUND(set->natoms * 0.05);
  }

  /* No ghosts at the beginning. */
  set->nghosts_min = set->nghosts_max = 0;

  set->pos.x = malloc(atom_set_size(set));
  if (!set->pos.x)
    return SOTL_OUT_OF_MEMORY;
  set->pos.y = set->pos.x + set->offset;
  set->pos.z = set->pos.y + set->offset;
  set->speed.dx = malloc(atom_set_size(set));
  if (!set->speed.dx) {
    atom_set_free(set);
    return SOTL_OUT_OF_MEMORY;
  }
  set->speed.dy = set->speed.dx + set->offset;
  set->speed.dz = set->speed.dy + set->offset;
#pragma omp parallel for schedule(runtime)
  for (unsigned i = 0 ; i < set->natoms ; i++)
    {
      set->pos.x[i] = 0;
      set->speed.dx[i] = 0;
      set->pos.y[i] = 0;
      set->speed.dy[i] = 0;
      set->pos.z[i] = 0;
      set->speed.dz[i] = 0;
    }
  return SOTL_SUCCESS;
}

void atom_set_print(const sotl_atom_set_t *set)
{
  sotl_log(DEBUG, "natoms = %d, current = %d, offset = %d\n",
	   set->natoms, set->current, set->offset);
  sotl_log(DEBUG, "pos.x = %p, pos.y = %p, pos.z = %p\n", set->pos.x,
	   set->pos.y, set->pos.z);
}

sotl_atom_set_t *get_global_atom_set()
{
  return &global_atom_set;
}

int atom_set_add(sotl_atom_set_t *set, const calc_t x, const calc_t y,
                 const calc_t z, const calc_t dx, const calc_t dy,
                 const calc_t dz)
{
  if (set->current >= set->natoms)
    return SOTL_INVALID_BUFFER_SIZE;

  set->pos.x[set->current] = x;
  set->pos.y[set->current] = y;
  set->pos.z[set->current] = z;

  set->speed.dx[set->current] = dx;
  set->speed.dy[set->current] = dy;
  set->speed.dz[set->current] = dz;

  set->current++;
  return SOTL_SUCCESS;
}

size_t atom_set_offset(const sotl_atom_set_t *set)
{
  return set->offset + set->offset_ghosts * 2;
}

size_t atom_set_size(sotl_atom_set_t *set)
{
  return sizeof(calc_t) * set->offset * 3;
}

size_t atom_set_border_size(const sotl_atom_set_t *set)
{
  return sizeof(calc_t) * set->offset_ghosts * 3;
}

size_t atom_set_begin(const sotl_atom_set_t *set)
{
  return set->offset_ghosts;
}

size_t atom_set_end(const sotl_atom_set_t *set)
{
  return atom_set_begin(set) + set->natoms;
}

#define FREE(x)					\
  do {						\
    free(x);					\
    x = NULL;					\
  } while (0)

void atom_set_free(sotl_atom_set_t *set)
{
  FREE(set->pos.x);
  FREE(set->speed.dx);
}

void atom_set_sort(sotl_atom_set_t *set)
{
  /* Sort atoms along z-axis. */
  heap_sort(set, set->natoms);
}

#ifdef HAVE_LIBGL
void atom_build (int natoms, sotl_atom_pos_t * pos_vec)
{
  int i;
  for (i = 0; i < natoms; i++)
    vbo_add_atom (pos_vec->x[i], pos_vec->y[i], pos_vec->z[i]);
}
#endif

int atom_get_num_box(const sotl_domain_t *dom, const calc_t x, const calc_t y,
                     const calc_t z, const calc_t rrc)
{
  int box_x, box_y, box_z;
  int box_id;

  box_x = (x - dom->min_border[0]) * rrc;
  box_y = (y - dom->min_border[1]) * rrc;
  box_z = (z - dom->min_border[2]) * rrc;

  box_id = box_z * dom->boxes[0] * dom->boxes[1] +
    box_y * dom->boxes[0] +
    box_x;

  assert(box_id >= 0 && (unsigned)box_id < dom->total_boxes);
  return box_id;
}

int *atom_set_box_count(const sotl_domain_t *dom, const sotl_atom_set_t *set)
{
  int *boxes = NULL;
  size_t size;

  size = dom->total_boxes * sizeof(int);
  if (!(boxes = calloc(1, size)))
    return NULL;

#pragma omp parallel for schedule(runtime)
  for (unsigned i = 0; i < set->natoms; i++) {
    int box_id = atom_get_num_box(dom, set->pos.x[i], set->pos.y[i], set->pos.z[i],
				  BOX_SIZE_INV);
#pragma omp atomic
    boxes[box_id]++;
  }

  return boxes;
}

int* atom_set_prefix_sum(const sotl_domain_t *dom, const sotl_atom_set_t *set)
{
  size_t size;
  int* boxes = atom_set_box_count(dom , set);
  int* prefix_array_sum = NULL;

  size = (dom->total_boxes + 1) * sizeof(int);
  if (!(prefix_array_sum = calloc(1, size)))
    return NULL;

  prefix_array_sum[0] = 0;
  for (unsigned i = 1 ; i <= dom->total_boxes ; i++)
    prefix_array_sum[i] = prefix_array_sum[i - 1] + boxes[i - 1];

  free(boxes);

  /* safe guard */
  assert((unsigned)prefix_array_sum[dom->total_boxes] == set->natoms);
  return prefix_array_sum;
}

void box_set_sort(const sotl_domain_t *dom, sotl_atom_set_t *set)
{
  // tableau contenant l'id de la boite de chaque atome
  int num_box_tab[set->natoms];
  int *sum_prefix = atom_set_prefix_sum(dom, set);
  
  int cpt = 0;

#pragma omp parallel for schedule(runtime)
  for (unsigned i = 0; i < set->natoms; i++)
    num_box_tab[i] = atom_get_num_box(dom, set->pos.x[i], set->pos.y[i], set->pos.z[i],
				      BOX_SIZE_INV);

  for(unsigned i = 0; i < dom->total_boxes; i++)
    for(unsigned j = 0+cpt; j < set->natoms; j++)
      {
	if(num_box_tab[j] == (int)i)
	  {
	    SWAP_ATOM(cpt, j);
	    SWAP_INT(num_box_tab , cpt , j);
	    cpt++;
	  }
	if(sum_prefix[i + 1] == cpt)
	  break;
      }
}
