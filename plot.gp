#! /usr/bin/gnuplot

set terminal pdf
#set yrange [0 : 25000]
set xlabel "Nb Workgroups"
set ylabel "Temps"

plot "perf/perf_203_opencl_choc3" using ($1):(47291/($2)) title 'choc3 en 203, avec 8 à 128 workgroups' with linespoints, \
     "perf/perf_203_opencl_opt_choc3" using ($1):(30927/($2)) title 'choc3 en 203, version optimisé avec 8 à 128 workgroups' with linespoints, \
     "perf/perf_203_opencl_choc2" using ($1):(20917/($2)) title 'choc2 en 203, avec 8 à 128 workgroups' with linespoints, \
     "perf/perf_203_opencl_opt_choc2" using ($1):(14321/($2)) title 'choc2 en 203, version optimisé avec 8 à 128 workgroups' with linespoints