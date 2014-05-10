#! /bin/bash

declare -a wg_size=("8" "16" "32" "64" "128")

help()
{
    echo "bash_mode.sh <NR_MAX_PROC> <NR_ATOMS> <NR_ITERATIONS> <CONF>"
    exit 1
}

if [ $# -lt 2 ]
then
    help
fi

#NR_PROC=$1
#NR_ATOMS=$1
NR_ITERATIONS=$1
CONF=$2
time=$(bin/atoms -v -s 1 -i $NR_ITERATIONS $CONF) #-n $NR_ATOMS 
seq_time=`awk '$1 == "[PERF]" && $2 ~ /[0-9]+/  {print $3}' <<EOF
$time
EOF`
printf "%s\n" $seq_time
#for i in `seq 1 $NR_PROC`
for i in ${wg_size[*]}
do
#    export OMP_PROC_BIND=true
#    export OMP_SCHEDULE=static,1
#    export GOMP_CPU_AFFINITY="0-$NR_PROC"
#    export OMP_NUM_THREADS=$i
    time=$(bin/atoms -v -d 1 -i $NR_ITERATIONS -w $i $CONF) #-n $NR_ATOMS
    omp_time=`awk '$1 == "[PERF]" && $2 ~ /[0-9]+/  {print $3}' <<EOF
$time
EOF`
    printf "%d %s\n" $i $omp_time 

done
