# tell OpenMP how many threads it should use
# I set it to the number of cores of my CPU. 
export OMP_NUM_THREADS=$(cat /proc/cpuinfo | grep processor | wc -l )

mkdir -p output
