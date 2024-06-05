make optimize

Mpi=(2 3 4)

N=(8 16 32 64)

omp=(2 8 16)

for x in "${Mpi[@]}"
do
    for y in "${omp[@]}"
    do
        for z in "${N[@]}"
        do
            mpiexec -n $x './main' '8*pi^2*sin(2*pi*x)*sin(2*pi*y)' $z $y
        done
    done
done