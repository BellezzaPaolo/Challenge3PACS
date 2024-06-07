# Remove the old report file and old matlab script
rm -f report.txt
rm -f plotData.m

cd ..
# Remove the old executable and object file
make clean
# Recompile the executable with optimization flag
make optimize


# Inizialise all possible value of mpi size, dimension of the mesh and number of Open MP threads
Mpi=(1 2 3 4)

N=(8 16 32 64 128 256)

omp=(2 8 16)

# This variable is used only to switch on or off the use of open MP parallelization
y=1

# Make all the simulations with every combination of above values
for x in "${Mpi[@]}"
do
    #for y in "${omp[@]}"
    #do
        printf "MPI processes: %s Open MP thread %2s:" "$x" "$y">>./test/report.txt
        for z in "${N[@]}"
        do
            # Run the program with the needed input (in order: right hand side, number of element for every row and column of the domain and number of open MP threads)
            exec=$(mpiexec -n $x './main' '8*pi^2*sin(2*pi*x)*sin(2*pi*y)' $z $y '0*x*y' 'sin(2*pi*x)*sin(2*pi*y)')
            
            # Extract error
	        error=$(echo "$exec" | grep -oP '(?<= error is )[\d.]+([eE][-+]?\d+)?')

            # Extract number of iterations
	        iter=$(echo "$exec" | grep -oP 'took \K\d+')

	        # Extract time
	        total_time=$(echo "$exec" | grep -oP 'needed \K[^ ]+(?: s)?')
	        # Remove ' s'
	        total_time=${total_time}

            # Save above values in the report file
	        printf "\n\t on a mesh %3s x %3s  Iteration count: %6s \t total_time %10s ms \t error: %8s" "$z" "$z" "$iter" "$total_time" "$error" >>./test/report.txt
            # Save above values in a Matlab file to make the plots
            printf "%s,%s,%s,%s;" "$x" "$z" "$error" "$total_time">>./test/plotData.m

        done
        printf "\n\n">>./test/report.txt
    #done
    printf "\n\n">>./test/report.txt
done
