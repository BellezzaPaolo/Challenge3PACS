# Challenge3PACS

This program computes the solution of the laplace differential problem with Dirichelet boundary condition with hibryd parallelism.
It uses MPI and Open MP. 

## Usage
To compile the code there are 2 ways:
-one is the classical:
```bash
make 
```
-the other one is to optimize compilation:
```bash
make optimize
```
and then:
```bash
mpiexec -n x ./main f y z bc uex
```
where:
-x is an integer that contains the number of MPI processes
-f is a string that contains the right hand side of the problem 
-y is an integer that contains the numebr of element in every row and column of the mesh 
-z is an integer that contains the number of open MP threads
-bc is a string that contains the boundary condtions (possibly as a function defined piecewise)
-uex is a string that contains the exact solution
To solve the proposed problem right hand side is '8* pi^2* sin(2* pi* x)* sin(2* pi* y)',the boundary is'0* x* y' 
and the exact solution is 'sin(2* pi* x)* sin(2* pi* y)'.

The makefile also has one target to clean all the object and executable:
```bash
make clean
```
It's also possible to start a test with the the command:
```bash
bash ./test/test.sh
```
That runs a test with:
- 1,2,3 and 4 MPI rank
- 8,16,32,64,128 and 256 element for every row and column of the mesh
- 2,8 and 16 open MP threads 
and every possible ombination of the above combinations. Saves the result in the report.txt file and in a matlab file called DataPlot.

The Matlab script "scriptPlot" imports the data in the file DataPlot and plot the error in relation to the refinement and also the time in relation to the refinement of the mesh.

## hardware Information
The hardware informations can be found in the file hw.info in the folder called test.

## Performance
The performance of the code can be valueted in many ways:
-file "Output.vtk" can be opened with ParaView and shows a 3D plot of the approximated solution computed by the program
-the Matlab script "scriptPlot" imports in Matlab WorkSpace the data obatined with test.sh. With that makes 2 plot in loglog scale:
    -one has the length of every cell (h) along the horizonatal axe and the error along the vertical axe. Shows how the error goes in the sequential code versus the parallel one (with 4 MPI cores)
    -one has the length of every cell (h) along the horizonatal axe and the time for the run along the vertical axe. Shows how the computation took time in the sequential code versus the parallel one (with 4 MPI cores)
-the file "report.txt" reports all the statistics in the computation like error, number of iterations and time to compute and how they change in relation to the mesh size, the number MPI rank and open MP threads.

As expected MPI parallelization decrease the time of the computation, the number of iterations and the error too.


## To improve
The assignation of boundary condition could be improved for example taking into account different type of boundary condition.
Also adding open MP parallelization penalise the velocity of the program (I don't understand why) so commneting this line of the code could be removed:
comment lines 211, 217, 227 and 229 of the file problem.cpp and lines 25, 26 and 51 of test.sh while uncomment line 20 of test.sh
So is possible to switch on and off open MP.