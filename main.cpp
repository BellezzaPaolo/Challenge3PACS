#include "problem.hpp"
#include "chrono"


int main(int argc, char *argv[]){
    MPI_Init(&argc,&argv);//initialization of MPI

    //check if the right input are given else the program stops
    if(argc!=6){
        std::cerr<<"Error:\ninput must be 3: right hand side,number of elements and number of parallel tasks"<<std::endl;
        exit(0);
    }
    if(atoi(argv[2])<0){
        std::cerr<<"Error:\nnumber of elements must be positive"<<std::endl;
        exit(0);
    }
    if(atoi(argv[3])<=0){
        std::cerr<<"Error:\nnumber of threads must be a positive integer"<<std::endl;
        exit(0);
    }

    int rank,size; //respectively are the number of this rank and the total number of MPI rank
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    std::vector<std::string> c; //for the boundary conditions
    c.push_back("d");
    problem P(0.0,1.0,0.0,1.0,atoi(argv[2]),argv[4],c,argv[1],argv[5]); //initailise the problem class
    std::ofstream file {"./test/Output.vtk"}; //name of the file for the print
    int iter=0; //number of iterate
    P.assignBoundaryCondition(); //assign the boundary conditions

    if(size==1){
        const auto start= std::chrono::steady_clock::now(); //start record the time

        int iter=P.SeqSolver(1e-7,1e5); //solves the problem sequentially

        const auto end= std::chrono::steady_clock::now(); //end record the time

        P.printPerformance(size, atoi(argv[2]),iter,std::chrono::duration<double>(end-start).count()); //print the preformance
        P.printSol(file); //print the solution in the vtk file
        
    }
    else{
        const auto start= std::chrono::steady_clock::now(); //start record the time

        iter=P.ParSolver(1e-7,1e5,atoi(argv[3])); //solves the problem in parallel

        const auto end= std::chrono::steady_clock::now(); //end record of time
        MPI_Barrier(MPI_COMM_WORLD); //sincronization of every MPI rank
        
        if(rank==0){
            P.printPerformance(size, atoi(argv[2]),iter,std::chrono::duration<double>(end-start).count()); //print the performance
            P.printSol(file); //print the solution in the vtk file
        }    
    }
  
    MPI_Finalize(); //end of the parallel
    file.close(); //close the file

    return 0;
}