#include "problem.hpp"
#include <iostream>


int main(int argc, char *argv[]){
    if(argc!=4){
        std::cerr<<"Error:\ninput must be 2: right hand side and number of elements"<<std::endl;
        exit(0);
    }
    if(atoi(argv[2])<0){
        std::cerr<<"Error:\nnumber of elements must be positive"<<std::endl;
        exit(0);
    }


    problem P(0.0,1.0,0.0,1.0,atoi(argv[2]),argv[1]);
    std::ofstream file {"Output.vtk"}; 
    /*P.SeqSolver(1e-7);
    std::cout<<"L'errore è "<<P.computeError()<<std::endl;
    P.EraseSol();*/
    MPI_Init(&argc,&argv);
    std::vector<double> Ulocal;
    P.split(Ulocal);
    P.ParSolver(1e-7,Ulocal);

    //P.printSol(file);

    

    MPI_Finalize();
    //std::cout<<"L' errore è "<<P.computeError()<<std::endl;
    file.close();

    return 0;
}