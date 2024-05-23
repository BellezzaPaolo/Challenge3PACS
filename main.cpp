#include "problem.hpp"
#include <iostream>


int main(int argc, char *argv[]){
    if(argc!=3){
        std::cerr<<"Error:\ninput must be 2: right hand side and number of elements"<<std::endl;
        exit(0);
    }
    if(atoi(argv[2])<0){
        std::cerr<<"Error:\nnumber of elements must be positive"<<std::endl;
        exit(0);
    }

    problem P(0.0,1.0,0.0,1.0,atoi(argv[2]),argv[1]);
    std::ofstream file {"Output.vtk"};

    //P.SeqSolver(1e-7);
    //P.printSol(file);

    //std::cout<<P.computeError()<<std::endl;

    file.close();

    return 0;
}