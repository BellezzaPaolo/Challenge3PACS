#include "problem.hpp"
#include <iostream>


int main(){
    problem P(0.0,1.0,0.0,1.0,32);

    P.SeqSolver(1e-7);
    //P.printSol();

    std::cout<<P.computeError()<<std::endl;

    return 0;
}