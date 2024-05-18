#include "problem.hpp"
#include <iostream>


int main(){
    std::cout<<"A";
    problem P(0.0,1.0,0.0,1.0);

    P.SeqSolver(10);
    P.printSol();

    return 0;
}