#include "problem.hpp"
#include <iostream>

std::vector<double> problem::SeqSolver(int N){
    double h=(D.x1-D.x0)/(N-1);
    std::vector<double> U(N*N,0.0);

    return U;
}

void problem::printSol() const{
    for(size_t i=0;i<N;++i){
        for(size_t j=0;j<N;++j){
            std::cout<<Uapproximate[i*N+j]<<"  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<"\n"<<std::endl;
    return;
}