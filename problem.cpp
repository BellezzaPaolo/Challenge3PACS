#include "problem.hpp"

std::vector<double> problem::SeqSolver(double toll){
    double h=(D.x1-D.x0)/(N-1);
    std::vector<double> U1(N*N,0.0),U0(N*N,0.0);
    double err=0;
    int iter=0;

    do{
        err=0;
        iter++;
        for(size_t i=1;i<N-1;++i){
            for(size_t j=1;j<N-1;++j){
                U1[i*N+j]= 0.25*(U0[(i-1)*N+j]+U0[(i+1)*N+j]+U0[i*N+j-1]+U0[i*N+j+1]+std::pow(h,2)*f(D.x0+h*(i),D.y0+h*(j)));
                err+=std::pow(U1[i*N+j]-U0[i*N+j],2);
            }
        }
        U0=U1;
    }while(std::sqrt(h*err)>toll);
    //std::cout<<iter<<std::endl;
    Uapproximate=U1;
    return U1;
}

double problem::computeError() const{
    double h=(D.x1-D.x0)/(N-1);
    double err=0;
    for(size_t i=0;i<N;++i){
        for(size_t j=0;j<N;++j){
            err+=std::pow(Uex(D.x0+i*h,D.y0+j*h)-Uapproximate[i*N+j],2);
        }
    }

    return std::sqrt(err*h);
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