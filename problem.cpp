#include "problem.hpp"

std::vector<double> problem::SeqSolver(double toll){
    std::vector<double> U1(N*N,0.0),U0(N*N,0.0);
    double err=0;
    int iter=0;

    do{
        err=0;
        iter++;
        for(size_t i=1;i<N-1;++i){
            for(size_t j=1;j<N-1;++j){
                U1[i*N+j]= 0.25*(U0[(i-1)*N+j]+U0[(i+1)*N+j]+U0[i*N+j-1]+U0[i*N+j+1]+std::pow(h,2)*f(D.x0+i*h,D.y0+j*h));//Eval(i,j,f));
                err+=std::pow(U1[i*N+j]-U0[i*N+j],2);
            }
        }
        U0=U1;
    }while(std::sqrt(h*err)>toll);
    std::cout<<iter<<std::endl;
    Uapproximate=U1;
    return U1;
}
/*
double problem::Eval(size_t i,size_t j,mu::Parser f) const{
    double x=D.x0+i*h, y=D.y0+j*h;

    f.DefineVar("x",&x);
    f.DefineVar("y",&y);

    return f.Eval();
}*/

double problem::computeError() const{
    double err=0;
    for(size_t i=0;i<N;++i){
        for(size_t j=0;j<N;++j){
            err+=std::pow(Uex(D.x0+i*h,D.y0+j*h)-Uapproximate[i*N+j],2);
        }
    }

    return std::sqrt(err*h);
}

void problem::printSol(std::ofstream& file) const{
    if(!file){
        std::cout<<"Error:\nCannot open the file"<<std::endl;
        exit(0);
    }

    file<<"# vtk DataFile Version 3.0"<<std::endl<<"Output challenge3"<<std::endl;
    file<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;
    file<<"DIMENSIONS "<<N<<" "<<N<<" 1"<<std::endl<<"POINTS "<<N*N<<" double"<<std::endl;
    
    for(size_t i=0;i<N;++i){
        for(size_t j=0;j<N;++j){
            file<<D.x0+i*h<<" "<<D.y0+j*h<<" "<<Uapproximate[i*N+j]<<std::endl;
        }
    }
    
    return;
}

void problem::split()const{
    int rank, size, residual=0;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    std::vector<int> local_n(size);

    if(rank==0){
        if(N%size!=0){
            residual=N%size;
        }
        for(size_t i=0;i<size;++i){
            
        }
    }

    return;
}

std::vector<double> problem::ParSolver(double toll){
    
}