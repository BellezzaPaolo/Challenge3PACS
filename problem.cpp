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
    //std::cout<<iter<<std::endl;
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

void problem::split(std::vector<double>& Ulocal){
    int rank, size;
    unsigned local_N=0;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    std::vector<int> sendcount(size);
    std::vector<int> display(size);

    if(rank==0){

        int residual=0;
        if(N%size!=0){
            residual=N%size;
        }
        for(size_t i=0;i<size;++i){
            sendcount[i]=(N/size+2+(residual>0)-(i==0 || i==size-1))*N;
            residual--;
        }
        display[0]=0;
        for(size_t i=1;i<size;++i){
            display[i]=display[i-1]+sendcount[i-1]-2*N;
        }

        for(size_t i=0;i<N;++i){
            for(size_t j=0;j<N;++j){
                for(size_t k=0;k<size;++k){
                    Uapproximate[i*N+j]+=((i*N+j)>=(display[k]+sendcount[k]-N));
                }
                std::cout<<Uapproximate[i*N+j]<<" ";
            }
            std::cout<<std::endl;
        }
        /*for(size_t i=0;i<size;++i){
            std::cout<<sendcount[i]<<" ";
        }*/
        /*std::cout<<"\n"<<std::endl;
        for(size_t i=0;i<size;++i){
            std::cout<<display[i]<<" ";
        }
        std::cout<<"\n"<<std::endl;*/
    }
    
    MPI_Scatter(sendcount.data(),1,MPI_UNSIGNED,&local_N,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    //std::cout<<"rank "<<rank<<" has "<<local_N<<std::endl;

    Ulocal.resize(local_N);

    MPI_Scatterv(Uapproximate.data(),sendcount.data(),display.data(),MPI_DOUBLE,Ulocal.data(),local_N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);


/*for(size_t j=0;j<size;++j){
if(rank==j){
    std::cout<<"\nrank "<<rank<<std::endl;
    for(size_t i=0;i<local_N/N;++i){
        for(size_t j=0;j<N;++j){
            std::cout<<Ulocal[i*N+j]<<" ";
        }
        std::cout<<std::endl;
    }
MPI_Barrier(MPI_COMM_WORLD);
}
}*/
   
    return;
}

std::vector<double>& problem::ParSolver(double toll,std::vector<double>& Ulocal){
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    unsigned n_local=Ulocal.size()/N;
    std::vector<double> UlocalV(n_local,0.0);
    double err=1;
    int converged=0;
    std::vector<int> sendcount(size,0);
    std::vector<int> display(size,0);
    int recvsize=0;
    int iter=0;

    /*for(size_t j=0;j<size;++j){
if(rank==j){
    std::cout<<"\nrank "<<rank<<" loc.size()"<<Ulocal.size()<<std::endl;
    for(size_t i=0;i<Ulocal.size();++i){
        std::cout<<" "<<Ulocal[i];
    }
    std::cout<<std::endl;
    for(size_t i=0;i<n_local/N;++i){
        for(size_t j=0;j<N;++j){
            std::cout<<Ulocal[i*N+j]<<" ";
        }
        std::cout<<std::endl;
    }
MPI_Barrier(MPI_COMM_WORLD);
}
MPI_Barrier(MPI_COMM_WORLD);
}*/

    if(rank==0){
        sendcount[rank+1]=N;
        display[rank+1]=(n_local-1)*N;
        recvsize=N;
    }
    else if(rank==size){
        sendcount[rank-1]=N;
        recvsize=N;
    }
    else{
        sendcount[rank-1]=N;
        sendcount[rank+1]=N;
        display[rank+1]=(n_local-1)*N;
        recvsize=2*N;
    }
    
    do{
        err=0;
        iter++;
        /*for(size_t i=1;i<n_local-1;++i){
            for(size_t j=1;j<N-1;++j){
                Ulocal[i*N+j]=0.25*(UlocalV[(i-1)*N+j]+UlocalV[(i+1)*N+j]+UlocalV[i*N+j-1]+UlocalV[i*N+1]+std::pow(h,2)*f(D.x0+h*i,D.y0+h*j+(n_local-2)*h*rank));
                //std::cout<<D.x0+h*i<<" "<<D.y0+h*j+n_local*h*rank<<" h "<<h<<" nloc "<<n_local<<" rank "<<rank<<" j "<<j<<std::endl;
                err+=std::pow(Ulocal[i*N+j]-UlocalV[i*N+j],2);
            }
        }
        err=std::sqrt(h*err);*/
        converged=err>toll;
        MPI_Allreduce(MPI_IN_PLACE,&converged,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
        UlocalV=Ulocal;

        /*std::cout<<"rank "<<rank<<": nloc "<<n_local<<" N "<<N<<" U.size()= "<<Ulocal.size()<<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);*/


        MPI_Scatterv(UlocalV.data(),sendcount.data(),display.data(),MPI_DOUBLE,UlocalV.data(),recvsize,MPI_DOUBLE,rank,MPI_COMM_WORLD);

        //std::cout<<"rank "<<rank<<": nloc "<<n_local<<" N "<<N<<" U.size()= "<<Ulocal.size()<<"at iter= "<<iter<<std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);

        /*for(size_t j=0;j<size;++j){
            if(rank==j){
                std::cout<<"\nrank "<<rank<<std::endl;
                for(size_t i=0;i<n_local;++i){
                    for(size_t j=0;j<N;++j){
                        std::cout<<Ulocal[i*N+j]<<" ";
                    }
                    std::cout<<std::endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        std::cout<<iter<<std::endl;*/


    }while(converged);

    Uapproximate=Ulocal;

    return Ulocal;
}

void problem::EraseSol(){
    Uapproximate.shrink_to_fit();
    Uapproximate.resize(N*N);
    return;
}