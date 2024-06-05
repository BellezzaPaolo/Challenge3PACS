#include "problem.hpp"

std::vector<double> problem::SeqSolver(double toll){
    std::vector<double> U1(N*N,0.0),U0(N*N,0.0);
    double err=0;
    int iter=0;
    double x=0,y=0;

    f.DefineVar("x",&x);
    f.DefineVar("y",&y);

    do{
        err=0;
        iter++;
        for(size_t i=1;i<N-1;++i){
            for(size_t j=1;j<N-1;++j){
                x=D.x0+i*h;
                y=D.y0+j*h;
                U1[i*N+j]= 0.25*(U0[(i-1)*N+j]+U0[(i+1)*N+j]+U0[i*N+j-1]+U0[i*N+j+1]+std::pow(h,2)*f.Eval());
                err+=std::pow(U1[i*N+j]-U0[i*N+j],2);
            }
        }
        std::swap(U0,U1);
    }while(std::sqrt(h*err)>toll);
    //std::cout<<iter<<std::endl;
    Uapproximate=U1;
    return U1;
}


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

void problem::split(std::vector<double>& Ulocal,int& precCell){
    int rank, size;
    unsigned local_N=0;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    std::vector<int> sendcount(size);
    std::vector<int> display(size,0);

    if(rank==0){

        int residual=0;
        if(N%size!=0){
            residual=N%size;
        }

        #pragma omp parallel for num_threads(size)
            for(size_t i=0;i<size;++i){
                sendcount[i]=(N/size+2+(i<residual)-(i==0 || i==size-1)-(i==0 && i==size-1))*N;
            }
        
        for(size_t i=1;i<size;++i){
            display[i]=display[i-1]+sendcount[i-1]-2*N+N*(rank==size-1);
        }


        /*for(size_t i=0;i<N;++i){
            std::cout<<i<<"  | ";
            for(size_t j=0;j<N;++j){
                Uapproximate[i*N+j]=0;
                Uapproximate[5*N+1]=1;
                Uapproximate[15*N+2]=1;
                std::cout<<Uapproximate[i*N+j]<<" ";
            }
            std::cout<<std::endl;
        }*/
        /*for(size_t i=0;i<size;++i){
            std::cout<<sendcount[i]<<" ";
        }*/
        //std::cout<<"\n"<<std::endl;
        /*for(size_t i=0;i<size;++i){
            std::cout<<display[i]<<" ";
        }
        std::cout<<"\n"<<std::endl;*/
    }
    
    MPI_Scatter(sendcount.data(),1,MPI_UNSIGNED,&local_N,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    MPI_Scatter(display.data(),1,MPI_INT,&precCell,1,MPI_INT,0,MPI_COMM_WORLD);

    Ulocal.resize(local_N);

    MPI_Scatterv(Uapproximate.data(),sendcount.data(),display.data(),MPI_DOUBLE,Ulocal.data(),local_N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);


/*for(size_t j=0;j<size;++j){
    if(rank==j){
        std::cout<<"\nsplit rank "<<rank<<std::endl;
        for(size_t i=0;i<local_N/N;++i){
            std::cout<<precCell/N+i<<"  | ";
            for(size_t j=0;j<N;++j){
                std::cout<<Ulocal[i*N+j]<<" ";
            }
            std::cout<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}*/
   
    return;
}

std::vector<double>& problem::ParSolver(double toll,int MaxIter,std::vector<double>& Ulocal, int& precCell,int n_threads){
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    unsigned n_local=Ulocal.size()/N;
    std::vector<double> UlocalV(n_local*N,0.0);
    std::vector<double> rhs(n_local*N,0);
    double err=0;
    int converged=0;
    int iter=0;
    double x=D.x0, y=D.y0;

    UlocalV=Ulocal;
    f.DefineVar("x",&x);
    f.DefineVar("y",&y);

    for(size_t i=1;i<n_local;++i){
        for(size_t j=1;j<N;++j){
            x=D.x0+(precCell/N+i)*h;
            y=D.y0+j*h;
            rhs[i*N+j]=f.Eval();
        }
    }



    //#pragma omp parallel num_threads(n_threads)
    {
        do{
            err=0;
            
            //#pragma omp parallel for 
                for(size_t i=1;i<n_local-1;++i){
                    //#pragma omp parallel for
                        for(size_t j=1;j<N-1;++j){
                            Ulocal[i*N+j]=0.25*(UlocalV[(i-1)*N+j]+UlocalV[(i+1)*N+j]+UlocalV[i*N+j-1]+UlocalV[i*N+j+1]+std::pow(h,2)*rhs[i*N+j]);
                            err+=std::pow(Ulocal[i*N+j]-UlocalV[i*N+j],2);
                        }
                }


            //#pragma omp barrier

            //#pragma omp single
            {
            err=std::sqrt(h*err);
            
            converged=(err>toll && iter<MaxIter);
            MPI_Allreduce(MPI_IN_PLACE,&converged,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
            

            if(rank!=0){
                MPI_Send(&Ulocal[N],N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD);
            }
            if(rank!=size-1){
                MPI_Recv(&Ulocal[(n_local-1)*N],N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            if(rank!=size-1){
                MPI_Send(&Ulocal[(n_local-2)*N],N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
            }
            if(rank!=0){
                MPI_Recv(Ulocal.data(),N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            std::swap(UlocalV,Ulocal);
            iter++;

            if(rank==0/* && iter%100==0*/){
                std::cout<<iter<<" ";
            }

            }

        }while(converged);
}
    if(rank==0){
        std::cout<<"\n"<<iter<<std::endl;
    }
    
    return Ulocal;
}

void problem::merge(std::vector<double>& Ulocal,int& precCell){
    int size,rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> recvcount(size,0);
    std::vector<int> display(size,0);
    int n_local=Ulocal.size()/N;
    int recv=Ulocal.size();

    /*for(size_t j=0;j<size;++j){
        if(rank==j){
            std::cout<<"\nmerge rank "<<rank<<std::endl;
            for(size_t i=0;i<n_local;++i){
                std::cout<<precCell/N+i<<"  | ";
                for(size_t j=0;j<N;++j){
                    std::cout<<Ulocal[i*N+j]<<" ";
                }
                std::cout<<std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    std::cout<<"\n\n"<<std::endl;*/


    MPI_Gather(&precCell,1,MPI_INT,display.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    /*if(rank==0){
        for(size_t i=0;i<display.size();++i){
            std::cout<<display[i]<<" ";
        }
        std::cout<<std::endl;
    }*/
    MPI_Gather(&recv,1,MPI_INT,recvcount.data(),1,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Gatherv(Ulocal.data(),Ulocal.size(),MPI_DOUBLE,Uapproximate.data(),recvcount.data(),display.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);

    /*if(rank==0){
        std::cout<<std::endl;
        for(size_t i=0;i<N;++i){
            std::cout<<precCell/N+i<<"  | ";
            for(size_t j=0;j<N;++j){
                std::cout<<Uapproximate[i*N+j]<<" ";
            }
            std::cout<<std::endl;
        }
    }*/

    return;
}

void problem::EraseSol(){
    Uapproximate.shrink_to_fit();
    Uapproximate.resize(N*N);
    return;
}