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
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    std::vector<double> Ulocal;
    int precCell;
/*{
    std::vector<int> send, sendcounts(size,0),display(size,0),recv(20,2);
    if(rank==0){
        send.resize(10);
        for(size_t i=0;i<send.size();++i){
            send[i]=i;
        }
        for(size_t i=1;i<size;++i){
            sendcounts[i]=10-i;
            display[i]=i;
        }
    }
    MPI_Scatterv(send.data(),sendcounts.data(),display.data(),MPI_INT,&recv[2],10,MPI_INT,0,MPI_COMM_WORLD);
    for(size_t i=0;i<size;++i){
        if(rank==i){
            std::cout<<"rank "<<rank<<std::endl;
            for(size_t j=0;j<20;++j){
                std::cout<<recv[j]<<" ";
            }
            std::cout<<"\n"<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}*/
    
    P.EraseSol();
    P.split(Ulocal,precCell);
    P.ParSolver(1e-7,Ulocal,precCell);
    P.merge(Ulocal,precCell);


        
    if(rank==0){
        std::cout<<"L' errore è "<<P.computeError()<<std::endl;
        P.printSol(file);
    }

    MPI_Finalize();
    file.close();

    return 0;
}