#include "problem.hpp"
#include "chrono"


int main(int argc, char *argv[]){
    MPI_Init(&argc,&argv);
    if(argc!=4){
        std::cerr<<"Error:\ninput must be 3: right hand side,number of elements and number of parallel tasks"<<std::endl;
        exit(0);
    }
    if(atoi(argv[2])<0){
        std::cerr<<"Error:\nnumber of elements must be positive"<<std::endl;
        exit(0);
    }
    if(atoi(argv[3])<=0){
        std::cerr<<"Error:\nnumber of threads must be a positive integer"<<std::endl;
        exit(0);
    }


    problem P(0.0,1.0,0.0,1.0,atoi(argv[2]),argv[1]);
    std::ofstream file {"Output.vtk"};
    
    /*const auto start= std::chrono::steady_clock::now();
    P.SeqSolver(1e-7);

    const auto end= std::chrono::steady_clock::now();

    std::cout<<"L'errore è "<<P.computeError()<<" with time "<<std::chrono::duration<double>(end-start).count()<<" s"<<std::endl;*/
    //P.EraseSol();
  
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    std::vector<double> Ulocal;
    int precCell;
    
    P.EraseSol();

    const auto start= std::chrono::steady_clock::now();

    P.split(Ulocal,precCell);
    P.ParSolver(1e-7,1e4,Ulocal,precCell,atoi(argv[3]));
    P.merge(Ulocal,precCell);

    const auto end= std::chrono::steady_clock::now();
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank==0){
        std::cout<<"L' errore è "<<P.computeError()<<std::endl;
        P.printSol(file);
        std::cout<<"\nCalcolation took "<<std::chrono::duration<double>(end-start).count()<<" s"<<std::endl;
        std::cout<<"With "<<size<<" MPI thread, "<<atoi(argv[3]) <<" open MP thread and "<< atoi(argv[2])<<"x"<<atoi(argv[2])<<" mesh"<<std::endl;
    }

    MPI_Finalize();
    file.close();

    return 0;
}