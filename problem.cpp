#include "problem.hpp"

int problem::SeqSolver(const double toll,const int MaxIter){
    /*
    * this function solves the differential problem sequentially using the jacobi algorithm combined with the finite difference method
    * the input are:
    * -toll that is the tollerance on the difference in norm L2 to stop the iterations
    * -MaxIter is the maximum number of iterations that the method can make
    * the output is:
    * -iter that is the number of iterations made by the method to reach convergence 
    */

    std::vector<double> Un(N*N,0.0),rhs(N*N); //respectively the new solution Un and the matrix for the right hand side
    double err=0; //the error between 2 iterate
    int iter=0; //the number of iterate
    double x=0,y=0; //store the values of the abscissa (x) and ordinata (y) 
    int precCell=0; //stores the offset of that specific rank (in the sequential method is useless but necessary to make the function evaluate work)

    evaluate(f,rhs,precCell); //valuates f in every point of my domain and saves every value in rhs

    //start of the real method
    do{
        err=0; //update of the error and iterate
        iter++;
        
        for(size_t i=1;i<N-1;++i){
            for(size_t j=1;j<N-1;++j){
                Un[i*N+j]= 0.25*(Uapproximate[(i-1)*N+j]+Uapproximate[(i+1)*N+j]+Uapproximate[i*N+j-1]+Uapproximate[i*N+j+1]+std::pow(h,2)*rhs[i*N+j]);//jacobi rule of update of the matrix
                err+=std::pow(Un[i*N+j]-Uapproximate[i*N+j],2); //computation of the error
            }
        }

        std::swap(Un,Uapproximate); //switch the new solution with the old one
    }while(std::sqrt(h*err)>toll && iter<MaxIter); //check of convergence

    Uapproximate=Un; //save of the solution in the class
    
    return iter;
}

int problem::ParSolver(const double toll,const int MaxIter,const int n_threads){
    /*
    * this function solves the differential problem in parallel using the jacobi algorithm combined with the finite difference method 
    * and MPI and/or open MP parallelisation.
    * the input are:
    * -toll that is the tollerance on the difference in norm L2 to stop the iterations
    * -MaxIter is the maximum number of iterations that the method can make
    * -n_threads is the number of thread of open MP for every rank of MPI
    * the output is:
    * -iter that is the number of iterations made by the method to reach convergence 
    */
    std::vector<double> Ulocal; //the local part of the solution relative to that rank
    int precCell; //number of cells given to previuos ranks

    split(Ulocal,precCell); //divides the matrix in parts and assign them to avery MPI rank
    int iter=HibridSolver(toll,MaxIter,Ulocal,precCell,n_threads); //the real parallel solver
    merge(Ulocal,precCell); //merges all the different part of the solution in the final one

    return iter;
}

void problem::evaluate(mu::Parser fun,std::vector<double>& funEval,const int& precCell)const{
    /*
    * evaluate the muParser given in input in every point and saves the values in the vector.
    * The inputs are:
    * -fun is the function that needs to be evaluated
    * -funEval is the vector that stores the value of fun in every point
    * -precCell is the number of cells given to previuos ranks and so the offeset for the coordinates
    */
    double x=D.x0, y=D.y0; //the abscissa (x) and the ordiante (y)
    size_t n_row=funEval.size()/N; //number of rows to evaluate the function

    fun.DefineVar("x",&x); //definition of the 2 variables of my function
    fun.DefineVar("y",&y);

    for(size_t i=0;i<n_row;++i){
        for(size_t j=0;j<N;++j){
            x=D.x0+(i+precCell)*h; //calculation of the coordinate realted to the (i+precCell)-th row
            y=D.y0+j*h; //calculation of the coordinate realted to the j-th column
            funEval[i*N+j]=fun.Eval(); //evaluation
        }
    }

    return;
}

double problem::computeError() const{
    /*
    * Computes the error in L2 norm between the true solution and the approximated one.
    * The output is:
    * -the error in L2 norm
    */
    double err=0;
    std::vector<double> Ues(N*N); //stores the evaluetion of the exact solution
    int precCell=0;//stores the offset of that specific rank (in this case is useless but necessary to make the function evaluate work)

    evaluate(Uex,Ues,precCell);

    //compute the L2 norm of the difference in every point
    #pragma omp parallel for
        for(size_t i=0;i<N;++i){
            for(size_t j=0;j<N;++j){
                err+=std::pow(Ues[i*N+j]-Uapproximate[i*N+j],2);
            }
        }

    return std::sqrt(err*h);
}

void problem::printSol(std::ofstream& file) const{
    /*
    * Prints the solution in a file that can be visualized.
    * The input is:
    * - file is the name of the vtk file where print the solution
    */

    //check if there is an error in input
    if(!file){
        std::cout<<"Error:\nCannot open the file"<<std::endl;
        exit(0);
    }

    //first lines of the vtk file needed to comunicate with ParaView
    file<<"# vtk DataFile Version 3.0"<<std::endl<<"Output challenge3"<<std::endl;
    file<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;
    file<<"DIMENSIONS "<<N<<" "<<N<<" 1"<<std::endl<<"POINTS "<<N*N<<" double"<<std::endl;
    
    //print of every point in the form: x y value
    for(size_t i=0;i<N;++i){
        for(size_t j=0;j<N;++j){
            file<<D.x0+i*h<<" "<<D.y0+j*h<<" "<<Uapproximate[i*N+j]<<std::endl;
        }
    }
    
    return;
}

void problem::split(std::vector<double>& Ulocal,int& precCell){
    /*
    * Divides the matrix in parts and assign them to avery MPI rank.
    * The inputs are:
    * -Ulocal: will store the part of my solution that this rank will compute
    * -precCell: will store the number of cells given to previous ranks
    */

    int rank, size;// respectively sotre the number of the rank and the total number of MPI ranks
    unsigned local_N=0; //will store the local number of elements
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    std::vector<int> sendcount(size);//will store the number of element that will be send to the i-th rank in the i-th entry
    std::vector<int> display(size,0); //will store the displaicement that will be used to send elements to the other ranks

    if(rank==0){
        //since the communications starts only from rank 0 it's the only one that initialise all the vectors
        int residual=0;
        if(N%size!=0){// to make the number of row to avery rank the most equal as possible I use this variable residual
            residual=N%size;
        }

        //this 2 for loop initialize the 2 vector properly
        #pragma omp parallel for num_threads(size)
            for(size_t i=0;i<size;++i){
                sendcount[i]=(N/size+2+(i<residual)-(i==0 || i==size-1)-(i==0 && i==size-1))*N;
            }
        
        for(size_t i=1;i<size;++i){
            display[i]=display[i-1]+sendcount[i-1]-2*N+N*(rank==size-1);
        }
    }
    
    //send all the data to the right MPI rank
    MPI_Scatter(sendcount.data(),1,MPI_UNSIGNED,&local_N,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    MPI_Scatter(display.data(),1,MPI_INT,&precCell,1,MPI_INT,0,MPI_COMM_WORLD);

    Ulocal.resize(local_N); //prepare the size of the local solution

    //send of the local solution
    MPI_Scatterv(Uapproximate.data(),sendcount.data(),display.data(),MPI_DOUBLE,Ulocal.data(),local_N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
    return;
}

int problem::HibridSolver(const double toll,const int MaxIter,std::vector<double>& Ulocal, int& precCell,const int n_threads){
    /*
    * The real parallel solver.
    * The inputs are:
    * -toll: that is the tollerance on the difference in norm L2 to stop the iterations
    * -MaxIter: is the maximum number of iterations that the method can make
    * -Ulocal: part of the solution that this rank must compute
    * -precCell: number of cells that preceeds this rank and so is the offeset of row in the calculation
    * -n_thread: number of open MP thread
    * the output is:
    * -iter that is the number of iterations made by the method to reach convergence 
    */

    int rank,size; //respectively the number of this rank and the total number of MPI rank
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    unsigned n_local=Ulocal.size()/N; //number of rows that this rank must compute
    std::vector<double> UlocalV=Ulocal; //previous iterate
    std::vector<double> rhs(n_local*N,0); //vector that will store the evaluation of the right hand side
    double err=0; //will store the error
    int converged=0; //boolean variable that will check the convergence
    int iter=0; //number of iterations

    evaluate(f,rhs,precCell); //evaluates the function f in every point of my domain


    //begin of the iterative solver
    #pragma omp parallel num_threads(n_threads)
    {
        do{
            err=0;
            
            //compute the new iterate and with that the L2 norm of the difference beetween the previous one
            #pragma omp parallel for 
                for(size_t i=1;i<n_local-1;++i){
                    //#pragma omp parallel for
                        for(size_t j=1;j<N-1;++j){
                            Ulocal[i*N+j]=0.25*(UlocalV[(i-1)*N+j]+UlocalV[(i+1)*N+j]+UlocalV[i*N+j-1]+UlocalV[i*N+j+1]+std::pow(h,2)*rhs[i*N+j]);
                            err+=std::pow(Ulocal[i*N+j]-UlocalV[i*N+j],2);
                        }
                }


            #pragma omp barrier

            #pragma omp single
            {
            err=std::sqrt(h*err);
            
            converged=(err>toll && iter<MaxIter); //checks the local convergence
            MPI_Allreduce(MPI_IN_PLACE,&converged,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD); //check the global convergence
            
            //communication of the first computated row to the previous rank
            if(rank!=0){
                MPI_Send(&Ulocal[N],N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD);
            }
            if(rank!=size-1){
                MPI_Recv(&Ulocal[(n_local-1)*N],N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            //communication of the last computed row to the next rank
            if(rank!=size-1){
                MPI_Send(&Ulocal[(n_local-2)*N],N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
            }
            if(rank!=0){
                MPI_Recv(Ulocal.data(),N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }

            //update the iteration and the solution
            std::swap(UlocalV,Ulocal);
            iter++;

            }

        }while(converged);
}

    return iter;
}

void problem::merge(std::vector<double>& Ulocal,int& precCell){
    /*
    * Merges all the different part of the solution in the final one.
    * The inputs are:
    * -Ulocal: part of the solution computed by that rank
    * -precCell: number of cells that preceeds that rank and so is the offeset of row in the merging
    */

    int size,rank; // number of MPI total rank and number of that MPI rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> recvcount(size,0);//will store the number of element that the rank will receive form the i-th rank in the i-th cell
    std::vector<int> display(size,0); //will store the offset to apply to every communication received
    int n_local=Ulocal.size()/N; //number of rows that the local solution has
    int recv=Ulocal.size(); //numebr of cells that the local solution has


    MPI_Gather(&precCell,1,MPI_INT,display.data(),1,MPI_INT,0,MPI_COMM_WORLD); //initialise the display vector

    MPI_Gather(&recv,1,MPI_INT,recvcount.data(),1,MPI_INT,0,MPI_COMM_WORLD); //initialize the recvcount vector

    MPI_Gatherv(Ulocal.data(),Ulocal.size(),MPI_DOUBLE,Uapproximate.data(),recvcount.data(),display.data(),MPI_DOUBLE,0,MPI_COMM_WORLD); //merges all the local solutions

    return;
}

void problem::printPerformance(const int size,const int n_thread,const int iter,const double time) const{
    /*
    * Prints the performance of the code at the terminal.
    * The inputs are:
    * -size: number of ranks of MPI
    * -n_thread: number of open MP thread
    * -iter: number of itereation made by my algorithm to compute the solution
    * -time: seconds usaed by the code to compute the solution
    */
   
    //prints everything
    std::cout<<"With "<<size<<" MPI thread, "<<n_thread<<" open MP thread and "<< N<<"x"<<N<<" mesh"<<std::endl;
    std::cout<<"The error is "<<computeError()<<" and to reach convergence took "<<iter<<" iterations"<<std::endl;
    std::cout<<"\nCalcolation needed "<<time<<" s"<<std::endl;

    return;
}

void problem::EraseSol(){
    /*
    * sets to 0 the approximate solution ready to make a new computation
    */
    Uapproximate.shrink_to_fit();
    Uapproximate.resize(N*N);
    return;
}

void problem::assignBoundaryCondition(){
    /*
    * Assigns the boundary conditions
    */
    double x,y;
    Bc.condition.DefineVar("x",&x);
    Bc.condition.DefineVar("y",&y);
    if(Bc.type.size()==1){ //circle on all the boundary
        for(size_t i=0;i<N;++i){
            x=D.x0+i*h;
            y=D.y0;
            Uapproximate[i*N+0]=Bc.condition.Eval();
            y=D.y1;
            Uapproximate[i*N+N-1]=Bc.condition.Eval();
        }
        for(size_t j=0;j<N;++j){
            y=D.y0+j*h;
            x=D.x0;
            Uapproximate[0*N+j]=Bc.condition.Eval();
            x=D.x1;
            Uapproximate[(N-1)*N+j]=Bc.condition.Eval();
        }
    }
    else{//not done yet
        for(size_t i=0;i<Bc.type.size();++i){
            std::cerr<<"\nError: not done yet, dupports only Dirichlet condition along all the boundary so a vector of one element is needed";
            exit(0);
        }
    }
    return;
}