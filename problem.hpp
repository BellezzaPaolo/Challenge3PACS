#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <iostream>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <mpi.h>
#include <omp.h>
#include <iomanip>
#include <muParser.h>

using funzione=std::function<double(double, double)>;

//definition of the struct domain
struct domain{
    double x0,x1,y0,y1; // represent the corners of my rectangular domain
    domain(double X0,double X1,double Y0,double Y1):x0(X0),x1(X1),y0(Y0),y1(Y1){}; 
};

struct BoundaryCondition{
    /*
    * struct with a vector of string (type) that each element indicates which type of boundary condition is applied:
    * -Dirichlet (d)
    * -Neumann (n)
    * -Robin (r)
    * if there is a sigle type along all the boundary is sufficient a vector of a single element.
    * the function mu::Parser indicates the value of the condition (the function could be defined piecewise).
    * This code implementation supports only Dirichelet boundary condition actually. 
    */
    mu::Parser condition;
    std::vector<std::string> type;
    //constructor
    BoundaryCondition(std::string bc,std::vector<std::string> c): type(c){
        condition.SetExpr(bc);
    }
};

class problem{
    private:
        domain D; //domain
        BoundaryCondition Bc; //boundary condition
        int N;  //number of cell in which is dividede every row and column
        double h=(D.x1-D.x0)/(N-1); //dimension of every cell
        mu::Parser f; //right hand side of my problem
        mu::Parser Uex; //true solution
        std::vector<double> Uapproximate; //solution approsimated

        void split(std::vector<double>& Ulocal,int& precCell); //divides the matrix in parts and assign them to avery MPI rank

        int HibridSolver(const double toll,const int MaxIter,std::vector<double>& Ulocal,int& precCell,const int n_threads); //the real parallel solver

        void merge(std::vector<double>& Ulocal,int& precCell); //merges all the different part of the solution in the final one

    public:
        //constructor
        problem(double x0, double x1,double y0,double y1,int n,std::string bc,std::vector<std::string> c,std::string F,std::string uesat): D(x0,x1,y0,y1), N(n),Uapproximate(N*N,0.0), Bc(bc,c){
            f.SetExpr(F);
            f.DefineConst("pi",M_PI);
            Uex.SetExpr(uesat);
            Uex.DefineConst("pi",M_PI);
            };

        int SeqSolver(const double toll,const int MaxIter); //solves the problem sequentially

        int ParSolver(const double toll,const int MaxIter,const int n_threads); //solves the problem in parallel using mpi and open MP

        void evaluate(mu::Parser fun,std::vector<double>& funEval,const int& precCell) const; //evaluate the muParser given in input in every point and saves the values in the vector

        double computeError() const; //computes the error in L2 norm between the true solution and the approximated one
        
        void printSol(std::ofstream& file) const; //prints the solution in a file that can be visualized

        void printPerformance(const int size,const int n_thread,const int iter,const double time) const; //prints the performance of the code at the terminal

        void EraseSol(); //cancel of the approximate solution

        void assignBoundaryCondition(); //asigns the boundary conditions
};



#endif