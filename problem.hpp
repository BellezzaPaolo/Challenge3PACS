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
//#include <muParser.h>
//#include <muParserDef.h>

using funzione=std::function<double(double, double)>;

struct domain{
    double x0,x1,y0,y1;
    domain(double X0,double X1,double Y0,double Y1):x0(X0),x1(X1),y0(Y0),y1(Y1){};
};

enum class BoundCond {Dirichlet,Neumann,Robin};

/*template <BoundCond B>
struct BoundaryCondition{
    std::array<BoundCond,4> tipe;
    mu::Parser Top;
    mu::Parser Right;
    mu::Parser Bottom;
    mu::Parser Left;
    BoundaryCondition(){};
};*/

class problem{
    private:
        domain D;
        //BoundaryCondition Bc;
        int N;
        double h=(D.x1-D.x0)/(N-1);
        //mu::Parser f;
        funzione f=[](double x, double y){return 8*std::pow(M_PI,2)*sin(2*M_PI*x)*sin(2*M_PI*y);};
        funzione Uex=[](double x,double y){return sin(2*M_PI*x)*sin(2*M_PI*y);};
        std::vector<double> Uapproximate;


    public:
        problem(double x0, double x1,double y0,double y1,int n,std::string F): D(x0,x1,y0,y1), N(n),Uapproximate(N*N,0.0){};//f.SetExpr(F);};

        std::vector<double> SeqSolver(double toll);

        /*double Eval(size_t i,size_t j,mu::Parser f) const;*/

        double computeError() const;
        
        void printSol(std::ofstream& file) const;

        void split(std::vector<double>& Ulocal);

        std::vector<double>& ParSolver(double toll,std::vector<double>& Ulocal);

        void EraseSol();
};



#endif