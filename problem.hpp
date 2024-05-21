#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include<vector>
#include <functional>
#include <cmath>
#include <iostream>

using funzione=std::function<double(double, double)>;

struct domain{
    double x0,x1,y0,y1;
    domain(double X0,double X1,double Y0,double Y1):x0(X0),x1(X1),y0(Y0),y1(Y1){};
};

struct BoundaryCondition{
    funzione Ovest;
    funzione Nord;
    funzione Est;
    funzione Sud;
    BoundaryCondition(){};
};

class problem{
    private:
        domain D;
        //BoundaryCondition Bc;
        int N;
        funzione f=[](double x, double y){return 8*std::pow(M_PI,2)*sin(2*M_PI*x)*sin(2*M_PI*y);};
        funzione Uex=[](double x,double y){return sin(2*M_PI*x)*sin(2*M_PI*y);};
        std::vector<double> Uapproximate;


    public:
        problem(double x0, double x1,double y0,double y1,int n): D(x0,x1,y0,y1), N(n){};

        std::vector<double> SeqSolver(double toll);

        double computeError() const;
        
        void printSol() const;
};



#endif