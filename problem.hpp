#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include<vector>
#include <functional>

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
        funzione f=0;
        funzione Uex=0;
        std::vector<double> Uapproximate;


    public:
        problem(double x0, double x1,double y0,double y1): D(x0,x1,y0,y1){};

        std::vector<double> SeqSolver(int N);
        
        void printSol() const;
};



#endif