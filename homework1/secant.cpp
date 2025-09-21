#include <cmath>
#include <iostream>
using namespace std;

class Function {
public:
    virtual double operator () (double x) = 0;
};

class SecantSolver {
protected:
    Function &F;
public:
    SecantSolver(Function & F) : F(F) {}
    double solve(double x0, double x1) {
        // 在这里实现你的割线法，将答案返回
        const double tolerance = 1e-12; 
        const int max_i = 30; 
        int i = 0;
        double f0 = F(x0); i++;
        double f1 = F(x1); i++;
        while (i < max_i) {
            if (abs(f1) < tolerance) {
                return x1; 
            }

            double x2 = x1 - f1 * (x1 - x0) / (f1 - f0); 
            x0 = x1;
            f0 = f1;
            x1 = x2;
            f1 = F(x1);
        }
        return x1;         
    }
};