#include <cmath>
#include <iostream>
using namespace std;

class Function {
public:
    virtual double operator () (double x) = 0;
};

class BisectionSolver {
protected:
    Function &F;
public:
    BisectionSolver(Function & F) : F(F) {}
    double solve(double l, double r) {
        const double tolerance = 1e-12;
        const int max_i = 40;
        int i = 0;
        double fl = F(l); i++;
        double fr = F(r); i++;
        if (fl * fr > 0) {
            return l;
        }
        while (i < max_i) {
            double mid = (l + r) / 2.0;
            double fmid = F(mid); i++;
            if (abs(fmid) < tolerance) {
                return mid;
            }
            if (fl * fmid < 0) {
                r = mid;
                fr = fmid;
            } else {
                l = mid;
                fl = fmid;
            }
        }
        return (l + r) / 2.0;
    }
};