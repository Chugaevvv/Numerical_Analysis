#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

class Function {
public:
    virtual double operator () (double x) = 0;
    virtual double d(double x) = 0;
};

class NewtonSolver {
protected:
    Function &F;
public:
    NewtonSolver(Function & F) : F(F) {}
    double solve(double x0) {
        // 在这里实现你的牛顿迭代法，将答案返回
        const int max_i = 5;

        double x = x0;
        for (int i = 0; i < max_i; i++) {
            double fx = F(x);
            double dfx = F.d(x);
            double x_new = x - fx / dfx;  // 牛顿迭代公式: x_{n+1} = x_n - f(x_n)/f'(x_n)
            x = x_new;
        }
        
        return x;
    }
};


// 你可以使用下面的代码在自己电脑上测试正确性，但提交时请勿包含以下代码
class Func1 : public Function {
public:
    double operator () (double x) {
        return x - tan(x);
    }
    double d(double x) {
        return 1 - 1 / pow(cos(x), 2);
    }
};

int main() {
    Func1 func;
    NewtonSolver solver(func);
    double x = solver.solve(4.5);
    cout << fixed << setprecision(10) << x << endl;

    double y = x - tan(x);
    cout << fixed << setprecision(10) << y << endl; // 输出函数值
    return 0;
}
