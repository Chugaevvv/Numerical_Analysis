#include <vector>
#include <string>
#include <Eigen/Eigen>
using namespace std;

#ifndef SPLINE_WRAPPER
#define SPLINE_WRAPPER
#include <memory>
// SplineWrapper 用于将您返回的自定义 class 包装为标准接口，可以接收的返回类型包括：普通类对象、派生类对象、lambda 函数。
class SplineWrapper {
    struct Concept { virtual ~Concept() = default; virtual double eval(double) const = 0; };

    template<class T>
    struct Model final : Concept {
        T impl;
        explicit Model(T x) : impl(std::move(x)) {}
        double eval(double x) const override { return impl(x); }
    };

    std::shared_ptr<const Concept> self;

public:
    SplineWrapper() = default;

    template<class T>
    SplineWrapper(T&& x)
    : self(std::make_shared<Model<std::decay_t<T>>>(std::forward<T>(x))) {}

    double operator()(double x) const { return self->eval(x); }
    explicit operator bool() const { return static_cast<bool>(self); }
};
#endif

SplineWrapper getCubicSpline(
    const vector<double> &t, 
    const vector<double> &f, 
    const string &bdType, 
    const vector<double> &bdVal) 
{
    int n = t.size();
    vector<double> h(n - 1);
    for (int i = 0; i < n - 1; ++i)
        h[i] = t[i + 1] - t[i];

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);

    // 内部节点方程
    for (int i = 1; i < n - 1; ++i) {
        A(i, i - 1) = h[i - 1];
        A(i, i)     = 2 * (h[i - 1] + h[i]);
        A(i, i + 1) = h[i];
        rhs(i) = 6.0 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]);
    }

    if (bdType == "Natural") {
        A(0, 0) = 1.0;
        A(n - 1, n - 1) = 1.0;
        rhs(0) = rhs(n - 1) = 0.0;
    }
    else if (bdType == "Complete") {
        double fp0 = bdVal[0], fpn = bdVal[1];

        A(0, 0) = 2 * h[0];
        A(0, 1) = h[0];
        rhs(0) = 6 * ((f[1] - f[0]) / h[0] - fp0);

        A(n - 1, n - 2) = h[n - 2];
        A(n - 1, n - 1) = 2 * h[n - 2];
        rhs(n - 1) = 6 * (fpn - (f[n - 1] - f[n - 2]) / h[n - 2]);
    }
    else if (bdType == "Second-Derivatives") {
        double s0 = bdVal[0], sn = bdVal[1];
        A(0, 0) = 1.0;
        A(n - 1, n - 1) = 1.0;
        rhs(0) = s0;
        rhs(n - 1) = sn;
    }
    else if (bdType == "Periodic") {
        // 周期条件：f、f'、f'' 连续
        for (int i = 1; i < n - 1; ++i) {
            A(i, i - 1) = h[i - 1];
            A(i, i)     = 2 * (h[i - 1] + h[i]);
            A(i, i + 1) = h[i];
            rhs(i) = 6 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]);
        }

        // 周期连接
        A(0, 0) = 2 * (h[0] + h[n - 2]);
        A(0, 1) = h[0];
        A(0, n - 2) = h[n - 2];
        rhs(0) = 6 * ((f[1] - f[0]) / h[0] - (f[n - 1] - f[n - 2]) / h[n - 2]);

        A(n - 1, 0) = 1.0;
        A(n - 1, n - 1) = -1.0;
        rhs(n - 1) = 0.0;
    }

    else if (bdType == "Not-A-Knot") {
        // Not-A-Knot 条件
        // 首端
        A(0, 0) = h[1];
        A(0, 1) = -(h[0] + h[1]);
        A(0, 2) = h[0];
        rhs(0) = 0.0;
        // 尾端
        A(n - 1, n - 3) = h[n - 2];
        A(n - 1, n - 2) = -(h[n - 3] + h[n - 2]);
        A(n - 1, n - 1) = h[n - 3];
        rhs(n - 1) = 0.0;
    }

    Eigen::PartialPivLU<Eigen::MatrixXd> solver(A);
    Eigen::VectorXd c_eig = solver.solve(rhs);

    vector<double> c(n);
    for (int i = 0; i < n; ++i)
        c[i] = c_eig(i);

    vector<double> a = f;
    vector<double> b(n - 1), d(n - 1);

    for (int i = 0; i < n - 1; ++i) {
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 6.0;
        d[i] = (c[i + 1] - c[i]) / (6.0 * h[i]);
    }

    struct Cubic {
        vector<double> t, a, b, c, d;
        double operator()(double x) const {
            int n = t.size() - 1;
            if (x <= t.front()) return a.front();
            if (x >= t.back()) return a.back();
            int i = 0;
            while (i < n - 1 && x > t[i + 1]) ++i;
            double dx = x - t[i];
            return a[i] + b[i]*dx + c[i]*dx*dx/2.0 + d[i]*dx*dx*dx;
        }
    };

    Cubic spline{t, a, b, c, d};
    return spline;
}