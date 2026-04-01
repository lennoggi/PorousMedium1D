#include <cassert>
#include <vector>

using namespace std;


void Thomas_tridiagonal_solver(const vector<double> &a,    // Sub-diagonal elements (untouched)
                                     vector<double> &b,    // Diagonal elements
                               const vector<double> &c,    // Above-diagonal elements (untouched)
                                     vector<double> &d,    // RHS elements
                                     vector<double> &x) {  // Solution vector
    const auto N = x.size();
    assert(N > 1);
    const auto Nm1 = N - 1;

    assert(a.size() == Nm1);
    assert(b.size() == N);
    assert(c.size() == Nm1);

    // Forward sweep
    for (auto i = decltype(N){1}; i < N; ++i) {
        const auto im1 = i - 1;
        const auto m = a[im1]/b[im1];
        b[i] -= m*c[im1];
        d[i] -= m*d[im1];
    }

    // Backward solve
    x[Nm1] = d[Nm1]/b[Nm1];
    const auto Nm2 = N - 2;

    for (auto i = decltype(N){Nm2}; i > 0; --i) {
        x[i] = (d[i] - c[i]*x[i+1]) / b[i];
    }

    return;
}
