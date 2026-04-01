#ifndef DECLARE_FUNCTIONS_HH
#define DECLARE_FUNCTIONS_HH

#include <vector>


double gaussian(const double &x);

void Thomas_tridiagonal_solver(const std::vector<double> &a,
                                     std::vector<double> &b,
                               const std::vector<double> &c,
                                     std::vector<double> &d,
                                     std::vector<double> &u);


#endif  // DECLARE_FUNCTIONS_HH
