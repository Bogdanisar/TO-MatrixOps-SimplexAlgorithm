#ifndef SIMPLEX_HEADER_GUARD
#define SIMPLEX_HEADER_GUARD


#include <utility>
#include "matrix.cpp"

enum class SimplexResult {
    NO_SOLUTION, OPTIMAL_SOLUTION, INFINITE_SOLUTION
};

enum class TypeOfObjective {
    MIN, MAX
};

struct SimplexState {
    Matrix<long double> BbA;
    Matrix<long double> z_c;
    long double z;
    vector<int> basis;
};

struct SimplexReturnType {
    SimplexResult result;
    vector<long double> vvb;
    SimplexState state;
};

// run the simplex algorithm without specifying a starting basis
SimplexReturnType runSimplex(
    Matrix<long double> A,
    Matrix<long double> b,
    Matrix<long double> c,
    TypeOfObjective obj
);

// run the simplex algorithm with a specified starting basis
SimplexReturnType runSimplex(
    Matrix<long double> A,
    Matrix<long double> b,
    Matrix<long double> c,
    TypeOfObjective obj,
    vector<int> basis
);




// run the dual simplex algorithm without specifying a starting basis
SimplexReturnType runDualSimplex(
    Matrix<long double> A,
    Matrix<long double> b,
    Matrix<long double> c,
    TypeOfObjective obj,
    vector<int> basis
);

// run the dual simplex algorithm with a specified starting basis
SimplexReturnType runDualSimplex(
    Matrix<long double> A,
    Matrix<long double> b,
    Matrix<long double> c,
    TypeOfObjective obj,
    vector<int> basis
);



#endif // SIMPLEX_HEADER_GUARD