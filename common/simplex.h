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

// afiseaza rezultatul executiei unui algoritm simplex;
void printSimplexReturnType(SimplexReturnType ret);


// executa algoritmul simplex (pe o functie obiectiv de minim sau maxim) fara a specifica o baza admisibila initiala
SimplexReturnType runSimplex(
    Matrix<long double> A,
    Matrix<long double> b,
    Matrix<long double> c,
    TypeOfObjective obj
);


// executa algoritmul simplex primal pe o baza admisibila specificata (si o functie obiectiv de minim sau maxim)
SimplexReturnType runSimplex(
    Matrix<long double> A,
    Matrix<long double> b,
    Matrix<long double> c,
    TypeOfObjective obj,
    vector<int> basis
);





// executa algoritmul simplex dual cu o baza dual admisibila specificata si cu o functie obiectiv de minim / maxim
SimplexReturnType runDualSimplexWithDualBasis(
    Matrix<long double> A,
    Matrix<long double> b,
    Matrix<long double> c,
    TypeOfObjective obj,
    vector<int> basis
);

// executa algoritmul simplex dual pe o baza oarecare si cu o functie obiectiv de minim sau de maxim
SimplexReturnType runDualSimplexWithAnyBasis(
    Matrix<long double> A,
    Matrix<long double> b,
    Matrix<long double> c,
    TypeOfObjective obj,
    vector<int> basis
);

#endif // SIMPLEX_HEADER_GUARD