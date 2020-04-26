#include <iostream>
#include <fstream>
#include "matrix.cpp"
#include "simplex.h"
#include "simplex.cpp"

#if 1
    #define pv(x) std::cerr<<#x<<" = "<<(x)<<"; ";std::cerr.flush()
    #define pn std::cerr<<std::endl
#else
    #define pv(x)
    #define pn
#endif

using namespace std;


// name=tema2; g++ $name.cpp -o $name.exe && ./$name.exe



void testSimplexWithBasis() {
    ifstream fin("simplexWithBasis.in");
    Matrix<long double> A = Matrix<long double>::readMatrixFromStream(fin);
    int M = A.get1Dim();
    int N = A.get2Dim();
    pv(M);pv(N);pn;

    vector<long double> b(M);
    for (long double& val : b) {
        fin >> val;
    }

    string type;
    TypeOfObjective obj;
    fin >> type;
    assert(type == "min" || type == "max");
    if (type == "min") {
        obj = TypeOfObjective::MIN;
    }
    else {
        obj = TypeOfObjective::MAX;
    }

    vector<long double> c(N);
    for (long double& val : c) {
        fin >> val;
    }

    vector<int> basis(M);
    for (int& val : basis) {
        fin >> val;
        val -= 1;
    }

    SimplexReturnType result = runSimplex(A, Matrix<long double>(b).getTranspose(), Matrix<long double>(c), obj, basis);
    pv((int)result.result);pn;
    pv(result.state.z);pn;

    for (int i = 0; i < result.vvb.size(); ++i) {
        pv(result.vvb[i]);pn;
    }
    for (int i = 0; i < result.state.basis.size(); ++i) {
        pv(result.state.basis[i]);pn;
    }

}


void testSimplexWithoutBasis() {
    ifstream fin("simplexWithoutBasis.in");
    Matrix<long double> A = Matrix<long double>::readMatrixFromStream(fin);
    int M = A.get1Dim();
    int N = A.get2Dim();

    vector<long double> b(M);
    for (long double& val : b) {
        fin >> val;
    }

    string type;
    TypeOfObjective obj;
    fin >> type;
    assert(type == "min" || type == "max");
    if (type == "min") {
        obj = TypeOfObjective::MIN;
    }
    else {
        obj = TypeOfObjective::MAX;
    }

    vector<long double> c(N);
    for (long double& val : c) {
        fin >> val;
    }

    SimplexReturnType result = runSimplex(A, Matrix<long double>(b).getTranspose(), Matrix<long double>(c), obj);

    pv((int)result.result);pn;
    pv(result.state.z);pn;
    for (int i = 0; i < result.vvb.size(); ++i) {
        pv(result.vvb[i]);pn;
    }
    for (int i = 0; i < result.state.basis.size(); ++i) {
        pv(result.state.basis[i]);pn;
    }
}


int main() {
    testSimplexWithoutBasis();

    return 0;
}