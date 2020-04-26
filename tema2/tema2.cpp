#include <iostream>
#include <fstream>
#include "../common/matrix.cpp"

#if 1
    #define pv(x) std::cerr<<#x<<" = "<<(x)<<"; ";std::cerr.flush()
    #define pn std::cerr<<std::endl
#else
    #define pv(x)
    #define pn
#endif

using namespace std;


// name=tema2; g++ $name.cpp -o $name.exe && ./$name.exe


using Numeric = long double;

void test1() {
    ifstream fin("matrix.in");
    Matrix<Numeric> matrix = Matrix<Numeric>::readMatrixFromStream(fin);
    Matrix<long double> m(vector<vector<long double>>());

    pv(matrix);pn;
    pv(matrix.getTranspose());pn;
    pv((-5) * matrix);pn;
    pv(matrix + (-1) * matrix);pn;
    pv(matrix.getSubmatrixWithoutRowAndColumn(1, 2));pn;
}

void testSquare() {
    ifstream fin("lemma.in");
    Matrix<Numeric> matrix = Matrix<Numeric>::readMatrixFromStream(fin);
    int index;
    fin >> index;

    int N = matrix.get1Dim();
    vector<long double> column(N);
    for (int i = 0; i < N; ++i) {
        fin >> column[i];
    }



    pv(matrix.getDeterminant());pn;
    pv(matrix.getCofactorMatrix());pn;
    pv(matrix.getAdjunctMatrix());pn;
    pv(matrix.getInverse());pn;
    pv((-2) * matrix + matrix);pn;
    pv(matrix * matrix.getInverse());pn;
    pv(matrix.getColumnCorrectMatrix(index, column));pn;
    pv(matrix.getInverseOfColumnCorrectedMatrix(index, column));pn;
    pv(matrix.getColumnCorrectMatrix(index, column) * matrix.getInverseOfColumnCorrectedMatrix(index, column));pn;
    pv(matrix | matrix);pn;
}

int main() {
    // test1();
    testSquare();

    return 0;
}