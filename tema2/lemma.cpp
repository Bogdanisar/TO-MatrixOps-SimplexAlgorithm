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


// name=lemma; g++ $name.cpp -o $name.exe && ./$name.exe


using Numeric = long double;


void showLemma() {
    ifstream fin("lemma.in");
    Matrix<Numeric> matrix = Matrix<Numeric>::readMatrixFromStream(fin);
    int index;
    fin >> index;

    int N = matrix.get1Dim();
    vector<long double> column(N);
    for (int i = 0; i < N; ++i) {
        fin >> column[i];
    }

    cout << "Matricea:\n";
    cout << matrix << "\n\n";
    
    cout << "Determinantul = " << matrix.getDeterminant() << "\n";

    cout << "Inversa:\n";
    cout << matrix.getInverse() << '\n';

    cout << "Matrice * Inversa:\n";
    cout << matrix * matrix.getInverse() << '\n';

    cout << "Matricea alterata:\n";
    cout << matrix.getColumnCorrectMatrix(index, column) << '\n';

    cout << "Inversa matricei alterate prin lema substitutiei:\n";
    cout << matrix.getInverseOfColumnCorrectedMatrix(index, column) << '\n';

    cout << "Matricea alterata * Inversa ei:\n";
    cout << matrix.getColumnCorrectMatrix(index, column) * matrix.getInverseOfColumnCorrectedMatrix(index, column) << '\n';
}

int main() {
    showLemma();

    return 0;
}