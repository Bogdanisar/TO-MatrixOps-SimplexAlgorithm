#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>


using namespace std;

// using NumberType = long double;
template<typename NumberType>
class Matrix {
public:
    static constexpr long double eps = 1e-5;

private:
    vector<vector<NumberType>> matrix;

    bool isMatrix(vector<vector<NumberType>> v) {
        if (v.size() == 0) {
            return false;
        }

        int N = v.size();
        int M = v[0].size();
        for (int i = 1; i < N; ++i) {
            if (v[i].size() != M) {
                return false;
            }
        }

        return true;
    }

public:
    Matrix() {
        this->matrix = {{1}};
    }

    Matrix(vector<vector<NumberType>> pmatrix) {
        assert(this->isMatrix(pmatrix));
        this->matrix = pmatrix;
    }

    Matrix(const Matrix& other) {
        this->matrix = other.matrix;
    }

    int get1Dim() const {
        return this->matrix.size();
    }

    int get2Dim() const {
        return this->matrix[0].size();
    }

    Matrix getTranspose() const {
        vector<vector<NumberType>> ans;

        int N = this->get1Dim();
        int M = this->get2Dim();
        for (int j = 0; j < M; ++j) {
            vector<NumberType> column;
            for (int i = 0; i < N; ++i) {
                column.push_back(this->matrix[i][j]);
            }

            ans.push_back(column);
        }

        return Matrix(ans);
    }

    static int getPermutationInversions(vector<int> permutation) {
        int inversions = 0;
        for (int i = 0; i < permutation.size(); ++i) {
            for (int j = i + 1; j < permutation.size(); ++j) {
                if (permutation[i] > permutation[j]) {
                    inversions += 1;
                }
            }
        }

        return inversions;
    }

    static int getPermutationSignature(vector<int> permutation) {
        int invs = Matrix::getPermutationInversions(permutation);
        return (invs % 2 == 0) ? 1 : -1;
    }

    bool isSquareMatrix() const {
        return this->get1Dim() == this->get2Dim();
    }

    NumberType getDeterminant() const {
        assert(this->isSquareMatrix());

        int N = this->get1Dim();

        vector<int> indexPermutation;
        indexPermutation.reserve(N);
        for (int i = 0; i < N; ++i) {
            indexPermutation.push_back(i);
        }

        NumberType result = 0;
        do {
            NumberType currentTerm = 1.0;
            for (int i = 0; i < N; ++i) {
                currentTerm *= this->matrix[i][indexPermutation[i]];
            }

            result += Matrix::getPermutationSignature(indexPermutation) * currentTerm;
        }
        while(next_permutation(indexPermutation.begin(), indexPermutation.end()));

        return result;
    }

    Matrix getSubmatrixWithoutRowAndColumn(int row, int column) const {
        vector<vector<NumberType>> ans;
        int N = this->get1Dim();
        int M = this->get2Dim();

        for (int i = 0; i < N; ++i) {
            if (i == row) {
                continue;
            }

            vector<NumberType> currentRow;
            for (int j = 0; j < M; ++j) {
                if (j == column) {
                    continue;
                }

                currentRow.push_back(this->matrix[i][j]);
            }

            ans.push_back(currentRow);
        }

        return Matrix(ans);
    }

    NumberType getCofactor(int i, int j) const {
        assert(this->isSquareMatrix());
        int N = this->get1Dim();
        assert(0 <= i && i < N && 0 <= j && j < N);

        Matrix subMatrix = this->getSubmatrixWithoutRowAndColumn(i, j);
        NumberType coef = ((i + 1 + j + 1) % 2 == 0) ? 1 : -1;
        return coef * subMatrix.getDeterminant();
    }

    Matrix getCofactorMatrix() const {
        assert(this->isSquareMatrix());

        int N = this->get1Dim();
        vector<vector<NumberType>> aux = vector<vector<NumberType>>(N, vector<NumberType>(N, 0));

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                aux[i][j] = this->getCofactor(i, j);
            }
        }

        return Matrix(aux);
    }

    Matrix getAdjunctMatrix() const {
        return this->getCofactorMatrix().getTranspose();
    }

    Matrix getInverse() const {
        NumberType determinant = this->getDeterminant();
        assert(("This is a singular matrix", abs(determinant) > 1e-10));

        Matrix inverse = this->getAdjunctMatrix();
        int N = inverse.get1Dim();
        int M = inverse.get2Dim();
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                inverse.matrix[i][j] /= determinant;
            }
        }

        return inverse;
    }

    Matrix getColumnCorrectMatrix(int columnIndex, vector<NumberType> column) const {
        int N = this->get1Dim();
        int M = this->get2Dim();
        assert(0 <= columnIndex && columnIndex < M);

        Matrix ret = *this;
        for (int i = 0; i < N; ++i) {
            ret.matrix[i][columnIndex] = column[i];
        }

        return ret;
    }

    // apply the substitution lemma
    Matrix getInverseOfColumnCorrectedMatrix(int columnIndex, vector<NumberType> column) const {
        assert(this->isSquareMatrix());

        int N = this->get1Dim();
        if (N != column.size()) {
            throw runtime_error("Column size is invalid in inverseOfColumnCorrectedMatrix");
        }
        if ( !(0 <= columnIndex && columnIndex < N) ) {
            throw runtime_error("Column index is not a valid value in inverseOfColumnCorrectedMatrix. Must be 0-indexed");
        }

        Matrix myInverse = this->getInverse();

        vector<NumberType> helperVector(N, 0);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                helperVector[i] += myInverse.matrix[i][j] * column[j];
            }
        }

        if (abs(helperVector[columnIndex]) < eps) {
            throw runtime_error("Can't apply lemma!. The index-th value in the Y vector is 0");
        }

        Matrix EMatrix = Matrix::getEye(N);
        for (int i = 0; i < N; ++i) {
            NumberType etaValue;

            if (i == columnIndex) {
                etaValue = 1 / helperVector[columnIndex];
            }
            else {
                etaValue = - helperVector[i] / helperVector[columnIndex];
            }

            EMatrix.matrix[i][columnIndex] = etaValue;
        }

        return EMatrix * myInverse;
    }

    Matrix getSubmatrixWithColumns(set<int> columns) const {
        int N = this->get1Dim();
        int M = this->get2Dim();

        for (int c : columns) {
            assert(0 <= c && c < M);
        }

        if (columns.size() == M) {
            return *this;
        }

        vector<vector<NumberType>> v(N, vector<NumberType>(columns.size(), 0));
        int currentColumn = 0;
        for (int j = 0; j < M; ++j) {
            if (columns.count(j) == 0) {
                continue;
            }

            for (int i = 0; i < N; ++i) {
                v[i][currentColumn] = this->matrix[i][j];
            }

            currentColumn += 1;
        }

        return Matrix(v);
    }

    Matrix getSubmatrixWithRows(set<int> rows) const {
        int N = this->get1Dim();
        int M = this->get2Dim();

        vector<vector<NumberType>> v;
        for (int i = 0; i < N; ++i) {
            if (rows.count(i) > 0) {
                v.push_back(this->matrix[i]);
            }
        }

        return Matrix(v);
    }








    static Matrix getEye(int N) {
        vector<vector<NumberType>> matrix(N, vector<NumberType>(N, 0));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    matrix[i][j] = 1;
                }
            }
        }

        return Matrix(matrix);
    }









    vector<NumberType>& operator [](int index) const {
        int N = this->get1Dim();
        assert(0 <= index && index < N);

        return this->matrix[index];
    }

    friend Matrix operator +(const Matrix& A, const Matrix& B) {
        int N = A.get1Dim();
        int M = A.get2Dim();
        assert(( "Matrixes must have the same dimensions", N == B.get1Dim() && M == B.get2Dim() ));

        Matrix ret = A;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                ret.matrix[i][j] += B.matrix[i][j];
            }
        }

        return ret;
    }

    friend Matrix operator *(NumberType scalar, Matrix sm) {
        int N = sm.get1Dim();
        int M = sm.get2Dim();

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                sm.matrix[i][j] *= scalar;
            }
        }

        return sm;
    }

    friend Matrix operator -(const Matrix& A, const Matrix& B) {
        int N = A.get1Dim();
        int M = A.get2Dim();
        assert(( "Matrixes must have the same dimensions", N == B.get1Dim() && M == B.get2Dim() ));

        return A + (NumberType)(-1) * B;
    }

    friend Matrix operator *(const Matrix& A, const Matrix& B) {
        int N = A.get1Dim();
        int M = A.get2Dim();
        int O = B.get2Dim();
        assert(( "Matrixes must have a common dimension in product", M == B.get1Dim() ));

        Matrix ret = A;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < O; ++j) {
                ret.matrix[i][j] = 0;

                for (int k = 0; k < M; ++k) {
                    ret.matrix[i][j] += A.matrix[i][k] * B.matrix[k][j];
                }
            }
        }

        return ret;
    }

    friend Matrix operator |(const Matrix& A, const Matrix& B) {
        int N = A.get1Dim();
        int M = A.get2Dim();
        int O = B.get2Dim();
        assert(( 
            "Can't append matrixes if they don't have the same number of rows",
            N == B.get1Dim()
        ));

        vector<vector<NumberType>> v(N, vector<NumberType>(M + O));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                v[i][j] = A.matrix[i][j];
            }

            for (int j = 0; j < O; ++j) {
                v[i][M + j] = B.matrix[i][j];
            }
        }

        return Matrix(v);
    }

    friend ostream& operator <<(ostream& out, const Matrix& matrix) {
        int N = matrix.get1Dim();
        int M = matrix.get2Dim();

        out << N << ' ' << M << '\n';
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                out << matrix.matrix[i][j] << ' ';
            }
            out << '\n';
        }

        return out;
    }

    
    static Matrix readMatrixFromStream(istream& in) {
        int N, M;
        in >> N >> M;
        vector<vector<NumberType>> v(N, vector<NumberType>(M, 0));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                in >> v[i][j];
            }
        }

        return Matrix(v);
    }
};

    