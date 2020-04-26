#include "simplex.h"
#include <limits>
#include <set>

using ld = long double;
using namespace std;

const long double EPS = 1e-10;

void printSimplexState(SimplexState state) {
    string sep = "==================================";
    cout << "\n\n" << sep << '\n';

    int N = state.z_c.get2Dim();

    cout << state.BbA;
    for (int v : state.basis) {
        cout << v << ' ';
    }
    cout << '\n' << state.z << ' ';
    for (int j = 0; j < N; ++j) {
        cout << state.z_c[0][j] << ' ';
    }
    cout << "\n";

    cout << sep << endl;
}

void assertState(const Matrix<ld>& A, const Matrix<ld>& c, const SimplexState& state) {
    int M = A.get1Dim();
    int N = A.get2Dim();

    assert(1 == c.get1Dim());
    assert(N == c.get2Dim());

    assert(M == state.BbA.get1Dim());
    assert(N + 1 == state.BbA.get2Dim());

    assert(1 == state.z_c.get1Dim());
    assert(N == state.z_c.get2Dim());

    assert(M == state.basis.size());
    for (int idx : state.basis) {
        assert(0 <= idx && idx < N);
    }
}


pair<SimplexResult, pair<int,int>> getPivot(const Matrix<ld>& A, const Matrix<ld>& c, const SimplexState& state) {
    assertState(A, c, state);

    int M = A.get1Dim();
    int N = A.get2Dim();

    set<int> basisSet(state.basis.begin(), state.basis.end());

    vector<int> R_plus;
    for (int j = 0; j < N; ++j) {
        if (basisSet.count(j) == 0 && state.z_c[0][j] > 0) {
            R_plus.push_back(j);
        }
    }

    if (R_plus.size() == 0) {
        return {
            SimplexResult::OPTIMAL_SOLUTION,
            {}
        };
    }

    for (int k : R_plus) {
        bool allNonPositive = true;
        for (int i = 0; i < M; ++i) {
            if (state.BbA[i][k + 1] > 0) {
                allNonPositive = false;
                break;
            }
        }

        if (allNonPositive) {
            return {
                SimplexResult::INFINITE_SOLUTION,
                {}
            };
        }
    }


    // apply Bland's rule:
    int k = R_plus[0];
    ld minumum_ratio = std::numeric_limits<ld>::infinity();
    int lineIndex = -1, varIndex = -1;
    for (int i = 0; i < M; ++i) {
        if (state.BbA[i][k + 1] <= 0) {
            continue;
        }

        ld ratio = state.BbA[i][0] / state.BbA[i][k + 1];
        if (minumum_ratio > ratio) {
            minumum_ratio = ratio;
            varIndex = state.basis[i];
            lineIndex = i;
        }
        else if (minumum_ratio == ratio && varIndex > state.basis[i]) {
            varIndex = state.basis[i];
            lineIndex = i;
        }
    }

    return {
        SimplexResult::NO_SOLUTION,
        {lineIndex, k + 1}
    };
}

// set state.z and state.z_c
void setBottomRowInState(const Matrix<ld>& A, const Matrix<ld>& c, SimplexState& state) {
    int M = A.get1Dim();

    // set state.z
    state.z = 0;
    for (int i = 0; i < M; ++i) {
        state.z += c[0][state.basis[i]] * state.BbA[i][0];
    }

    // set state.z_c
    vector<vector<ld>> v_cB(M);
    for (int i = 0; i < M; ++i) {
        v_cB[i].push_back(c[0][state.basis[i]]);
    }
    Matrix<ld> cB(v_cB);
    state.z_c = (state.BbA.getTranspose() * cB).getSubmatrixWithoutRows({0}).getTranspose();
    state.z_c = state.z_c - c;
}

SimplexState changeStateWithPivot(const Matrix<ld>& A, const Matrix<ld>& c, SimplexState state, pair<int,int> pivot) {
    assertState(A, c, state);

    int M = A.get1Dim();
    int N = A.get2Dim();
    assert(0 <= pivot.first && pivot.first < M);
    assert(0 <= pivot.second && pivot.second < N + 1);

    SimplexState newState;

    // set newState.BbA
    newState.BbA = state.BbA;
    for (int i = 0; i < M; ++i) {
        if (i == pivot.first) {
            for (int j = 0; j < N + 1; ++j) {
                newState.BbA[i][j] = state.BbA[i][j] / state.BbA[pivot.first][pivot.second];
            }
        }
        else {
            for (int j = 0; j < N + 1; ++j) {
                newState.BbA[i][j] = state.BbA[i][j] * state.BbA[pivot.first][pivot.second];
                newState.BbA[i][j] -= state.BbA[pivot.first][j] * state.BbA[i][pivot.second];
                newState.BbA[i][j] /= state.BbA[pivot.first][pivot.second];
            }
        }
    }

    // set newState.basis
    newState.basis = state.basis;
    newState.basis[pivot.first] = pivot.second - 1;

    // set newState.z and newState.z_c;
    setBottomRowInState(A, c, newState);

    assertState(A, c, newState);
    return newState;
}

vector<ld> getVVB(const Matrix<ld>& BbA) {
    int M = BbA.get1Dim();

    vector<ld> vvb;
    for (int i = 0; i < M; ++i) {
        vvb.push_back(BbA[i][0]);
    }

    return vvb;
}



// run the simplex algorithm with a specified starting basis and a min objective function
SimplexReturnType runSimplexWithMinObjective(
    Matrix<ld> A,
    Matrix<ld> b,
    Matrix<ld> c,
    vector<int> basis
) {
    int M = A.get1Dim();
    int N = A.get2Dim();

    // check that the basis is correct
    assert(M == basis.size());
    set<int> columns;
    for (int idx : basis) {
        columns.insert(idx);
    }
    Matrix<ld> B = A.getSubmatrixWithColumns(columns);
    assert(( "The given variables don't form a basis", abs(B.getDeterminant()) > EPS ));

    // set initial SimplexState
    SimplexState state;
    state.BbA = B.getInverse() * (b | A);
    state.basis = basis;
    setBottomRowInState(A, c, state);

    assertState(A, c, state);

    cout << "Starting simplex with the following state: \n";

    const int iterationLimit = 1e5;
    int currentIteration = 0;
    while (currentIteration++ < iterationLimit) {

        printSimplexState(state);

        auto temp = getPivot(A, c, state);
        SimplexResult result = temp.first;
        pair<int,int> pivot = temp.second;

        cout << "Pivot: " << pivot.first + 1 << ' ' << pivot.second << "\n\n";

        if (result == SimplexResult::INFINITE_SOLUTION || result == SimplexResult::OPTIMAL_SOLUTION) {
            return SimplexReturnType{
                result,
                getVVB(state.BbA),
                state
            };
        }

        state = changeStateWithPivot(A, c, state, pivot);
    }

    assert(( "Stopping simplex algorithm after reaching the iteration limit. Did it cycle?", false ));
}



// run the simplex algorithm with a specified starting basis
SimplexReturnType runSimplex(
    Matrix<ld> A,
    Matrix<ld> b,
    Matrix<ld> c,
    TypeOfObjective obj,
    vector<int> basis
) {
    int M = A.get1Dim();
    int N = A.get2Dim();
    assert(M <= N);

    if (obj == TypeOfObjective::MIN) {
        return runSimplexWithMinObjective(A, b, c, basis);
    }
    else {
        SimplexReturnType ret = runSimplexWithMinObjective(A, b, (-1) * c, basis);
        ret.state.z *= (-1);
        return ret;
    }

}



// run the simplex algorithm without specifying a starting basis
SimplexReturnType runSimplex(
    Matrix<ld> A,
    Matrix<ld> b,
    Matrix<ld> c,
    TypeOfObjective obj
) {
    int M = A.get1Dim();
    int N = A.get2Dim();

    assert(1 == c.get1Dim());
    assert(N == c.get2Dim());
    assert(M == b.get1Dim());
    assert(1 == b.get2Dim());

    for (int i = 0; i < M; ++i) {
        if (b[i][0] < 0) {
            b[i][0] *= - 1;
            for (int j = 0; j < N; ++j) {
                A[i][j] *= (-1);
            }
        }
    }

    Matrix<ld> auxA = (A | Matrix<ld>::getEye(M));

    vector<vector<ld>> v_auxC(1, vector<ld>(N + M, 0));
    for (int j = N; j < N + M; ++j) {
        v_auxC[0][j] = 1;
    }
    Matrix<ld> auxC = Matrix<ld>(v_auxC);

    vector<int> auxBasis;
    for (int j = N; j < N + M; ++j) {
        auxBasis.push_back(j);
    }

    SimplexReturnType ret = runSimplexWithMinObjective(auxA, b, auxC, auxBasis);
    assert(ret.result != SimplexResult::NO_SOLUTION && ret.result != SimplexResult::INFINITE_SOLUTION);

    cout << "The first simplex of the two phase method has this result:\n";
    printSimplexState(ret.state);
    
    if (abs(ret.state.z) > EPS) {
        return SimplexReturnType{
            SimplexResult::NO_SOLUTION,
            {},
            {}
        };
    }

    bool allArtificialNotInBasis = true;
    for (int idx : ret.state.basis) {
        if (idx >= N) {
            allArtificialNotInBasis = false;
            break;
        }
    }

    if (allArtificialNotInBasis) {
        cout << "No artificial variables found in basis!\n";
        cout << "Rerunning simplex with the found basis...\n";
        return runSimplex(A, b, c, obj, ret.state.basis);
    }

    set<int> artificialThatCouldntBeRemoved;
    for (int i = 0; i < M; ++i) {
        if (ret.state.basis[i] < N) {
            continue;
        }

        // we have an artificial variable in the basis
        // try to remove it by doing a pivot
        int k = -1;
        for (int j = 1; j < N + 1; ++j) {
            if (abs(ret.state.BbA[i][j]) > EPS) {
                k = j - 1;
                break;
            }
        }

        if (k == -1) {
            artificialThatCouldntBeRemoved.insert(ret.state.basis[i]);
        }
        else {
            ret.state = changeStateWithPivot(auxA, auxC, ret.state, {i, k + 1});
        }
    }

    if (artificialThatCouldntBeRemoved.size() == 0) {
        cout << "All artificial variables where removed from the basis!\n";
        cout << "Rerunning simplex with the found basis...\n";
        return runSimplex(A, b, c, obj, ret.state.basis);
    }

    cout << "Some artificial variables couldn't be removed\n";

    vector<int> newBasis;
    for (int idx : ret.state.basis) {
        if (idx < N) {
            newBasis.push_back(idx);
        }
    }




    // build remainingRows;
    set<int> remainingRows;
    for (int i = 0; i < M; ++i) {
        remainingRows.insert(i);
    }
    for (int x : artificialThatCouldntBeRemoved) {
        cout << "Removing row " << x - N + 1 << "...\n";
        remainingRows.erase(x - N);
    }

    // build newA (remove the rows);
    Matrix<ld> newA = A.getSubmatrixWithRows(remainingRows);

    // build newB (remove the rows);
    vector<vector<ld>> v_newB;
    for (int i = 0; i < M; ++i) {
        if (remainingRows.count(i) > 0) {
            vector<ld> aux;
            aux.push_back(b[i][0]);
            v_newB.push_back(aux);
        }
    }
    Matrix<ld> newB = Matrix<ld>(v_newB);


    return runSimplex(newA, newB, c, obj, newBasis);
}