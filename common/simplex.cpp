#include "simplex.h"
#include <limits>
#include <set>

using ld = long double;
using namespace std;

const long double EPS = 1e-10;
const long double LARGE_VALUE = 1e4;

enum class PivotResult {
    NO_SOLUTION, OPTIMAL_SOLUTION, INFINITE_SOLUTION, FOUND_PIVOT
};





// metode comune si pentru simplex primal si pentru simplex dual

void printSimplexState(SimplexState state) {
    string sep = "==================================";
    cout << "\n\n" << sep << '\n';

    int N = state.z_c.get2Dim();

    cout << "VVB | Y:\n";
    printMatrix(cout, state.BbA);
    cout << '\n';

    cout << "z: " << state.z << ' ';
    for (int j = 0; j < N; ++j) {
        cout << state.z_c[0][j] << ' ';
    }
    cout << "\n";

    cout << "basis: ";
    for (int v : state.basis) {
        cout << v << ' ';
    }
    cout << '\n';

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



// seteaza state.z (valoarea functiei obiectiv) si state.z_c (zj - cj)
void setBottomRowInState(const Matrix<ld>& A, const Matrix<ld>& c, SimplexState& state) {
    int M = A.get1Dim();

    // seteaza state.z
    state.z = 0;
    for (int i = 0; i < M; ++i) {
        state.z += c[0][state.basis[i]] * state.BbA[i][0];
    }

    // seteaza state.z_c
    vector<vector<ld>> v_cB(M);
    for (int i = 0; i < M; ++i) {
        v_cB[i].push_back(c[0][state.basis[i]]);
    }
    Matrix<ld> cB(v_cB);
    state.z_c = (state.BbA.getTranspose() * cB).getSubmatrixWithoutRows({0}).getTranspose();
    state.z_c = state.z_c - c;
}

// pivotul ar trebui dat ca o pozitie in matricea BbA (dar nu pe prima coloana)
SimplexState changeStateWithPivot(const Matrix<ld>& A, const Matrix<ld>& c, SimplexState state, pair<int,int> pivot) {
    assertState(A, c, state);

    int M = A.get1Dim();
    int N = A.get2Dim();
    assert(0 <= pivot.first && pivot.first < M);
    assert(0 < pivot.second && pivot.second < N + 1);

    SimplexState newState;

    // seteaza newState.BbA
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

    // seteaza newState.basis
    newState.basis = state.basis;
    newState.basis[pivot.first] = pivot.second - 1;

    // seteaza newState.z si newState.z_c;
    setBottomRowInState(A, c, newState);

    assertState(A, c, newState);
    return newState;
}

SimplexState buildInitialState(const Matrix<ld>& A, const Matrix<ld>& b, const Matrix<ld>& c, const vector<int>& basis) {
    assert(set<int>(basis.begin(), basis.end()).size() == basis.size());
    assert(basis.size() == A.get1Dim());

    Matrix<ld> B = A.getSubmatrixWithColumns(basis);
    assert(( "The given variables don't form a basis", abs(B.getDeterminant()) > EPS ));

    SimplexState state;
    state.BbA = B.getInverse() * (b | A);
    state.basis = basis;
    setBottomRowInState(A, c, state);

    assertState(A, c, state);

    return state;
}

vector<ld> getVVB(const Matrix<ld>& BbA) {
    int M = BbA.get1Dim();

    vector<ld> vvb;
    for (int i = 0; i < M; ++i) {
        vvb.push_back(BbA[i][0]);
    }

    return vvb;
}

























// metode pentru algoritmul simplex primal

pair<PivotResult, pair<int,int>> getPrimalPivot(const Matrix<ld>& A, const Matrix<ld>& c, const SimplexState& state) {
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
            PivotResult::OPTIMAL_SOLUTION,
            {}
        };
    }

    for (int k : R_plus) { // k e indexul unei variabile si e indexat de la 0
        bool allNonPositive = true;
        for (int i = 0; i < M; ++i) {
            if (state.BbA[i][k + 1] > 0) {
                allNonPositive = false;
                break;
            }
        }

        if (allNonPositive) {
            return {
                PivotResult::INFINITE_SOLUTION,
                {}
            };
        }
    }


    // aplica regula lui Bland
    int k = R_plus[0]; // k e indexul unei variabile si e indexat de la 0
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
        PivotResult::FOUND_PIVOT,
        {lineIndex, k + 1} // intoarce pivotul, dar ca o pozitie in matricea BbA;
    };
}



// executa algoritmul simplex primal pe o baza admisibila specificata si cu o functie obiectiv de minim
SimplexReturnType runSimplexWithMinObjective(
    Matrix<ld> A,
    Matrix<ld> b,
    Matrix<ld> c,
    vector<int> basis
) {
    int M = A.get1Dim();
    int N = A.get2Dim();

    // verifica ca baza e corecta
    assert(basis.size() == M);
    assert(set<int>(basis.begin(), basis.end()).size() == basis.size());

    Matrix<ld> B = A.getSubmatrixWithColumns(basis);
    assert(( "The given variables don't form a basis", abs(B.getDeterminant()) > EPS ));

    // seteaza starea simplex initiala (tabelul)
    SimplexState state = buildInitialState(A, b, c, basis);

    cout << "Starting min-objective primal simplex with the following state: \n";

    const int iterationLimit = 1e5;
    int currentIteration = 0;
    while (currentIteration++ < iterationLimit) {

        printSimplexState(state);

        auto temp = getPrimalPivot(A, c, state);
        PivotResult result = temp.first;
        assert(result != PivotResult::NO_SOLUTION);
        pair<int,int> pivot = temp.second;

        cout << "Pivot: " << pivot.first + 1 << ' ' << pivot.second << "\n\n";

        if (result == PivotResult::INFINITE_SOLUTION) {
            return SimplexReturnType{
                SimplexResult::INFINITE_SOLUTION,
                getVVB(state.BbA),
                state
            };
        }
        if (result == PivotResult::OPTIMAL_SOLUTION) {
            return SimplexReturnType{
                SimplexResult::OPTIMAL_SOLUTION,
                getVVB(state.BbA),
                state
            };
        }

        state = changeStateWithPivot(A, c, state, pivot);
    }

    assert(( "Stopping the primal simplex algorithm after reaching the iteration limit. Did it cycle?", false ));
}



// executa algoritmul simplex primal pe o baza admisibila specificata (si o functie obiectiv de minim sau maxim)
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



// executa algoritmul simplex (pe o functie obiectiv de minim sau maxim) fara a specifica o baza admisibila initiala
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

        // avem o variabila artificiala in baza;
        // se incearca eliminarea ei printr-o pivotare;
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




    // construieste remainingRows;
    set<int> remainingRows;
    for (int i = 0; i < M; ++i) {
        remainingRows.insert(i);
    }
    for (int x : artificialThatCouldntBeRemoved) {
        cout << "Removing row " << x - N + 1 << "...\n";
        remainingRows.erase(x - N);
    }

    // construieste newA (eliminand din A unele randuri);
    Matrix<ld> newA = A.getSubmatrixWithRows(remainingRows);

    // construieste newB (eliminand din b aceleasi randuri);
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

















// metode pentru algoritmul simplex dual

pair<PivotResult, pair<int,int>> getDualPivot(const Matrix<ld>& A, const Matrix<ld>& c, const SimplexState& state) {
    assertState(A, c, state);

    int M = A.get1Dim();
    int N = A.get2Dim();

    int L = -1;
    long double minVal;
    for (int i = 0; i < M; ++i) {
        if (state.BbA[i][0] < 0) {
            if (L == -1 || minVal > state.BbA[i][0]) {
                L = i;
                minVal = state.BbA[i][0];
            }
        }
    }

    if (L == -1) { // toate valorile au fost >= 0
        return {
            PivotResult::OPTIMAL_SOLUTION,
            {}
        };
    }

    int k = -1; // k e un index de variabila si e indexat de la 0
    long double minRatio;
    for (int j = 1; j < N + 1; ++j) {
        if (state.BbA[L][j] < 0) {
            long double ratio = abs(state.z_c[0][j - 1] / state.BbA[L][j]);
            if (k == -1 || minRatio > ratio) {
                k = j - 1;
                minRatio = ratio;
            }
        }
    }

    if (k == -1) {
        return {
            PivotResult::NO_SOLUTION,
            {}
        };
    }

    return {
        PivotResult::FOUND_PIVOT,
        {L, k + 1} // intoarce pivotul, dar ca o pozitie in matricea BbA;
    };
}


// executa algoritmul simplex dual cu o baza dual admisibila specificata si cu o functie obiectiv de minim
SimplexReturnType runDualSimplexWithDualBasisAndMinObjective(
    Matrix<ld> A,
    Matrix<ld> b,
    Matrix<ld> c,
    vector<int> basis
) {
    int M = A.get1Dim();
    int N = A.get2Dim();

    // verifica ca baza e corecta
    assert(basis.size() == M);
    assert(set<int>(basis.begin(), basis.end()).size() == basis.size());

    Matrix<ld> B = A.getSubmatrixWithColumns(basis);
    assert(( "The given variables don't form a basis", abs(B.getDeterminant()) > EPS ));

    // seteaza starea simplex initiala (tabelul)
    SimplexState state = buildInitialState(A, b, c, basis);

    cout << "Starting min-objective dual simplex with the following state: \n";

    const int iterationLimit = 1e5;
    int currentIteration = 0;
    while (currentIteration++ < iterationLimit) {

        printSimplexState(state);

        auto temp = getDualPivot(A, c, state);
        PivotResult result = temp.first;
        assert(result != PivotResult::INFINITE_SOLUTION);
        pair<int,int> pivot = temp.second;

        cout << "Pivot: " << pivot.first + 1 << ' ' << pivot.second << "\n\n";

        if (result == PivotResult::OPTIMAL_SOLUTION) {
            return SimplexReturnType{
                SimplexResult::OPTIMAL_SOLUTION,
                getVVB(state.BbA),
                state
            };
        }

        if (result == PivotResult::NO_SOLUTION) {
            return SimplexReturnType{
                SimplexResult::NO_SOLUTION,
                getVVB(state.BbA),
                state
            };
        }

        state = changeStateWithPivot(A, c, state, pivot);
    }

    assert(( "Stopping the dual simplex algorithm after reaching the iteration limit. Did it cycle?", false ));
}



// executa algoritmul simplex dual cu o baza dual admisibila specificata si cu o functie obiectiv de minim / maxim
SimplexReturnType runDualSimplexWithDualBasis(
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
        return runDualSimplexWithDualBasisAndMinObjective(A, b, c, basis);
    }
    else {
        SimplexReturnType ret = runDualSimplexWithDualBasisAndMinObjective(A, b, (-1) * c, basis);
        ret.state.z *= (-1);
        return ret;
    }
}




// executa algoritmul simplex dual cu o baza oarecare specificata si cu o functie obiectiv de minim
SimplexReturnType runDualSimplexWithAnyBasisAndMinObjective(
    Matrix<ld> A,
    Matrix<ld> b,
    Matrix<ld> c,
    vector<int> basis
) {
    int M = A.get1Dim();
    int N = A.get2Dim();

    assert(1 == c.get1Dim());
    assert(N == c.get2Dim());
    assert(M == b.get1Dim());
    assert(1 == b.get2Dim());

    // verifica ca baza e corecta
    assert(basis.size() == M);
    assert(set<int>(basis.begin(), basis.end()).size() == basis.size());
    for (int idx : basis) {
        assert(0 <= idx && idx < N);
    }

    Matrix<ld> B = A.getSubmatrixWithColumns(basis);
    assert(( "The given variables don't form a basis", abs(B.getDeterminant()) > EPS ));
    set<int> basisSet(basis.begin(), basis.end()); // baza ca un set; va contine indici indexati de la 0


    // verifica daca baza este chiar dual admisibila;
    // daca nu, obtine coloana k din algoritmul de determinare a unei baze dual admisibile.
    SimplexState aux_state = buildInitialState(A, b, c, basis);

    int k = -1; // k e un index de variabila si este indexat de la 0
    for (int j = 0; j < N; ++j) {
        if (basisSet.count(j) == 0) { // variabila j nu este in baza
            if (k == -1 || aux_state.z_c[0][k] < aux_state.z_c[0][j]) {
                k = j;
            }
        }
    }

    if (aux_state.z_c[0][k] <= 0) { // baza data este dual admisibila
        cout << "The given basis is actually dually feasible. Applying dual simplex...\n";
        return runDualSimplexWithDualBasisAndMinObjective(A, b, c, basis);
    }

    cout << "The initial k (1-indexed) needed to find the dual basis is " << k + 1 << '\n';

    // construieste A1 pentru problema liniara asociata
    Matrix<ld> A1(M + 1, N + 1);
    A1[0][0] = 1;
    for (int i = 1; i < M + 1; ++i) {
        A1[i][0] = 0;
    }

    for (int j = 0; j < N; ++j) {
        if (basisSet.count(j) == 0) { // variabila nu este in baza, deci pe primul rand introdus punem valoarea 1
            A1[0][j + 1] = 1;
        }
    }

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            A1[i + 1][j + 1] = A[i][j];
        }
    }

    assert(A1.getSubmatrixWithoutRowAndColumn(0, 0) == A);


    // construieste b1 pentru problema liniara asociata
    Matrix<ld> b1(M + 1, 1);
    b1[0][0] = LARGE_VALUE; // o valoare suficient de mare (M)
    for (int i = 0; i < M; ++i) {
        b1[i + 1][0] = b[i][0];
    }


    // construieste c1 pentru problema liniara asociata
    Matrix<ld> c1(1, N + 1);
    c1[0][0] = 0; // costul lui x0 (variabila artificiala introdusa)
    for (int j = 0; j < N; ++j) {
        c1[0][j + 1] = c[0][j];
    }



    // construieste o noua baza adaugand variabila k la inceputul bazei date
    vector<int> newBasis = basis;
    newBasis.push_back(-100);
    for (int i = newBasis.size() - 1; i > 0; --i) {
        newBasis[i] = newBasis[i - 1];
    }
    newBasis[0] = k;

    // toate variabilele din newBasis sunt neartificiale, dar incep de la 0,
    // asa ca le crestem valoarea cu 1
    for (int& val : newBasis) {
        val += 1;
    }
    

    cout << "newBasis: " << '\n';
    for (int val : newBasis) {
        cout << val << ' ';
    }
    cout << "\n\n";


    SimplexReturnType ret = runDualSimplexWithDualBasisAndMinObjective(A1, b1, c1, newBasis);
    assert(ret.result != SimplexResult::INFINITE_SOLUTION);

    cout << "The result of the dual simplex on the extended problem is:\n";
    printSimplexState(ret.state);

    if (ret.result == SimplexResult::NO_SOLUTION) {
        cout << "The extended problem did not have a solution so the initial problem does not have a solution\n";
        return {
            SimplexResult::NO_SOLUTION,
            {},
            {}
        };
    }

    // am gasit o solutie optima pentru problema asociata


    set<int> retBasisSet(ret.state.basis.begin(), ret.state.basis.end());
    set<int> retBasisSetWithout0 = retBasisSet;
    retBasisSetWithout0.erase(0);


    if (retBasisSet.count(0) > 0) { // variabila artificiala x0 este in baza gasita

        // atunci restul bazei este o baza optima pentru problema initiala
        vector<int> basisForInitialProblem;
        for (int idx : retBasisSetWithout0) {
            basisForInitialProblem.push_back(idx - 1);
        }

        SimplexState resultingState = buildInitialState(A, b, c, basisForInitialProblem);

        cout << "The artificial variable was in the returned basis.\n";
        cout << "so the rest of the basis is an optimal basis for the initial problem.\n";
        return {
            SimplexResult::OPTIMAL_SOLUTION,
            getVVB(resultingState.BbA),
            resultingState
        };
    }



    // gaseste coeficientul lui LARGE_VALUE (M) din expresia functiei obiectiv
    long double coefficient;
    {
        // obtine coeficientii din c corespunzatori variabilelor din baza
        Matrix<ld> cB1W(M + 1, 1);
        for (int i = 0; i < M + 1; ++i) {
            cB1W[i][0] = c1[0][ret.state.basis[i]];
        }

        Matrix<ld> B1W = A1.getSubmatrixWithColumns( vector<int>(retBasisSet.begin(), retBasisSet.end()) );

        coefficient = (cB1W.getTranspose() * B1W)[0][0];
    }

    cout << "M's coefficient is " << coefficient << '\n';
    
    if (abs(coefficient) < EPS) { // coeficientul lui M este 0

        vector<int> basisForInitialProblem;
        for (int idx : retBasisSet) {
            basisForInitialProblem.push_back(idx - 1);
        }

        cout << "The coefficient was 0,\n";
        cout << "so the found optimal value of the extended problem\n";
        cout << "is also an optimal value for the initial problem\n";

        return {
            SimplexResult::OPTIMAL_SOLUTION,
            getVVB(ret.state.BbA),
            ret.state
        };
    }
    else {
        return {
            SimplexResult::INFINITE_SOLUTION,
            {},
            {}
        };
    }
}



// executa algoritmul simplex dual pe o baza oarecare si cu o functie obiectiv de minim sau de maxim
SimplexReturnType runDualSimplexWithAnyBasis(
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
        return runDualSimplexWithAnyBasisAndMinObjective(A, b, c, basis);
    }
    else {
        SimplexReturnType ret = runDualSimplexWithAnyBasisAndMinObjective(A, b, (-1) * c, basis);
        ret.state.z *= (-1);
        return ret;
    }
}
