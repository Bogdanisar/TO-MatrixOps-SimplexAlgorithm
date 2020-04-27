Pentru tema 3, am implementat:
a) algoritmul simplex primal pentru situatia cand se cunoaste o baza admisibila;
b) algoritmul simplex primal pentru situatia cand nu se cunoaste o astfel de baza; (asa cum este descris in cursul 5)

Implementarile propriu-zise se afla in fisierul ../common/simplex.cpp
Fisierul simplexWithBasis.cpp rezolva situatia a) si citeste din simplexWithBasis.in
Fisierul simplexWithoutBasis.cpp rezolva situatia b) si citeste din simplexWithoutBasis.in
Toate problemele liniare se citesc in forma standard


Pentru compilare, se executa din directorul tema3:
g++ simplexWithBasis.cpp -o simplexWithBasis.exe
g++ simplexWithoutBasis.cpp -o simplexWithoutBasis.exe

Apoi, pentru executare:
./simplexWithBasis.exe
./simplexWithoutBasis.exe




========================== Formatul fisierului simplexWithBasis.in =========================
M (numarul de ecuatii) N (numarul de variabile)
* matricea A a formei standard * (de dimensiune MxN)
* b * (ca vector coloana cu M elemente)
"min" sau "max"
* c * (ca vector linie cu N elemente)
* baza initiala * (ca vector linie cu M elemente)

Baza trebuie sa contina M variabile distincte (indexate de la 1) si sa fie admisibila


========================== Formatul fisierului simplexWithoutBasis.in =========================
M (numarul de ecuatii) N (numarul de variabile)
* matricea A a formei standard * (de dimensiune MxN)
* b * (ca vector coloana cu M elemente)
"min" sau "max"
* c * (ca vector linie cu N elemente)
