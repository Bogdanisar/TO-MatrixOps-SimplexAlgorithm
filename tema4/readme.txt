Pentru tema 4, am implementat:
a) algoritmul simplex dual pentru situatia cand se cunoaste o baza dual admisibila;
b) algoritmul simplex dual pentru situatia cand se cunoaste o baza oarecare; (asa cum este descris in "Baza dual admisibila.pdf")

Implementarile propriu-zise se afla in fisierul ../common/simplex.cpp
Fisierul dualWithDualBasis.cpp rezolva situatia a) si citeste din dualWithDualBasis.in
Fisierul dualWithAnyBasis.cpp rezolva situatia b) si citeste din dualWithAnyBasis.in
Toate problemele liniare se citesc in forma standard


Pentru compilare, se executa din directorul tema4:
g++ dualWithDualBasis.cpp -o dualWithDualBasis.exe
g++ dualWithAnyBasis.cpp -o dualWithAnyBasis.exe

Apoi, pentru executare:
./dualWithDualBasis.exe
./dualWithAnyBasis.exe




========================== Formatul fisierului dualWithDualBasis.in =========================
M (numarul de ecuatii) N (numarul de variabile)
* matricea A a formei standard * (de dimensiune MxN)
* b * (ca vector coloana cu M elemente)
"min" sau "max"
* c * (ca vector linie cu N elemente)
* baza initiala * (ca vector linie cu M elemente)

Baza trebuie sa contina M variabile distincte (indexate de la 1) si sa fie dual admisibila


========================== Formatul fisierului dualWithAnyBasis.in =========================
M (numarul de ecuatii) N (numarul de variabile)
* matricea A a formei standard * (de dimensiune MxN)
* b * (ca vector coloana cu M elemente)
"min" sau "max"
* c * (ca vector linie cu N elemente)
* baza initiala * (ca vector linie cu M elemente)

Baza trebuie sa contina M elemente distincte (indexate de la 1)
