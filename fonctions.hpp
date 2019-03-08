#ifndef LDPC_FONCTIONS_HPP
#define LDPC_FONCTIONS_HPP

#include <vector>
#include <iostream>

using Row = std::vector< int >;
using Matrice = std::vector< Row >;
using RowIt = Row::iterator;

/// @struct Params
/// @var n nombre de colonnes (largeur)
/// @var i nombre de 1 par colonne
/// @var j nombre de 1 par ligne
/// @var k nombre de lignes (hauteur)
struct Params {
    Params(int n, int j, int k) : n(n), i(k*j/n), j(j), k(k) {}
    Params(int n, int k, float density) : n(n), i(density*k), j(density*n), k(k) {}
    int n,i,j,k;
};

std::ostream& operator<< (std::ostream& os, const Matrice& mat);
Matrice operator* (const Matrice& lhs, const Matrice& rhs);
Matrice toFull(const Matrice& mat, int rowSize);
int randomInt(int min, int max);
void shuffle_n(Row& row, int n);
bool estValide(const Matrice& mat, Params p);

void swapColumns(Matrice& mat, int r, int x);
void transvection(Row& row_r, Row& row_x);
void annulation(Matrice& mat, int r);

Matrice gauss_jordan(Matrice mat);
bool isFullRank(const Matrice& mat);

Matrice gallager(Params p);
Matrice macKayNeal(Params p);

Matrice transpose(const Matrice& mat);
void comp2(Matrice& mat);

Matrice getP(Matrice H, Params p);

#endif //LDPC_FONCTIONS_HPP
