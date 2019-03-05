#ifndef LDPC_FONCTIONS_HPP
#define LDPC_FONCTIONS_HPP

#include <vector>
#include <iostream>

using Row = std::vector< int >;
using Matrice = std::vector< Row >;
using RowIt = Row::iterator;

/// @struct Params
/// @var n nombre de colonnes
/// @var i nombre de 1 par colonne
/// @var j nombre de 1 par ligne
/// @var k nombre de lignes
struct Params {
    int n,i,j,k;
};

std::ostream& operator<< (std::ostream& os, const Matrice& mat);
Matrice operator* (const Matrice& lhs, const Matrice& rhs);
Matrice toFull(const Matrice& mat, int rowSize);
int randomInt(int min, int max);

void swapColumns(Matrice& mat, int r, int x);
void transvection(Row& row_r, Row& row_x);
void annulation(Matrice& mat, int r);

Matrice gauss_jordan(Matrice mat);
bool isFullRank(const Matrice& mat);

Matrice gallager(Params p);
Matrice macKayNeal(Params p);

Matrice transpose(const Matrice& mat);
void oppose(Matrice& mat);

Matrice getP(Matrice H);

#endif //LDPC_FONCTIONS_HPP
