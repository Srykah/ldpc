#ifndef LDPC_FONCTIONS_HPP
#define LDPC_FONCTIONS_HPP

#include <vector>
#include <iostream>

using Row = std::vector< int >;
using Matrice = std::vector< Row >;
using RowIt = Row::iterator;

std::ostream& operator<< (std::ostream& os, const Matrice& mat);

void swapColumns(Matrice& mat, int r, int x);
void transvection(Row& row_r, Row& row_x);
void annulation(Matrice& mat, int r);

void gauss_jordan(Matrice& mat);

#endif //LDPC_FONCTIONS_HPP
