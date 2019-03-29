#ifndef LDPC_FONCTIONS_HPP
#define LDPC_FONCTIONS_HPP

#include "Params.hpp"
#include <vector>
#include <iostream>

class SparseMatrix;
class FullMatrix;

int randomInt(int min, int max);
template <typename T>
void shuffle_n(std::vector<T>& vec, int n);
bool symError(float proba);
int poids(const std::vector<bool>& vec);

SparseMatrix gauss_jordan(SparseMatrix H);

SparseMatrix gallager(Params p);
SparseMatrix macKayNeal(Params p);

FullMatrix getP(SparseMatrix G);

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> vec);

#include "fonctions.inl"

#endif //LDPC_FONCTIONS_HPP
