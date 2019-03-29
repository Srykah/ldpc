//
// Created by Thomas on 08/03/2019.
//

#ifndef LDPC_SPARSEMATRIX_HPP
#define LDPC_SPARSEMATRIX_HPP

#include "Params.hpp"
#include "fonctions.hpp"
#include <vector>
#include <ostream>
#include <functional>

class FullMatrix;

using SparseRow = std::vector< int >;

class SparseMatrix {
public:
    explicit SparseMatrix(Params p);
    explicit SparseMatrix(Params p, std::vector< SparseRow >&& init);
    SparseMatrix(const SparseMatrix& other) = default;
    SparseMatrix& operator=(const SparseMatrix& other) = default;
    SparseMatrix(SparseMatrix&& other) = default;
    SparseMatrix& operator=(SparseMatrix&& other) = default;

    void swapRows(int r1, int r2);
    void swapColumns(int c1, int c2);
    /// @brief M[y] += M[r]
    /// \param r index of row to add
    /// \param y index of row to modifiy
    void transvection(int r, int y);
    void annulation(int r);
    FullMatrix toFull() const;

    SparseMatrix& transpose();

    const SparseRow& operator[](size_t row) const;
    Params getParams() const;
    bool estValide() const;
    bool isFullRank() const;

private:
    Params p;
    std::vector< SparseRow > matrix;

    friend std::ostream& operator<< (std::ostream& os, const SparseMatrix& mat);
    friend SparseMatrix gallager(Params p);
    friend bool helper(SparseMatrix& mat, int c, std::vector<int>& liste, const std::function<bool(int)>& pred);
    friend FullMatrix getP(SparseMatrix H);
};

std::ostream& operator<< (std::ostream& os, const SparseMatrix& mat);
std::vector<bool> operator*(const std::vector<bool>& row, const SparseMatrix& mat);

#endif //LDPC_SPARSEMATRIX_HPP
