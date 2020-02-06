//
// Created by Thomas on 08/03/2019.
//

#include "SparseMatrix.hpp"
#include "FullMatrix.hpp"
#include <algorithm>

SparseMatrix::SparseMatrix(Params p)
: p(p)
, matrix(p.k) {}

SparseMatrix::SparseMatrix(Params p, std::vector< SparseRow >&& init)
: p(p)
, matrix(init) {

}

void SparseMatrix::swapRows(int r1, int r2) {
    matrix[r1].swap(matrix[r2]);
}

void SparseMatrix::swapColumns(int c1, int c2) {
    for (SparseRow& row : matrix) { // parallelisable
        auto itr = std::lower_bound(row.begin(), row.end(), c1);
        auto itx = std::lower_bound(row.begin(), row.end(), c2);
        bool r_exists(itr != row.end() && *itr == c1);
        bool x_exists(itx != row.end() && *itx == c2);
        if (x_exists && !r_exists) {
            row.erase(itx);
            row.insert(itr, c1);
        } else if (r_exists && !x_exists) {
            int distr = std::distance(row.begin(), itr); // reallocation
            row.insert(itx, c2);
            row.erase(row.begin()+distr);
        }
    }
}

FullMatrix SparseMatrix::toFull() const {
    FullMatrix result(p.k, p.n);
    for (int i(0); i < p.k; i++) { // parallelisable
        for (int val : matrix[i]) // parallelisable
            result[i][val] = true;
    }
    return result;
}

SparseMatrix& SparseMatrix::transpose() {
    SparseMatrix next(Params(p.k, p.i, p.n));
    for (int r = 0; r < p.k; r++) {
        for (int col : matrix[r]) {
            next.matrix[col].push_back(r);
        }
    }
    matrix = std::move(next.matrix);
    return *this;
}

Params SparseMatrix::getParams() const {
    return p;
}

bool SparseMatrix::estValide() const {
    std::vector<int> count(p.n, 0);
    for (const SparseRow& row : matrix) {
        if (row.size() != p.j)
            return false;
        for (int c : row) {
            if (c >= p.n)
                return false;
            count[c]++;
        }
    }
    for (int c : count) {
        if (c != p.i)
            return false;
    }
    return true;
}

std::ostream& operator<< (std::ostream& os, const SparseMatrix& mat) {
    for (const SparseRow& row : mat.matrix) {
        os << "{";
        for (int x : row) {
            os << '\t' << x;
        }
        os << "}\n";
    }

    return os;
}

std::vector<bool> operator*(const std::vector<bool>& row, const SparseMatrix& mat) {
    std::vector<bool> result(mat.getParams().n, false);
    for (int r = 0; r < row.size(); r++) {
        if (row[r]) {
            for (int indice : mat[r])
                result[r] = !result[r];
        }
    }
    return result;
}

void SparseMatrix::transvection(int r, int y) {
    int ir = matrix[r].size()-1;
    int iy = matrix[y].size()-1;
    while (ir >= 0 && iy >= 0) {
        int valr = matrix[r][ir], valx = matrix[y][iy];
        if (valr == valx) {
            matrix[y].erase(matrix[y].begin()+iy);
            iy--;
            ir--;
        } else if (valr < valx) {
            iy--;
        } else { // valr > valx
            matrix[y].insert(matrix[y].begin()+iy+1, valr);
            ir--;
        }
    }
}

void SparseMatrix::annulation(int r) {
    for (int y(0); y < p.k; y++) { // parallelisable
        if (y != r && std::find(matrix[y].begin(), matrix[y].end(), r) != matrix[y].end()) {
            transvection(r, y);
        }
    }
}

const SparseRow& SparseMatrix::operator[](size_t row) const {
    return matrix[row];
}

bool SparseMatrix::isFullRank() const {
    for (int i = 0; i < p.k; i++) { // parallelisable
        if (matrix[i].empty() || matrix[i][0] != i ||(matrix[i].size() > 1 && matrix[i][1] < p.k))
            return false;
    }
    return true;
}
