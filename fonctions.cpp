#include "fonctions.hpp"
#include <algorithm>

#include <random>
#include <iostream>

std::ostream& operator<< (std::ostream& os, const Matrice& mat) {
    for (const Row& row : mat) {
        for (int x : row) {
            os << x << '\t';
        }
        os << '\n';
    }

    return os;
}

namespace {
    std::random_device rd;
    std::mt19937 mt(rd());
    int randint(int min, int max) {
        return std::uniform_int_distribution<>(min, max)(mt);
    }
    auto pred = [](auto& elem){ return true; };
}

void swapColumns(Matrice& mat, int r, int x) {
    for (Row& row : mat) {
        RowIt itr = std::lower_bound(row.begin(), row.end(), r);
        RowIt itx = std::lower_bound(row.begin(), row.end(), x);
        bool r_exists(itr != row.end() && *itr == r);
        bool x_exists(itx != row.end() && *itx == x);
        if (x_exists && !r_exists) {
            row.erase(itx);
            row.insert(itr, r);
        } else if (r_exists && !x_exists) {
            int distr = std::distance(row.begin(), itr); // reallocation
            row.insert(itx, x);
            row.erase(row.begin()+distr);
        }
    }
}

void transvection(Row& row_r, Row& row_x) {
    int ir = row_r.size()-1;
    int ix = row_x.size()-1;
    while (ir >= 0 && ix >= 0) {
        int valr = row_r[ir], valx = row_x[ix];
        if (valr == valx) {
            row_x.erase(row_x.begin()+ix);
            ix--;
            ir--;
        } else if (valr < valx) {
            ix--;
        } else { // valr > valx
            row_x.insert(row_x.begin()+ix+1, valr);
            ir--;
        }
    }
}

void annulation(Matrice& mat, int r) {
    for (int y(r+1); y < mat.size(); y++) {
        if (std::find(mat[y].begin(), mat[y].end(), r) != mat[y].end()) {
            transvection(mat[r], mat[y]);
        }
    }
}

void gauss_jordan(Matrice& mat) {
    for (int r(0); r < mat.size(); r++) {
        std::cout << "r = " << r << std::endl;
        RowIt it = std::find_if(mat[r].begin()+r, mat[r].end(), pred);
        int y(r);
        while (y < mat.size()-1 && it == mat[y].end()) {
            y++;
            it = std::find_if(mat[y].begin()+r, mat[y].end(), pred);
        }
        if (it != mat[y].end()) {
            int x = *it;
            std::cout << "\t(x,y) = (" << x << ',' << y << ')' << std::endl;
            if (y != r) {
                std::cout << "\tswapping lines " << r << " and " << y << std::endl;
                std::swap(mat[y], mat[r]);
            }
            if (x != r) {
                std::cout << "\tswapping columns " << r << " and " << x << std::endl;
                swapColumns(mat, r, x);
            }
            annulation(mat, r);
        }
    }
}
