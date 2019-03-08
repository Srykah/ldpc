#include "fonctions.hpp"
#include <algorithm>
#include <random>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <functional>

namespace {
    const size_t NB_ESSAIS = 5;
    std::mt19937 mt(std::chrono::system_clock::now().time_since_epoch().count());
}

std::ostream& operator<< (std::ostream& os, const Matrice& mat) {
    for (const Row& row : mat) {
        os << "{";
        for (int x : row) {
            os << x << '\t';
        }
        os << "}\n";
    }

    return os;
}

Matrice operator* (const Matrice& lhs, const Matrice& rhs) {
    size_t hauteur = lhs.size(), largeur = rhs[0].size();
    Matrice result(hauteur, Row(largeur, 0));
    for (int x = 0; x < largeur; x++) { // parallelisable
        for (int y = 0; y < hauteur; y++) { // parallelisable
            for (int z = 0; z < rhs.size(); z++) {
                result[y][x] += lhs[y][z] * rhs[z][x];
            }
            result[y][x] %= 2;
        }
    }
    return result;
}

Matrice toFull(const Matrice& mat, int rowSize) {
    Matrice result(mat.size(), Row(rowSize, 0));
    for (int i(0); i < mat.size(); i++) { // parallelisable
        const Row& row = mat[i];
        for (int val : row) // parallelisable
            result[i][val] = 1;
    }

    return result;
}

int randomInt(int min, int max) {
    return std::uniform_int_distribution<int>(min,max)(mt);
}

void shuffle_n(Row& row, int n) {
    int i = 0;
    while (n--) {
        std::swap(row[i], row[randomInt(i, row.size()-1)]);
        i++;
    }
}

bool estValide(const Matrice& mat, Params p) {
    if (mat.size() != p.k)
        return false;
    Row count(p.n, 0);
    for (const Row& row : mat) {
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

void swapColumns(Matrice& mat, int r, int x) {
    for (Row& row : mat) { // parallelisable
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
    for (int y(0); y < mat.size(); y++) { // parallelisable
        if (y != r && std::find(mat[y].begin(), mat[y].end(), r) != mat[y].end()) {
            transvection(mat[r], mat[y]);
        }
    }
}

namespace {
    auto getPivot(Matrice& mat, int r) {
        int y(r-1);
        RowIt it;
        do {
            y++;
            it = std::lower_bound(mat[y].begin(), mat[y].end(), r);
        } while (y < mat.size()-1 && it == mat[y].end());
        return std::pair((it == mat[y].end() ? -1 : *it), y);
    }
}

Matrice gauss_jordan(Matrice mat) {
    for (int r(0); r < mat.size(); r++) {
        auto [x, y] = getPivot(mat, r);
        if (x != -1) {
            if (y != r)
                std::swap(mat[y], mat[r]);
            if (x != r)
                swapColumns(mat, r, x);
            annulation(mat, r);
        }
    }
    return mat;
}

bool isFullRank(const Matrice& mat) {
    size_t hauteur = mat.size();
    for (int i = 0; i < hauteur; i++) { // parallelisable
        if (mat[i].empty() || mat[i][0] != i ||(mat[i].size() > 1 && mat[i][1] < hauteur))
            return false;
    }
    return true;
}

Matrice gallager(Params p) {
    if (p.n/float(p.j) != p.n/p.j)
        throw std::domain_error("Gallager : invalid arguments");

    Matrice matrice;
    matrice.resize(p.k);

    std::vector<int> perm(p.n, 0);
    std::iota(perm.begin(), perm.end(), 0); // perm[x] = x;

    for (int r = 0; r < p.n/p.j; r++) { // w = indice de la ligne dans la sous-matrice
        for (int c = r*p.j; c < (r+1)*p.j; c++) {// s = indice de la colonne
            matrice[r].push_back(c);
        }
    }
    for (int s = 1; s < p.n/p.j; s++) { // s = indice de la sous-matrice
        std::shuffle(perm.begin(), perm.end(), mt);
        for (int c = 0; c < p.n; c++) {
            int r = s*(p.n/p.j) + int(perm[c] / p.j);
            matrice[r].push_back(c);
        }
    }

    return matrice;
}

namespace {
    bool helper(Matrice& mat, const Params& p, int c, Row& liste, const std::function<bool(int)>& pred) {
        if (c == p.n)
            return liste.empty();
        if (liste.size() < p.i)
            return false;

        for (int essai = 0; essai < NB_ESSAIS; essai++) { // parallelisable
            shuffle_n(liste, p.i);
            Row lignesChoisies(liste.begin(), liste.begin()+p.i);
            for (int x : lignesChoisies) // parallelisable
                mat[x].push_back(c);
            Row lignesPleines(p.i, 0);
            lignesPleines.erase(std::copy_if(liste.begin(), liste.end(), lignesPleines.begin(), pred), lignesPleines.end());
            liste.erase(std::remove_if(liste.begin(), liste.end(), pred), liste.end());

            if (helper(mat, p, c+1, liste, pred))
                return true;

            for (int x : lignesChoisies) // parallelisable
                mat[x].pop_back();
            liste.insert(liste.end(), lignesPleines.begin(), lignesPleines.end());
        }
        return false;
    }
}

Matrice macKayNeal(Params p) {
    Matrice mat;
    mat.resize(p.k);

    std::vector<int> liste(p.k, 0);
    std::iota(liste.begin(), liste.end(), 0); // liste[x] = x;

    auto pred = [&p, &mat](int x){ return mat[x].size() >= p.j; };

    helper(mat, p, 0, liste, pred);
    return mat;
}

Matrice transpose(const Matrice& mat) {
    unsigned long long int hauteur = mat.size(), largeur = mat[0].size();
    Matrice result(largeur, Row(hauteur, 0));
    for (int y = 0; y < hauteur; y++) { // parallelisable
        for (int  x = 0; x < largeur; x++) { // parallelisable
            result[x][y] = mat[y][x];
        }
    }
    return result;
}

void comp2(Matrice& mat) {
    unsigned long long int hauteur = mat.size(), largeur = mat[0].size();
    for (int y = 0; y < hauteur; y++) { // parallelisable
        for (int  x = 0; x < largeur; x++) { // parallelisable
            mat[y][x] = !mat[y][x];
        }
    }
}

Matrice getP(Matrice H, Params p) {
    gauss_jordan(H);
    for (auto& row : H) { // parallelisable
        row.erase(row.begin());
        for (auto& elem : row) // parallelisable
            elem -= p.k;
    }
    return transpose(toFull(H, p.n - p.k));
}
