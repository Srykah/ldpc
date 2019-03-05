#include "fonctions.hpp"
#include <algorithm>
#include <random>
#include <iostream>
#include <stdexcept>
#include <chrono>

namespace {
    const size_t NB_ESSAIS = 10;
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
    for (int x = 0; x < largeur; x++) {
        for (int y = 0; y < hauteur; y++) {
            for (int z = 0; z < rhs.size(); z++) {
                result[y][x] += lhs[y][z] * rhs[z][x];
                result[y][x] %= 2;
            }
        }
    }
    return result;
}

Matrice toFull(const Matrice& mat, int rowSize) {
    Matrice result(mat.size(), Row(rowSize, 0));
    for (int i(0); i < mat.size(); i++) {
        const Row& row = mat[i];
        for (int val : row)
            result[i][val] = 1;
    }

    return result;
}

int randomInt(int min, int max) {
    return std::uniform_int_distribution<int>(min,max)(mt);
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
    for (int y(0); y < mat.size(); y++) {
        if (y != r && std::find(mat[y].begin(), mat[y].end(), r) != mat[y].end()) {
            transvection(mat[r], mat[y]);
        }
    }
}

Matrice gauss_jordan(Matrice mat) {
    for (int r(0); r < mat.size(); r++) {
        int y(r);
        RowIt it = std::find_if(mat[y].begin(), mat[y].end(), [r](int x){ return x >= r; });
        while (y < mat.size()-1 && it == mat[y].end()) {
            y++;
            it = std::find_if(mat[y].begin(), mat[y].end(), [r](int x){ return x >= r; });
        }
        if (it != mat[y].end()) {
            int x = *it;
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
    for (int i = 0; i < hauteur; i++) {
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
    bool helper(Matrice& mat, const Params& p, int c, Row& liste) {
        if (c == p.n)
            return liste.empty();
        if (liste.size() < p.i)
            return false;

        for (int _ = 0; _ < NB_ESSAIS; _++) {
            std::shuffle(liste.begin(), liste.end(), mt);
            Row debut(liste.begin(), liste.begin()+p.i);
            for (int x : debut)
                mat[x].push_back(c);
            std::sort(liste.begin(), liste.end(), [&p, &mat](int x, int y){
                return mat[x].size() < p.j;
            });
            RowIt it = std::find_if(liste.begin(), liste.end(), [&p, &mat](int x){ return mat[x].size() >= p.j; });
            Row aSupprimer(it, liste.end());
            liste.erase(it, liste.end());
            if (helper(mat, p, c+1, liste))
                return true;
            // else
            for (int x : debut)
                mat[x].pop_back();
            liste.insert(liste.end(), aSupprimer.begin(), aSupprimer.end());
        }
        return false;
    }
}

Matrice macKayNeal(Params p) {
    Matrice mat;
    mat.resize(p.k);

    std::vector<int> liste(p.k, 0);
    std::iota(liste.begin(), liste.end(), 0); // liste[x] = x;

    helper(mat, p, 0, liste);
    return mat;
}

Matrice transpose(const Matrice& mat) {
    unsigned long long int hauteur = mat.size(), largeur = mat[0].size();
    Matrice result(largeur, Row(hauteur, 0));
    for (int y = 0; y < hauteur; y++) {
        for (int  x = 0; x < largeur; x++) {
            result[x][y] = mat[y][x];
        }
    }
    return result;
}

void oppose(Matrice& mat) {
    unsigned long long int hauteur = mat.size(), largeur = mat[0].size();
    for (int y = 0; y < hauteur; y++) {
        for (int  x = 0; x < largeur; x++) {
            mat[y][x] = !mat[y][x];
        }
    }
}

Matrice getP(Matrice H) {
    size_t largeur = 0;
    for (Row& row : H) {
        int tmp = *std::max_element(row.begin(), row.end());
        if (tmp > largeur)
            largeur = tmp;
    }
    largeur++;
    size_t hauteur = H.size();
    gauss_jordan(H);
    for (auto& row : H) {
        row.erase(row.begin());
        for (auto& elem : row)
            elem -= hauteur;
    }
    Matrice P = transpose(toFull(H, largeur - hauteur));
    oppose(P);
    return P;
}
