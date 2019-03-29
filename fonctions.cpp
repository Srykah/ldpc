#include "fonctions.hpp"
#include "SparseMatrix.hpp"
#include "FullMatrix.hpp"
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

int randomInt(int min, int max) {
    return std::uniform_int_distribution<int>(min,max)(mt);
}

bool symError(float proba) {
    return std::uniform_real_distribution<float>(0.f, 1.f)(mt) <= proba;
}

int poids(const std::vector<bool>& vec) {
    int res = 0;
    for (int x : vec) {
        if (x)
            res++;
    }
    return res;
}

namespace {
    auto getPivot(SparseMatrix& mat, int r) {
        int y(r-1);
        SparseRow::const_iterator it;
        do {
            y++;
            it = std::lower_bound(mat[y].begin(), mat[y].end(), r);
        } while (y < mat.getParams().k-1 && it == mat[y].end());
        return std::pair((it == mat[y].end() ? -1 : *it), y);
    }
}

SparseMatrix gauss_jordan(SparseMatrix H) {
    for (int r(0); r < H.getParams().k; r++) {
        auto [x, y] = getPivot(H, r);
        if (x != -1) {
            if (y != r)
                H.swapRows(r, y);
            if (x != r)
                H.swapColumns(r, x);
            H.annulation(r);
        }
    }
    return H;
}

SparseMatrix gallager(Params p) {
    if (p.n/float(p.j) != p.n/p.j)
        throw std::domain_error("Gallager : invalid arguments");

    SparseMatrix result(p);

    std::vector<int> perm(p.n, 0);
    std::iota(perm.begin(), perm.end(), 0); // perm[x] = x;

    for (int r = 0; r < p.n/p.j; r++) { // r = indice de la ligne dans la sous-matrice // parallelisable
        for (int c = r*p.j; c < (r+1)*p.j; c++) { // c = indice de la colonne
            result.matrix[r].push_back(c);
        }
    }
    for (int s = 1; s < p.n/p.j; s++) { // s = indice de la sous-matrice // parallelisable
        std::shuffle(perm.begin(), perm.end(), mt);
        for (int c = 0; c < p.n; c++) {
            int r = s*(p.n/p.j) + perm[c] / p.j;
            result.matrix[r].push_back(c);
        }
    }

    return result;
}

bool helper(SparseMatrix& mat, int c, std::vector<int>& liste, const std::function<bool(int)>& pred) {
    if (c % 1000 == 0)
        std::cerr << c << ' ';
    if (c == mat.getParams().n)
        return liste.empty();
    if (liste.size() < mat.getParams().i)
        return false;

    for (int essai = 0; essai < NB_ESSAIS; essai++) { // parallelisable
        shuffle_n(liste, mat.getParams().i);
        std::vector<int> lignesChoisies(liste.begin(), liste.begin()+mat.getParams().i);
        for (int x : lignesChoisies) // parallelisable
            mat.matrix[x].push_back(c);
        std::vector<int> lignesPleines(mat.getParams().i, 0);
        lignesPleines.erase(std::copy_if(liste.begin(), liste.end(), lignesPleines.begin(), pred), lignesPleines.end());
        liste.erase(std::remove_if(liste.begin(), liste.end(), pred), liste.end());

        if (helper(mat, c+1, liste, pred))
            return true;

        for (int x : lignesChoisies) // parallelisable
            mat.matrix[x].pop_back();
        liste.insert(liste.end(), lignesPleines.begin(), lignesPleines.end());
    }
    return false;
}

SparseMatrix macKayNeal(Params p) {
    SparseMatrix mat(p);

    std::vector<int> liste(p.k, 0);
    std::iota(liste.begin(), liste.end(), 0); // liste[x] = x;

    auto pred = [&p, &mat](int x){ return mat[x].size() >= p.j; };

    helper(mat, 0, liste, pred);
    return mat;
}

FullMatrix getP(SparseMatrix G) {
    for (auto& row : G.matrix) { // parallelisable
        row.erase(row.begin());
        for (auto& elem : row) // parallelisable
            elem -= G.p.k;
    }
    G.p.n -= G.p.k;
    return G.toFull().transpose();
}
