#include "fonctions.hpp"
#include "SparseMatrix.hpp"
#include "FullMatrix.hpp"
#include "decode.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

int main()
{
    Params p(12,4,9);
    float proba = 0.2;
    SparseMatrix H = macKayNeal(p);
    std::cout << "H:" << H << std::endl;
    SparseMatrix G = gauss_jordan(H);
    std::cout << "G:" << G << std::endl;
    if (G.isFullRank()) {
        FullMatrix P = getP(G);
        std::cout << "P:" << P << std::endl;
        std::vector<bool> code(p.n, false);
        while (poids(code) == 0) {
            canalSym(code, proba);
        }
        std::cout << "code:" << code;
        std::vector<bool> resultatGallager = algoAGallager(code, H);
        std::cout << "\tGallager:" << resultatGallager;
        std::vector<float> gamma = LLR_Sym(code, proba);
        std::vector<bool> resultatSPA = SPA(gamma, H);
        std::cout << "\tSPA:" << resultatSPA;
    }
    return 0;
}
