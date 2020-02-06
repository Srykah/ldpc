#include "fonctions.hpp"
#include "SparseMatrix.hpp"
#include "FullMatrix.hpp"
#include "decode.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

#define NB_MAT 5
#define NB_CODES 20

int main()
{
    std::vector<Params> params = {
      {12,4,9},
      {1000,4,750},
      {1000,12,750},
      {10000,5,6000},
      {10000,10,6000},
      {64800,4,48600}
    };
    std::vector<float> probas = {
      0.01,
      0.05,
      0.1,
      0.33,
    };
    for (auto p : params) {
      std::cout << "Params = (" << p.n << ',' << p.j << ',' << p.k << ')' << std::endl;
      std::clog << '|';
      int numFullRank = 0;
      for (int matNum = 0; matNum < NB_MAT; matNum++) {
        std::cout << "\tMatrice n°" << matNum << std::endl;
        std::clog << '#';
        SparseMatrix H = macKayNeal(p);
        SparseMatrix G = gauss_jordan(H);
        if (G.isFullRank()) {
          numFullRank++;
          FullMatrix P = getP(G);
          for (float proba : probas) {
            float spaOutputErrorRatio = 0.f;
            float gallagerOutputErrorRatio = 0;
            for (int numCodes = 0; numCodes < NB_CODES; numCodes++) {
              std::vector<bool> code(p.n, false);
              canalSym(code, proba);
              int wCode = poids(code);
              if (wCode != 0.f) {
                std::vector<float> gamma = LLR_Sym(code, proba);
                std::vector<bool> resultatSPA = SPA(gamma, H);
                std::vector<bool> resultatGallager = algoAGallager(code, H);
                spaOutputErrorRatio += float(poids(resultatSPA))/float(wCode);
                gallagerOutputErrorRatio += float(poids(resultatGallager))/float(wCode);
              }
              std::clog << '.';
            }
            spaOutputErrorRatio /= NB_CODES;
            gallagerOutputErrorRatio /= NB_CODES;
            std::cout << "\t\tProba = " << proba << ", TES SPA : " << spaOutputErrorRatio << ", TES Gallager : " << gallagerOutputErrorRatio << std::endl;
          }
        }
      }
    }

    return 0;
}
