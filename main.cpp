#include "fonctions.hpp"
#include <iostream>
#include <algorithm>

/** TODO :
 * Faire la diff entre MatricePleine et MatriceCreuse
 * Checker les dimensions dans le produit
 * Passer en POO
**/

int main()
{
    Params p {12, 3, 4, 9};
    Matrice H;
    Matrice gj;
    do {
        H = macKayNeal(p);
        gj = gauss_jordan(H);
    } while (!isFullRank(gj));
    std::cout << "H:\n" << H << std::endl;

    Matrice P = getP(gj);
    std::cout << "P:\n" << P << std::endl;

    Matrice u {Row(p.n-p.k)};
    std::generate(u[0].begin(), u[0].end(), [](){ return randomInt(0, 1); });
    std::cout << "u :\t" << u << std::endl;

    Matrice code = u * P;
    code[0].insert(code[0].end(), u[0].begin(), u[0].end());
    std::cout << "code :\t" << code << std::endl;

    return 0;
}
