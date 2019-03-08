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
    Params p(12, 8, 9);
    /*
    Matrice H;
    Matrice G;
    do {
        H = macKayNeal(p);
        G = gauss_jordan(H);
    } while (!isFullRank(G));*/
    Matrice H {
        {0,1,2,3},
        {4,5,6,7},
        {8,9,10,11},
        {0,3,6,9},
        {1,4,7,10},
        {2,5,8,11},
        {0,1,6,7},
        {2,3,8,9},
        {4,5,10,11}
    };
    Matrice G = gauss_jordan(H);
    std::cerr << "H :\n" << H << std::endl;
    std::cerr << "H est valide ? " << estValide(H, p) << std::endl;
    std::cerr << "G :\n" << G << std::endl;
    std::cerr << "H est full-rank ? " << isFullRank(G) << std::endl;

    Matrice P = getP(G, p);
    std::cerr << "P:\n" << P << std::endl;

    Matrice u {Row(p.n-p.k)};
    std::generate(u[0].begin(), u[0].end(), [](){ return randomInt(0, 1); });
    std::cerr << "u :\t" << u << std::endl;

    Matrice code = u * P;
    code[0].insert(code[0].end(), u[0].begin(), u[0].end());
    std::cerr << "code :\t" << code << std::endl;

    return 0;
}
