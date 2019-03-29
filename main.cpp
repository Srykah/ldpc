#include "fonctions.hpp"
#include "SparseMatrix.hpp"
#include "FullMatrix.hpp"
#include "decode.hpp"
#include <iostream>
#include <algorithm>
#include <vector>

int main()
{
    /*Params p(12, 4, 9);
    std::vector< std::vector<int> > exempleMat({
        {0,5,7,9},
        {0,3,4,10},
        {1,4,6,8},
        {2,5,10,11},
        {2,6,7,11},
        {1,4,8,10},
        {0,3,6,9},
        {1,5,7,9},
        {2,3,8,11}
    });
    SparseMatrix H(p, std::move(exempleMat));*/

    Params p(1000, 4, 750);
    float proba = 0.05;
    SparseMatrix H = macKayNeal(p);
    SparseMatrix G = gauss_jordan(H);
    //std::cerr << "H :\n" << H << std::endl;
    std::cerr << "H est valide ? " << H.estValide() << std::endl;
    //std::cerr << "G :\n" << G << std::endl;
    bool fullRank = G.isFullRank();
    std::cerr << "H est full-rank ? " << fullRank << std::endl;

    if (fullRank) {
        FullMatrix P = getP(G);
        //std::cerr << "P:\n" << P << std::endl;

        std::vector<bool> message(p.n-p.k, false);
        //std::generate(message.begin(), message.end(), [](){ return symError(0.5); });
        //std::cerr << "message :\t" << message << std::endl;

        std::vector<bool> code = (FullMatrix(message, true) * P)[0];
        code.insert(code.end(), message.begin(), message.end());
        //std::cerr << "code :\t" << code << std::endl;

        canalSym(code, proba);
        std::cerr << "poids du code bruite :\t" << poids(code) << std::endl;

        std::vector<float> gamma = LLR_Sym(code, proba);
        //std::cerr << "gamma :\t" << gamma << std::endl;

        std::vector<bool> resultat = SPA(code, H, gamma);
        std::cerr << "poids du resultat :\t" << poids(resultat) << std::endl;
    }

    return 0;
}
