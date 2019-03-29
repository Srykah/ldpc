//
// Created by Thomas on 08/03/2019.
//

#ifndef LDPC_PARAMS_HPP
#define LDPC_PARAMS_HPP

/// @struct Params
/// @var n nombre de colonnes (largeur)
/// @var i nombre de 1 par colonne
/// @var j nombre de 1 par ligne
/// @var k nombre de lignes (hauteur)
struct Params {
    Params(int n, int k) : n(n), i(0), j(0), k(k) {}
    Params(int n, int j, int k) : n(n), i(k*j/n), j(j), k(k) {}
    Params(int n, int k, float density) : n(n), i(density*k), j(density*n), k(k) {}
    int n,i,j,k;
};

#endif //LDPC_PARAMS_HPP
