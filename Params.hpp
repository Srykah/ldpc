//
// Created by Thomas on 08/03/2019.
//

#ifndef LDPC_PARAMS_HPP
#define LDPC_PARAMS_HPP

/// @struct Params
/// @var n nombre de colonnes (largeur) (nombre de VN)
/// @var i nombre de 1 par colonne (nombre de CN par VN)
/// @var j nombre de 1 par ligne (nombre de VN par CN)
/// @var k nombre de lignes (hauteur) (nombre de CN)
struct Params {
    Params(int n, int k) : n(n), i(0), j(0), k(k) {}
    Params(int n, int j, int k) : n(n), i(k*j/n), j(j), k(k) {}
    Params(int n, int k, float density) : n(n), i(density*k), j(density*n), k(k) {}
    int n,i,j,k;
};

#endif //LDPC_PARAMS_HPP
