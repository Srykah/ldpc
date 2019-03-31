//
// Created by Thomas on 15/03/2019.
//

#include "decode.hpp"
#include "FullMatrix.hpp"
#include <cmath>
#include <algorithm>

const int NB_ITER = 50;

void canalSym(std::vector<bool>& code, float proba) {
    for (auto elem : code) {
        if (symError(proba))
            elem = !elem;
    }
}

void BIAWGNC(std::vector<bool>& code, float sigma) {
    float proba = 0.5*std::erfc(1.f/(sigma*std::sqrt(2)));
    canalSym(code, proba);
}

std::vector<float> LLR_Sym(const std::vector<bool>& code, float proba) {
    std::vector<float> result(code.size());
    float llr = std::log(proba/(1-proba));
    for (int i = 0; i < code.size(); i++) {
        result[i] = (code[i]?1.f:-1.f)*llr;
    }
    return result;
}

std::vector<float> LLR_Gauss(const std::vector<bool>& code, float sigma) {
    return LLR_Sym(code, 0.5*std::erfc(1.f/(sigma*std::sqrt(2))));
}

std::vector<bool> algoAGallager(std::vector<bool> code, const SparseMatrix& H) {
    auto p = H.getParams();
    auto Ht = H;
    Ht.transpose();
    FullMatrix VtoC(p.j, p.n);
    FullMatrix CtoV(p.i, p.k);

    for (int v = 0; v < p.n; v++) { // pour chaque variable
        for (int c_v = 0; c_v < p.i; c_v++) { // pour chaque check attachée à cette variable
            VtoC[v][c_v] = code[v]; // le message initial de la variable à la check est le bit du code associé à la variable
        }
    }

    for (int iter = 0; iter < NB_ITER; iter++) {
        CtoV = FullMatrix(p.i, p.k);
        for (int v = 0; v < p.n; v++) { // pour chaque variable
            for (int c_v = 0; c_v < p.i; c_v++) { // pour chaque check attachée à cette variable
                int c = Ht[v][c_v];
                int v_c0 = std::distance(H[c].begin(), std::lower_bound(H[c].begin(), H[c].end(), v));
                for (int v_c = 0; v_c < p.j; v_c++) { // pour chaque variable attachée à cette check...
                    if (v_c != v_c0) // ...exceptée la variable envoyant le message
                        CtoV[c][v_c] = CtoV[c][v_c] ^ VtoC[v][c_v]; // on propage le message
                }
            }
        }

        for (int v = 0; v < p.n; v++) {
            for (int c_v0 = 0; c_v0 < p.i; c_v0++) {
                bool inverse = !code[v];
                bool allInverse = true;
                for (int c_v = 0; c_v < p.i; c_v++) {
                    if (c_v != c_v0) {
                        int c = Ht[v][c_v];
                        int v_c = std::distance(H[c].begin(), std::lower_bound(H[c].begin(), H[c].end(), v));
                        if (CtoV[c][v_c] != inverse) {
                            allInverse = false;
                            break;
                        }
                    }
                }
                if (allInverse) {
                    VtoC[v][c_v0] = inverse;
                } else {
                    VtoC[v][c_v0] = code[v];
                }
            }
        }
    }

    for (int v = 0; v < p.n; v++) {
        int count = 0;
        for (int c_v = 0; c_v < p.i; c_v++) {
            count += VtoC[v][c_v]; // on compte le nombre de messages à 1 envoyés depuis la variable
        }
        code[v] = (count > 0.5 * p.i); // le seuil de décision est la moitié des CN par VN.
    }

    std::cerr << "Code apres correction : " << code << std::endl;
    return std::vector<bool>(code.begin()+p.k, code.end());
}

std::vector<bool> SPA(std::vector<bool> code, const SparseMatrix& H, std::vector<float> gamma) {
    auto p = H.getParams();
    auto Ht = H;
    Ht.transpose();
    // preliminary test : check syndrome
    for (int j = 0; j < p.n; j++) {
        code[j] = gamma[j] < 0;
    }
    auto syndrome = code * Ht;
    //std::cerr << "Test preliminaire : syndrome = " << syndrome << std::endl;
    if (syndrome == std::vector<bool>(p.k, false))
        return code;

    // initialisation
    std::vector<std::vector<float>> phi(p.n, std::vector<float>(p.i));
    std::vector<std::vector<float>> psi(p.k, std::vector<float>(p.j));
    for (int j = 0; j < p.n; j++) {
        for (int m = 0; m < p.i; m++) {
            phi[j][m] = gamma[j];
        }
    }

    for (int essai = 0; essai < NB_ITER; essai++) {
        // step 1: compute messages from CN to VN (in parallel)
        for (int i = 0; i < p.k; i++) {
            float prod = 1.f;
            for (int m = 0; m < p.j; m++) {
                prod *= std::tanh(0.5*phi[m][i]);
            }
            for (int m = 0; m < p.j; m++) {
                psi[i][m] = 2*std::atanh(prod/std::tanh(0.5*phi[m][i]));
            }
        }
        // stop criterion : check syndrome
        // step : compute beliefs
        std::vector<float> lambda = gamma;
        for (int j = 0; j < p.n; j++) {
            for (int u = 0; u < p.i; u++) {
                lambda[j] += psi[Ht[j][u]][u];
            }
            code[j] = lambda[j] < 0;
        }
        syndrome = code * Ht;
        //std::cerr << "Boucle principale, essai " << essai << " : syndrome = " << syndrome << std::endl;
        if (syndrome == std::vector<bool>(p.k, false))
            return code;
        // step 2: compute messages from VN to CN (in parallel)
        for (int j = 0; j < p.n; j++) {
            for (int u = 0; u < p.i; u++) {
                phi[j][u] = lambda[j] - phi[Ht[j][u]][j];
            }
        }
    }
    return code;
}