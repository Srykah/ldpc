//
// Created by Thomas on 15/03/2019.
//

#ifndef LDPC_DECODE_HPP
#define LDPC_DECODE_HPP

#include "SparseMatrix.hpp"
#include <vector>

void canalSym(std::vector<bool>& code, float proba);
void BIAWGNC(std::vector<bool>& code, float sigma);
std::vector<float> LLR_Sym(const std::vector<bool>& code, float proba);
std::vector<float> LLR_Gauss(const std::vector<bool>& code, float sigma);
std::vector<bool> algoAGallager(std::vector<bool> code, const SparseMatrix& H);
std::vector<bool> SPA(std::vector<float> gamma, const SparseMatrix& H);

#endif //LDPC_DECODE_HPP
