//
// Created by Thomas on 08/03/2019.
//

#ifndef LDPC_FULLMATRIX_HPP
#define LDPC_FULLMATRIX_HPP

#include <iostream>
#include <vector>

using FullRow = std::vector<bool>;

class FullMatrix {
public:
    FullMatrix(size_t height, size_t width);
    FullMatrix(const FullRow& vec, bool horizontal);
    explicit FullMatrix(std::vector<FullRow>&& init);

    FullRow& operator[](size_t row);
    const FullRow& operator[](size_t row) const;
    size_t width() const;
    size_t height() const;

    FullMatrix& transpose();
    FullMatrix& comp2();

private:
    std::vector<FullRow> matrix;

    friend FullMatrix operator* (const FullMatrix& lhs, const FullMatrix& rhs);
    friend std::ostream& operator<< (std::ostream& os, const FullMatrix& mat);
};

std::ostream& operator<< (std::ostream& os, const FullMatrix& mat);

#endif //LDPC_FULLMATRIX_HPP
