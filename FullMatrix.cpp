//
// Created by Thomas on 08/03/2019.
//

#include "FullMatrix.hpp"

FullMatrix::FullMatrix(size_t height, size_t width)
: matrix(height, FullRow(width, false)) {}

FullMatrix::FullMatrix(const FullRow &vec, bool horizontal) {
    if (horizontal) {
        matrix.resize(1);
        matrix[0] = vec;
    } else {
        matrix.resize(vec.size());
        for (int i = 0; i < vec.size(); i++)
            matrix[i].push_back(vec[i]);
    }
}

FullMatrix::FullMatrix(std::vector<FullRow>&& init)
: matrix(init) {}

FullRow& FullMatrix::operator[](size_t row) {
    return matrix[row];
}

const FullRow& FullMatrix::operator[](size_t row) const {
    return matrix[row];
}

size_t FullMatrix::width() const {
    return (matrix.empty() ? 0 : matrix[0].size());
}

size_t FullMatrix::height() const {
    return matrix.size();
}

FullMatrix& FullMatrix::transpose() {
    const size_t initialWidth(width()), initialHeight(height());
    const size_t finalWidth(initialHeight), finalHeight(initialWidth);
    if (initialWidth >= initialHeight) {
        matrix.resize(finalHeight, FullRow(finalWidth));
        for (int y = 0; y < initialHeight; y++) { // parallelisable
            for (int x = y+1; x < initialWidth; x++) { // parallelisable
                std::swap(matrix[x][y], matrix[y][x]);
            }
            matrix[y].erase(matrix[y].begin()+finalWidth, matrix[y].end());
        }
    } else {
        for (size_t y = 0; y < finalHeight; y++) {
            matrix[y].resize(finalWidth);
        }
        for (int x = 0; x < initialWidth; x++) { // parallelisable
            for (int y = x+1; y < initialHeight; y++) { // parallelisable
                std::swap(matrix[x][y], matrix[y][x]);
            }
        }
        matrix.resize(finalHeight);
    }
    return *this;
}

FullMatrix& FullMatrix::comp2() {
    for (int y = 0; y < height(); y++) { // parallelisable
        for (int  x = 0; x < width(); x++) { // parallelisable
            matrix[y][x] = !matrix[y][x];
        }
    }

    return *this;
}

FullMatrix operator* (const FullMatrix& lhs, const FullMatrix& rhs) {
    FullMatrix result(lhs.height(), rhs.width());
    for (int x = 0; x < rhs.width(); x++) { // parallelisable
        for (int y = 0; y < lhs.height(); y++) { // parallelisable
            for (int z = 0; z < lhs.width(); z++) { // more or less parallelisable
                result[y][x] = result[y][x] ^ (lhs[y][z] && rhs[z][x]);
            }
        }
    }
    return result;
}

std::ostream& operator<< (std::ostream& os, const FullMatrix& mat) {
    for (const FullRow& row : mat.matrix) {
        os << "{ ";
        for (int x : row) {
            os << ' ' << x;
        }
        os << " }\n";
    }

    return os;
}
