#pragma once

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>
#include <type_traits>

#include "../mode/mode.h"
#include "buffer.h"

namespace Linear {
using pair = std::pair<size_t, size_t>;

template <typename T>
class Matrix;

namespace Determinant {
enum class Type {
    STD = 0,
    GAUSS = 1,
};

template <typename T>
T std_alg(const Linear::Matrix<T> &matrix);

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
T gauss(const Linear::Matrix<T> &matrix);
}  // namespace Determinant

template <typename T>
class Matrix final : private detail::MatrixBuffer<T> {
   private:
    size_t nRows_ = 0;
    size_t nCols_ = 0;

    using detail::MatrixBuffer<T>::size_;
    using detail::MatrixBuffer<T>::data_;

    void reverse_gauss();

   public:
    Matrix(size_t rows, size_t cols, T value = T{});

    Matrix(size_t rows = 0) : Matrix(rows, rows) {}

    template <class InputIt1, class InputIt2>
    Matrix(size_t rows, size_t cols, InputIt1 begin_, InputIt2 end_);

    Matrix(const Matrix &rhs);

    Matrix(Matrix &&rhs) noexcept;

    ~Matrix() = default;

    void swap(Matrix &rhs) noexcept;

    Matrix &operator=(const Matrix &rhs);
    Matrix &operator=(Matrix &&rhs) noexcept;
    bool operator==(const Matrix &rhs) const;
    bool operator!=(const Matrix &rhs) const;
    Matrix &operator*=(const T number);
    Matrix &operator*=(const Matrix &rhs);
    Matrix &operator+=(const Matrix &rhs);
    Matrix &operator-=(const Matrix &rhs);

    operator std::vector<T>();

    void resize(pair shape);
    void transpose();
    void negate();
    void clear();
    void diagonalize();

    int direct_gauss();

    void swap_rows(int lhs, int rhs);
    void add_rows(int src, int dest, T factor);
    void swap_cols(int lhs, int rhs);
    void add_cols(int src, int dest, T factor);
    const T *data() const noexcept;

    static Matrix zeros(int n);
    static Matrix eye(int n);

    pair shape() const;
    int size() const noexcept;
    T trace() const;
    const T &at(size_t i, size_t j) const;
    void dump(std::ostream &stream) const;

    T &at(size_t i, size_t j);

    T determinant(Determinant::Type type = Determinant::Type::STD) const;
};

template <typename T>
std::istream &operator>>(std::istream &stream, Matrix<T> &rhs);
template <typename T>
std::ostream &operator<<(std::ostream &stream, const Matrix<T> &rhs);

template <typename T>
Matrix<T> operator+(const Matrix<T> &lhs, const Matrix<T> &rhs);
template <typename T>
Matrix<T> operator-(const Matrix<T> &lhs, const Matrix<T> &rhs);
template <typename T>
Matrix<T> operator*(const T number, const Matrix<T> &rhs);
template <typename T>
Matrix<T> operator*(const Matrix<T> &lhs, const T number);
template <typename T>
Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs);
}  // namespace Linear

template <typename T>
T Linear::Determinant::std_alg(const Linear::Matrix<T> &matrix) {
    auto shape = matrix.shape();
    size_t nRows = shape.first;
    size_t nCols = shape.second;

    if (nRows != nCols)
        throw std::invalid_argument("Attempt to calculate the determinant of a non-square matrix.");

    if (nRows == 1) {
        return matrix.at(0, 0);
    } else {
        T ans{};
        for (size_t i = 0; i < nCols; ++i) {
            Linear::Matrix<T> tmp{nRows - 1, nCols - 1};

            for (size_t j = 0; j < nRows - 1; ++j) {
                for (size_t k = 0; k < nCols; ++k) {
                    if (k < i)
                        tmp.at(j, k) = matrix.at(j + 1, k);
                    else if (k == i)
                        continue;
                    else
                        tmp.at(j, k - 1) = matrix.at(j + 1, k);
                }
            }
            ans += (i % 2 == 0 ? 1 : -1) * matrix.at(0, i) * Determinant::std_alg(tmp);
        }

        return ans;
    }
}

template <typename T, typename>
T Linear::Determinant::gauss(const Linear::Matrix<T> &matrix) {
    auto shape = matrix.shape();
    int nRows = shape.first;
    int nCols = shape.second;

    if (nRows != nCols)
        throw std::invalid_argument("Attempt to calculate the determinant of a non-square matrix.");

    if (nRows == 1) {
        return matrix.at(0, 0);
    } else {
        double ans = 1.0;
        Linear::Matrix<double> tmp = matrix;

        ans *= tmp.direct_gauss();
        if (std::fabs(ans) <= ACR) return 0;

        for (int i = 0; i < nRows; ++i) ans *= tmp.at(i, i);

        ans = (std::fabs(ans) <= ACR) ? 0 : ans;

        return ans;
    }
}

template <typename T>
void Linear::Matrix<T>::reverse_gauss() {
    size_t columnStartValue = nCols_ - 1;
    for (size_t i = std::min<size_t>(nRows_ - 1, columnStartValue); i >= 0; --i) {
        if (std::fabs(at(i, i)) >= ACR) {
            for (size_t j = i - 1; j >= 0; --j) add_rows(i, j, -(at(j, i) / at(i, i)));
        }
    }
}

template <typename T>
Linear::Matrix<T>::Matrix(size_t rows, size_t cols, T value)
    : detail::MatrixBuffer<T>(rows * cols), nRows_(rows), nCols_(cols) {
    for (size_t i = 0; i < size_; ++i) {
        new (data_ + i) T{value};
    }
}

template <typename T>
template <class InputIt1, class InputIt2>
Linear::Matrix<T>::Matrix(size_t rows, size_t cols, InputIt1 begin_, InputIt2 end_)
    : Matrix(rows, cols) {
    size_t size = 0;
    for (InputIt1 itt = begin_; itt != end_; size++, itt++);

    if (size != size_) throw std::invalid_argument("Vector size is not matched with matrix size.");

    for (size_t pos = 0; pos < std::min<size_t>(size_, size); ++pos)
        at(pos / nCols_, pos % nCols_) = *(begin_ + pos);
}

template <typename T>
Linear::Matrix<T>::Matrix(const Matrix &rhs)
    : detail::MatrixBuffer<T>(rhs), nRows_(rhs.nRows_), nCols_(rhs.nCols_) {}

template <typename T>
void Linear::Matrix<T>::swap(Matrix &rhs) noexcept {
    std::swap(nRows_, rhs.nRows_);
    std::swap(nCols_, rhs.nCols_);
    std::swap(size_, rhs.size_);

    detail::MatrixBuffer<T>::swap(rhs);
}

template <typename T>
Linear::Matrix<T> &Linear::Matrix<T>::operator=(const Matrix &rhs) {
    if (this != &rhs) {
        Matrix tmp{rhs};
        this->swap(tmp);
    }
    return *this;
}

template <typename T>
Linear::Matrix<T>::Matrix(Matrix &&rhs) noexcept {
    this->swap(rhs);
}

template <typename T>
Linear::Matrix<T> &Linear::Matrix<T>::operator=(Matrix &&rhs) noexcept {
    if (this != &rhs) this->swap(rhs);

    return *this;
}

template <typename T>
bool Linear::Matrix<T>::operator==(const Matrix &rhs) const {
    if (shape() != rhs.shape()) return false;

    return std::equal(data_, data_ + size(), rhs.data_, rhs.data_ + rhs.size());
}

template <typename T>
bool Linear::Matrix<T>::operator!=(const Matrix &rhs) const {
    return !(*this == rhs);
}

template <typename T>
Linear::Matrix<T> &Linear::Matrix<T>::operator*=(const T number) {
    for (size_t i = 0; i < nRows_; ++i)
        for (size_t j = 0; j < nCols_; ++j) at(i, j) *= number;

    return *this;
}

template <typename T>
Linear::Matrix<T> &Linear::Matrix<T>::operator*=(const Matrix &rhs) {
    if (nCols_ != rhs.nRows_) {
        throw std::invalid_argument("Attempt of multiplication of matrix with inappropriate size.");
    } else {
        Matrix<T> ans{nRows_, rhs.nCols_};

        for (size_t i = 0; i < nRows_; ++i)
            for (size_t j = 0; j < rhs.nCols_; ++j)
                for (size_t k = 0; k < nCols_; ++k) ans.at(i, j) += at(i, k) * rhs.at(k, j);

        *this = std::move(ans);
    }

    return *this;
}

template <typename T>
Linear::Matrix<T> &Linear::Matrix<T>::operator+=(const Matrix &rhs) {
    if (shape() != rhs.shape()) {
        throw std::invalid_argument("Matrix sizes do not match.");
    } else {
        for (size_t i = 0; i < nRows_; ++i)
            for (size_t j = 0; j < nCols_; ++j) at(i, j) += rhs.at(i, j);
    }

    return *this;
}

template <typename T>
Linear::Matrix<T> &Linear::Matrix<T>::operator-=(const Matrix &rhs) {
    if (shape() != rhs.shape()) {
        throw std::invalid_argument("Matrix sizes do not match.");
    } else {
        for (size_t i = 0; i < nRows_; ++i)
            for (size_t j = 0; j < nCols_; ++j) at(i, j) -= rhs.at(i, j);
    }

    return *this;
}

template <typename T>
Linear::Matrix<T>::operator std::vector<T>() {
    std::vector<T> ans{};
    for (size_t i = 0; i < nRows_; ++i)
        for (size_t j = 0; j < nCols_; ++j) ans.push_back(at(i, j));

    return ans;
}

template <typename T>
void Linear::Matrix<T>::resize(pair shape) {
    Matrix<T> ans{shape.first, shape.second};

    for (size_t i = 0; i < std::min<size_t>(shape.first, nRows_); ++i)
        for (size_t j = 0; j < std::min<size_t>(shape.second, nCols_); ++j) ans.at(i, j) = at(i, j);

    *this = ans;
}

template <typename T>
void Linear::Matrix<T>::transpose() {
    Matrix<T> tmp{nCols_, nRows_};

    for (size_t i = 0; i < nRows_; ++i)
        for (size_t j = 0; j < nCols_; ++j) tmp.at(j, i) = at(i, j);

    *this = std::move(tmp);
}

template <typename T>
void Linear::Matrix<T>::negate() {
    for (size_t i = 0; i < nRows_; ++i)
        for (size_t j = 0; j < nCols_; ++j) at(i, j) = -1 * at(i, j);
}

template <typename T>
void Linear::Matrix<T>::clear() {
    *this = std::move(Matrix<T>{});
}

template <typename T>
void Linear::Matrix<T>::diagonalize() {
    direct_gauss();
    reverse_gauss();
}

template <typename T>
int Linear::Matrix<T>::direct_gauss() {
    int gaussFactor = 1;
    for (size_t i = 0; i < std::min<size_t>(nRows_, nCols_) - 1; ++i) {
        T max_elem = at(i, i);
        size_t max_indx = i;
        for (size_t j = i + 1; j < nRows_; ++j) {
            if (std::fabs(at(j, i)) > std::fabs(max_elem)) {
                max_elem = at(j, i);
                max_indx = j;
            }
        }

        if (std::fabs(max_elem - T{}) <= ACR) {
            return 0;
        } else if (max_indx != i) {
            gaussFactor *= -1;
            swap_rows(i, max_indx);
        }

        for (size_t j = i + 1; j < nRows_; ++j) {
            double factor = -(at(j, i) / at(i, i));
            add_rows(i, j, factor);
        }
    }

    return gaussFactor;
}

template <typename T>
void Linear::Matrix<T>::swap_rows(int lhs, int rhs) {
    if (lhs == rhs) {
        return;
    } else if (lhs >= nRows_ || rhs >= nRows_) {
        throw std::invalid_argument("Incorrect row index value.");
    } else {
        for (size_t i = 0; i < nCols_; ++i) std::swap(at(lhs, i), at(rhs, i));
    }
}

template <typename T>
void Linear::Matrix<T>::add_rows(int src, int dest, T factor) {
    if (src >= nRows_ || dest >= nRows_ || src == dest) {
        throw std::invalid_argument("Incorrect src/dest index value.");
    } else {
        for (size_t i = 0; i < nCols_; ++i) at(dest, i) += at(src, i) * factor;
    }
}

template <typename T>
void Linear::Matrix<T>::swap_cols(int lhs, int rhs) {
    if (lhs >= nCols_ || rhs >= nCols_ || lhs == rhs) {
        throw std::invalid_argument("Incorrect col index value.");
    } else {
        for (size_t i = 0; i < nRows_; ++i) std::swap(at(i, lhs), at(i, rhs));
    }
}

template <typename T>
void Linear::Matrix<T>::add_cols(int src, int dest, T factor) {
    if (src >= nCols_ || dest >= nCols_ || src == dest) {
        throw std::invalid_argument("Incorrect src/dest index value.");
    } else {
        for (size_t i = 0; i < nRows_; ++i) at(i, dest) += at(i, src) * factor;
    }
}

template <typename T>
const T *Linear::Matrix<T>::data() const noexcept {
    return data_;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::zeros(int n) {
    Matrix<T> tmp{n};

    return tmp;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::eye(int n) {
    Matrix<T> tmp{n};
    for (size_t i = 0; i < n; ++i) tmp.at(i, i) = static_cast<T>(1);

    return tmp;
}

template <typename T>
Linear::pair Linear::Matrix<T>::shape() const {
    return pair{nRows_, nCols_};
}

template <typename T>
int Linear::Matrix<T>::size() const noexcept {
    return size_;
}

template <typename T>
T Linear::Matrix<T>::trace() const {
    T ans{};
    for (size_t i = 0; i < std::min<size_t>(nRows_, nCols_); ++i) ans += at(i, i);

    return ans;
}

template <typename T>
const T &Linear::Matrix<T>::at(size_t i, size_t j) const {
    if (i >= nRows_ || j >= nCols_)
        throw std::invalid_argument("Incorrect rows/cols number value.");
    else
        return data_[i * nCols_ + j];
}

template <typename T>
void Linear::Matrix<T>::dump(std::ostream &stream) const {
    stream.precision(2);
    for (size_t i = 0; i < nRows_; ++i) {
        for (size_t j = 0; j < nCols_; ++j) {
            stream << std::left << std::setw(10) << at(i, j);
            if (!stream.good()) throw std::runtime_error("stream error: invalid element");
        }
        stream << std::endl;
    }
}

template <typename T>
T &Linear::Matrix<T>::at(size_t i, size_t j) {
    if (i >= nRows_ || j >= nCols_)
        throw std::invalid_argument("Incorrect rows/cols number value.");
    else
        return data_[i * nCols_ + j];
}

template <typename T>
T Linear::Matrix<T>::determinant(Determinant::Type type) const {
    switch (type) {
        case Determinant::Type::STD:
            return Determinant::std_alg(*this);
        case Determinant::Type::GAUSS:
            return Determinant::gauss(*this);
        default:
            throw std::runtime_error("Invalid determinant type.");
            return 0;
    }
}

template <typename T>
std::istream &Linear::operator>>(std::istream &stream, Matrix<T> &rhs) {
    int rows = 0, cols = 0;
    stream >> rows;
    if (!stream.good()) throw std::runtime_error("stream error: invalid argument");

    stream >> cols;
    if (!stream.good()) throw std::runtime_error("stream error: invalid argument");

    std::vector<T> tmp(rows * cols);
    for (T &element : tmp) {
        stream >> element;
        if (!stream.good()) throw std::runtime_error("stream error: invalid element");
    }

    rhs = {rows, cols, tmp};
    return stream;
}

template <typename T>
std::ostream &Linear::operator<<(std::ostream &stream, const Matrix<T> &rhs) {
    rhs.dump(stream);
    return stream;
}

template <typename T>
Linear::Matrix<T> Linear::operator+(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    Matrix<T> tmp{lhs};
    tmp += rhs;

    return tmp;
}

template <typename T>
Linear::Matrix<T> Linear::operator-(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    Matrix<T> tmp{lhs};
    tmp -= rhs;

    return tmp;
}

template <typename T>
Linear::Matrix<T> Linear::operator*(const T number, const Matrix<T> &rhs) {
    Matrix<T> tmp{rhs};
    tmp *= number;

    return tmp;
}

template <typename T>
Linear::Matrix<T> Linear::operator*(const Matrix<T> &lhs, const T number) {
    return number * lhs;
}

template <typename T>
Linear::Matrix<T> Linear::operator*(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    Matrix<T> tmp{lhs};
    tmp *= rhs;

    return tmp;
}
