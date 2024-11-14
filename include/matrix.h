#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <exception>

#include <fstream>

#include "buffer.h"
#include "../mode/mode.h"

namespace Linear
{
	using int_pair = std::pair<int, int>;

	template <typename T>
	class Matrix;

	namespace Determinant
    {
		enum class Type
        {
			STD   = 0,
			GAUSS = 1,
		};

		template <typename T>
		T      std_alg (const Linear::Matrix<T>& matrix);
	    double gauss   (const Linear::Matrix<double>& matrix);
	}

	template <typename T>
	class Matrix : private MatrixBuffer<T>
    {
		private:
			int nRows_ = 0;
            int nCols_ = 0;

			using MatrixBuffer<T>::size_;
			using MatrixBuffer<T>::used_;
			using MatrixBuffer<T>::data_;

			void reverse_gauss ();

		public:
			Matrix (int rows, int cols, T value = T {});

			Matrix (int rows = 0) : Matrix (rows, rows) {}

			~Matrix () = default;

			Matrix (int rows, int cols, std::vector<T>& vec);

			Matrix (const Matrix& rhs);

			Matrix (Matrix&& rhs);

			Matrix& operator=  (const Matrix& rhs);
			Matrix& operator=  (Matrix&& rhs);
			bool 	operator== (const Matrix& rhs) const;
			bool 	operator!= (const Matrix& rhs) const;
			Matrix& operator*= (const T number);
			Matrix& operator*= (const Matrix& rhs);
			Matrix& operator+= (const Matrix& rhs);
			Matrix& operator-= (const Matrix& rhs);

			operator std::vector<T> ();

			void resize 	 (int_pair shape);
			void transpose 	 ();
			void negate 	 ();
			void clear 		 ();
			void diagonalize ();

			int  direct_gauss ();
			void make_eye 	  ();

			void swap_rows 	 (int lhs, int rhs);
			void add_rows	 (int src, int dest, T factor);
			void append_rows (Matrix <T>& additional, bool inFront = false);
			void swap_cols 	 (int lhs, int rhs);
			void add_cols 	 (int src, int dest, T factor);
			void append_cols (Matrix <T>& additional, bool inFront = false);


			static Matrix zeros (int n);
			static Matrix eye 	(int n);

			int_pair shape () const;
			int 	 size  () const;
			T 		 trace () const;
			const T& at    (int i, int j) const;
			void 	 dump  (std::ostream& stream) const;

			T& at (int i, int j);

			T   determinant (Determinant::Type type = Determinant::Type::STD) const;
			int rank        () const;
	};

	template <typename T>
	std::istream& operator>> (std::istream& stream, Matrix<T>& rhs);
	template <typename T>
	std::ostream& operator<< (std::ostream& stream, const Matrix<T>& rhs);

	template <typename T>
	Matrix<T> operator+ (const Matrix <T>& lhs, const Matrix <T>& rhs);
	template <typename T>
	Matrix<T> operator- (const Matrix <T>& lhs, const Matrix <T>& rhs);
	template <typename T>
	Matrix<T> operator* (const T number, const Matrix <T>& rhs);
	template <typename T>
	Matrix<T> operator* (const Matrix <T>& lhs, const T number);
	template <typename T>
	Matrix<T> operator* (const Matrix <T>& lhs, const Matrix <T>& rhs);
}


template <typename T>
T Linear::Determinant::std_alg (const Linear::Matrix <T>& matrix)
{
	auto shape = matrix.shape();
	int nRows  = shape.first;
    int nCols  = shape.second;

	if (nRows != nCols) throw std::invalid_argument("Attempt to calculate the determinant of a non-square matrix.");

	if (nRows == 1) return matrix.at(0, 0);
	else
    {
		T ans {};
		for (int i = 0; i < nCols; ++i)
        {
			Linear::Matrix <T> tmp {nRows - 1, nCols - 1};

			for (int j = 0; j < nRows - 1; ++j)
            {
				for (int k = 0; k < nCols; ++k)
                {
					if      (k < i)  tmp.at(j, k) = matrix.at(j + 1, k);
					else if (k == i) continue;
					else             tmp.at(j, k - 1) = matrix.at(j + 1, k);
				}
			}
			ans += (i % 2 == 0 ? 1 : - 1) * matrix.at(0, i) * Determinant::std_alg(tmp);
		}

		return ans;
	}
}

double Linear::Determinant::gauss (const Linear::Matrix<double>& matrix)
{
	auto shape = matrix.shape();
	int  nRows = shape.first;
	int  nCols = shape.second;

	if (nRows != nCols) throw std::invalid_argument("Attempt to calculate the determinant of a non-square matrix.");

	if (nRows == 1) return matrix.at(0, 0);
	else
	{
		double ans = 1.0;
		Linear::Matrix<double> tmp = matrix;

		ans *= tmp.direct_gauss();

		for (int i = 0; i < nRows; ++i) ans *= tmp.at(i, i);

		ans = (std::fabs(ans) <= ACR) ? 0 : ans;

		return ans;
	}
}

template <typename T>
void Linear::Matrix<T>::reverse_gauss ()
{
	int columnStartValue = nCols_ - 1;
    for (int i = std::min(nRows_ - 1, columnStartValue); i >= 0; --i)
        if (std::fabs(at(i, i)) >= ACR)
            for (int j = i - 1; j >= 0; --j) add_rows(i, j, - (at(j, i) / at(i, i)));
}

template <typename T>
Linear::Matrix <T>::Matrix (int rows, int cols, T value) : MatrixBuffer<T>(rows * cols),
	                                                       nRows_(rows),
	                                                       nCols_(cols)
{
	if (nRows_ < 0 || nCols_ < 0) throw std::invalid_argument("Incorrect number of rows or columns of the matrix.");

	else if (nRows_ * nCols_ == 0)
	{
		nRows_ = 0;
		nCols_ = 0;
	}
	else
	{
		for (int i = 0; i < nRows_ * nCols_; ++i)
		{
			new (data_ + i) T {value};
			++size_;
		}
	}
}

template <typename T>
Linear::Matrix <T>::Matrix (int rows, int cols, std::vector<T> &vec) : Matrix (rows, cols)
{
	if (vec.size() != nRows_ * nCols_) throw std::invalid_argument("Vector size is not matched with matrix size.");

	for (int pos = 0; pos < std::min<int>(nRows_ * nCols_, vec.size()); ++pos)
		at(pos / nCols_, pos % nCols_) = vec[pos];
}

template <typename T>
Linear::Matrix<T>::Matrix (const Matrix& rhs) : MatrixBuffer<T>(rhs.nRows_ * rhs.nCols_),
												nRows_(rhs.nRows_),
												nCols_(rhs.nCols_)
{
	for (int i = 0; i < nRows_ * nCols_; ++i)
	{
		new (data_ + i) T {rhs.data_[i]};
		++size_;
	}
}

template <typename T>
Linear::Matrix<T>::Matrix (Matrix&& rhs)
{
	std::swap(nRows_, rhs.nRows_);
	std::swap(nCols_, rhs.nCols_);
	MatrixBuffer <T>::swap(rhs);
}

template <typename T>
Linear::Matrix<T>& Linear::Matrix <T>::operator= (const Matrix& rhs)
{
	if (this != &rhs)
	{
		Matrix tmp {rhs};
		std::swap(*this, tmp);
	}
	return *this;
}

template <typename T>
Linear::Matrix<T>& Linear::Matrix <T>::operator= (Matrix&& rhs)
{
	if (this != &rhs)
	{
		std::swap(nRows_, rhs.nRows_);
		std::swap(nCols_, rhs.nCols_);
		MatrixBuffer <T>::swap(rhs);
	}
	return *this;
}

template <typename T>
bool Linear::Matrix<T>::operator== (const Matrix& rhs) const
{
	if (shape() != rhs.shape()) return false;

	for (int i = 0; i < nRows_; ++i)
		for (int j = 0; j < nCols_; ++j)
			if (at(i, j) != rhs.at(i, j)) return false;

	return true;
}

template <typename T>
bool Linear::Matrix<T>::operator!= (const Matrix& rhs) const { return !(*this == rhs); }

template <typename T>
Linear::Matrix <T>& Linear::Matrix<T>::operator*= (const T number)
{
	for (int i = 0; i < nRows_; ++i)
		for (int j = 0; j < nCols_; ++j)
			at(i, j) *= number;

	return *this;
}

template <typename T>
Linear::Matrix <T>& Linear::Matrix<T>::operator*= (const Matrix& rhs)
{
	if (this->nCols_ != rhs.nRows_) throw std::invalid_argument("Attempt of multiplication of matrix with inappropriate size.");
	else
	{
		Matrix <T> ans {this->nRows_, rhs.nCols_};

		for (int i = 0; i < this->nRows_; ++i)
			for (int j = 0; j < rhs.nCols_; ++j)
				for (int k = 0; k < this->nCols_; ++k)
					ans.at(i, j)+= at (i, k) * rhs.at (k, j);

		*this = std::move(ans);
	}

	return *this;
}

template <typename T>
Linear::Matrix<T>& Linear::Matrix<T>::operator+= (const Matrix& rhs)
{
	if (shape() != rhs.shape()) throw std::invalid_argument("Matrix sizes do not match.");
	else
	{
		for (int i = 0; i < nRows_; ++i)
			for (int j = 0; j < nCols_; ++j)
				at(i, j) += rhs.at(i, j);
	}

	return *this;
}

template <typename T>
Linear::Matrix<T>& Linear::Matrix<T>::operator-= (const Matrix& rhs)
{
	if (shape() != rhs.shape()) throw std::invalid_argument("Matrix sizes do not match.");
	else
	{
		for (int i = 0; i < nRows_; ++i)
			for (int j = 0; j < nCols_; ++j)
				at(i, j) -= rhs.at(i, j);
	}

	return *this;
}

template <typename T>
Linear::Matrix<T>::operator std::vector<T> ()
{
	std::vector<T> ans {};
	for (int i = 0; i < nRows_; ++i)
		for (int j = 0; j < nCols_; ++j)
			ans.push_back(at(i, j));

	return ans;
}

template <typename T>
void Linear::Matrix<T>::resize (int_pair shape)
{
	Matrix<T> ans {shape.first, shape.second};
	for (int i = 0; i < std::min<int>(shape.first, nRows_); ++i)
		for (int j = 0; j < std::min<int>(shape.second, nCols_); ++j)
			ans.at(i, j) = at(i, j);

	*this = ans;
}

template <typename T>
void Linear::Matrix<T>::transpose ()
{
	Matrix <T> tmp {nCols_, nRows_};
	for (int i = 0; i < nRows_; ++i)
		for (int j = 0; j < nCols_; ++j)
			tmp.at(j, i) = at(i, j);

	*this = std::move(tmp);
}

template <typename T>
void Linear::Matrix<T>::negate ()
{
	for (int i = 0; i < nRows_; ++i)
		for (int j = 0; j < nCols_; ++j)
			at(i, j) = - 1 * at(i, j);
}

template <typename T>
void Linear::Matrix<T>::clear () { *this = std::move(Matrix <T> {}); }

template <typename T>
void Linear::Matrix<T>::diagonalize ()
{
	direct_gauss();
	reverse_gauss();
}

template <typename T>
int Linear::Matrix<T>::direct_gauss ()
{
	int gaussFactor = 1;
	for (int i = 0; i < std::min<int>(nRows_, nCols_) - 1; ++i)
	{
		T max_elem = at(i, i);
		int max_indx = i;
		for (int j = i + 1; j < nRows_; ++j)
		{
			if (std::fabs(at(j, i)) > std::fabs(max_elem))
			{
				max_elem = at(j, i);
				max_indx = j;
			}
		}

		if (std::fabs(max_elem - T{}) <= ACR) return 0;
		else if (max_indx != i)
		{
			gaussFactor *= -1;
			swap_rows(i, max_indx);
		}

		for (int j = i + 1; j < nRows_; ++j)
		{
			double factor = - (at(j, i) / at(i, i));
			add_rows(i, j, factor);
		}
	}

	return gaussFactor;
}

template <typename T>
void Linear::Matrix<T>::make_eye ()
{
	diagonalize();

	int columnStartValue = nCols_ - 1;
    for (int i = std::min(nRows_ - 1, columnStartValue); i >= 0; --i)
	{
        if (std::fabs(at(i, i)) >= ACR)
		{
            T divisor = at(i, i);
            for (int j = nCols_ - 1; j >= i; --j) at(i, j) /= divisor;
        }
    }
}

template <typename T>
void Linear::Matrix<T>::swap_rows (int lhs, int rhs)
{
	if (lhs == rhs) return;
	else if (lhs >= nRows_ || rhs >= nRows_) throw std::invalid_argument("Incorrect row index value.");
	else for (int i = 0; i < nCols_; ++i) std::swap(at(lhs, i), at(rhs, i));
}

template <typename T>
void Linear::Matrix<T>::add_rows (int src, int dest, T factor)
{
	if (src >= nRows_ || dest >= nRows_ || src == dest) throw std::invalid_argument("Incorrect src/dest index value.");
	else
	{
		for (int i = 0; i < nCols_; ++i) at(dest, i) += at(src, i) * factor;
	}
}

template <typename T>
void Linear::Matrix<T>::append_rows (Matrix<T>& add_matrix, bool inFront)
{
	if (add_matrix.nCols_ != nCols_) throw std::invalid_argument("Colomns numbers do not match.");

    auto matrix_vec = static_cast<std::vector<T>>(*this);
	auto add_matrix_vec = static_cast<std::vector<T>>(add_matrix);
	Linear::Matrix <T> result {};
	if (inFront)
	{
		add_matrix_vec.insert(add_matrix_vec.end(), matrix_vec.begin(), matrix_vec.end());
		result = {nRows_ + add_matrix.nRows_, nCols_, add_matrix_vec};
	}
	else
	{
		matrix_vec.insert(matrix_vec.end(), add_matrix_vec.begin(), add_matrix_vec.end());
		result = {nRows_ + add_matrix.nRows_, nCols_, matrix_vec};
	}
	*this = result;
}

template <typename T>
void Linear::Matrix<T>::swap_cols (int lhs, int rhs)
{
	if (lhs >= nCols_ || rhs >= nCols_ || lhs == rhs) throw std::invalid_argument("Incorrect col index value.");
	else for (int i = 0;  i < nRows_; ++i) std::swap(at(i, lhs), at(i, rhs));
}

template <typename T>
void Linear::Matrix<T>::add_cols (int src, int dest, T factor)
{
	if (src >= nCols_ || dest >= nCols_ || src == dest) throw std::invalid_argument("Incorrect src/dest index value.");
	else for (int i = 0; i < nRows_; ++i) at(i, dest) += at(i, src) * factor;
}

template <typename T>
void Linear::Matrix<T>::append_cols (Matrix<T>& add_matrix, bool inFront)
{
	if (add_matrix.nRows_ != this->nRows_) throw std::invalid_argument("Rows numbers do not match.");

    transpose();
    auto matrix_vec     = static_cast<std::vector <T>>(*this);
	auto add_matrix_vec = static_cast<std::vector <T>>(add_matrix);
	transpose();

	Linear::Matrix<T> result {};

	if (inFront)
	{
		add_matrix_vec.insert(add_matrix_vec.end(), matrix_vec.begin(), matrix_vec.end());
		result = {nCols_ + add_matrix.nCols_, nRows_, add_matrix_vec};
	}
	else
	{
		matrix_vec.insert(matrix_vec.end(), add_matrix_vec.begin(), add_matrix_vec.end());
		result = {nCols_ + add_matrix.nCols_, nRows_, matrix_vec};
	}
	result.transpose();
	*this = result;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::zeros (int n)
{
	Matrix <T> tmp {n};

	return tmp;
}

template <typename T>
Linear::Matrix<T> Linear::Matrix<T>::eye (int n)
{
	Matrix <T> tmp {n};
	for (int i = 0; i < n; ++i) tmp.at(i, i) = static_cast<T>(1);

	return tmp;
}

template <typename T>
Linear::int_pair Linear::Matrix<T>::shape () const
{
	return int_pair {nRows_, nCols_};
}

template <typename T>
int Linear::Matrix <T>::size () const
{
	return nRows_ * nCols_;
}

template <typename T>
T Linear::Matrix<T>::trace () const
{
	T ans {};
	for (int i = 0; i < std::min<int>(nRows_, nCols_); ++i) ans += at(i, i);

	return ans;
}

template <typename T>
const T& Linear::Matrix<T>::at (int i, int j) const
{
	if (i >= nRows_ || j >= nCols_) throw std::invalid_argument("Incorrect rows/cols number value.");
	else 							return data_[i * nCols_ + j];
}

template <typename T>
void Linear::Matrix<T>::dump (std::ostream& stream) const
{
	stream.precision(2);
	for (int i = 0; i < nRows_; ++i)
	{
		for (int j = 0; j < nCols_; ++j) stream << std::left << std::setw(10) << at(i, j);
		stream << std::endl;
	}
}

template <typename T>
T& Linear::Matrix<T>::at (int i, int j)
{
	return const_cast<T&>(static_cast<const Matrix<T>*>(this)->at(i, j));
}

template <typename T>
T Linear::Matrix<T>::determinant (Determinant::Type type) const
{
	switch (type)
	{
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
int Linear::Matrix<T>::rank () const
{
	Matrix<T>tmp = *this;
	tmp.make_eye(false);
	int ans = 0;
	for (int i = 0; i < std::min<int>(nRows_, nCols_); ++i)
	{
		if (std::fabs(at(i, i)) > ACR)
		{
			++ans;
			continue;
		}
	}

	return ans;
}

template <typename T>
std::istream& Linear::operator>> (std::istream& stream, Matrix <T>& rhs)
{
	int rows = 0, cols = 0;
	stream >> rows >> cols;

	std::vector <T> tmp (rows * cols);
	for (T& element : tmp) stream >> element;

	rhs = {rows, cols, tmp};
	return stream;
}

template <typename T>
std::ostream& Linear::operator<< (std::ostream& stream, const Matrix<T>& rhs)
{
	rhs.dump(stream);
	return stream;
}

template <typename T>
Linear::Matrix<T> Linear::operator+ (const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix <T> tmp {lhs};
	tmp += rhs;

	return tmp;
}

template <typename T>
Linear::Matrix<T> Linear::operator- (const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix <T> tmp {lhs};
	tmp -= rhs;

	return tmp;
}

template <typename T>
Linear::Matrix<T> Linear::operator* (const T number, const Matrix<T>& rhs)
{
	Matrix <T> tmp {rhs};
	tmp *= number;

	return tmp;
}

template <typename T>
Linear::Matrix<T> Linear::operator* (const Matrix<T>& lhs, const T number)
{
	return number * lhs;
}

template <typename T>
Linear::Matrix<T> Linear::operator* (const Matrix<T>& lhs, const Matrix<T>& rhs)
{
	Matrix <T> tmp {lhs};
	tmp *= rhs;

	return tmp;
}
