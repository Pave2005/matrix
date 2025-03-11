#include "matrix.h"

#include <gtest/gtest.h>

TEST(matrix, copy_ctor)
{
	std::vector<double> vec = {1, 2, 3, 4, 6, 7, 8, 12, 1};

	Linear::Matrix<double> matrix {3, 3, vec.cbegin(), vec.cend()};
	Linear::Matrix<double> lhs = matrix;

	EXPECT_TRUE(matrix.size() == lhs.size());
	size_t size = matrix.size();
	for (size_t i = 0; i < size; ++i)
	{
		EXPECT_EQ(matrix.data()[i], lhs.data()[i]);
	}
}

TEST(matrix, move_ctor)
{
    std::vector<double> vec = {-1, 2, 3, -4, 6, 7, 8, -12, 1};

	Linear::Matrix<double> matrix {3, 3, vec.cbegin(), vec.cend()};

    std::vector<std::vector<double>> data = {};
    for (size_t i = 0; i < 3; ++i) {
        std::vector<double> row = {};
        for (auto j = 0; j < 3; ++j) row.push_back(matrix.at(i, j));
        data.push_back(row);
    }
    Linear::Matrix<double> lhs = std::move(matrix);

    for (size_t i = 0; i < 3; ++i)
        for (auto j = 0; j < 3; ++j) EXPECT_EQ(data[i][j], lhs.at(i, j));
}

TEST(matrix, equal_true)
{
	std::vector<double> vec = {1, 2, 3, 4, 6, 7, 8, 12, 1};

	Linear::Matrix<double> matrix_1 {3, 3, vec.cbegin(), vec.cend()};
	Linear::Matrix<double> matrix_2 {3, 3, vec.cbegin(), vec.cend()};

	EXPECT_TRUE(matrix_1 == matrix_2);
}

TEST(matrix, equal_false)
{
	std::vector<double> vec_1 = {1, 2, 3, 4, 6, 7, 8, 12, 1};
	std::vector<double> vec_2 = {1, 5, 3, 4, 22, 7, 81, 12, 1};

	Linear::Matrix<double> matrix_1 {3, 3, vec_1.cbegin(), vec_1.cend()};
	Linear::Matrix<double> matrix_2 {3, 3, vec_2.cbegin(), vec_2.cend()};

	EXPECT_TRUE(matrix_1 != matrix_2);
}
