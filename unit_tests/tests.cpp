#include "matrix.h"
#include "funcs.h"

#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <vector>

TEST(test1, size100)
{
	test_funcs::run_test("/test1/size100");
}

TEST(test2, size150)
{
	test_funcs::run_test("/test2/size150");
}

TEST(test3, size150)
{
	test_funcs::run_test("/test3/size150");
}

TEST(test4, size4)
{
	test_funcs::run_test("/test4/size4");
}

TEST(test5, size100)
{
	test_funcs::run_test("/test5/size100");
}

TEST(test6, copyConstrTest)
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

TEST(test7, equalOpTest)
{
	std::vector<double> vec = {1, 2, 3, 4, 6, 7, 8, 12, 1};

	Linear::Matrix<double> matrix_1 {3, 3, vec.cbegin(), vec.cend()};
	Linear::Matrix<double> matrix_2 {3, 3, vec.cbegin(), vec.cend()};

	EXPECT_TRUE(matrix_1 == matrix_2);
}

TEST(test8, unequalOpTest)
{
	std::vector<double> vec_1 = {1, 2, 3, 4, 6, 7, 8, 12, 1};
	std::vector<double> vec_2 = {1, 5, 3, 4, 22, 7, 81, 12, 1};

	Linear::Matrix<double> matrix_1 {3, 3, vec_1.cbegin(), vec_1.cend()};
	Linear::Matrix<double> matrix_2 {3, 3, vec_2.cbegin(), vec_2.cend()};

	EXPECT_TRUE(matrix_1 != matrix_2);
}
