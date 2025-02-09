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

TEST(test6, unit_test6)
{
	std::vector<double> vec = {1, 2, 3, 4, 6, 7, 8, 12, 1};
	Linear::Matrix<double> matrix = {3, 3, vec};
	Linear::Matrix<double> lhs = matrix;

	EXPECT_TRUE(matrix.size() == lhs.size());
	size_t size = matrix.size();
	for (size_t i = 0; i < size; ++i)
	{
		EXPECT_EQ(matrix.data()[i], lhs.data()[i]);
	}
}
