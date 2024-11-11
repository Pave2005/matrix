#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "matrix.h"
#include "funcs.h"

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
