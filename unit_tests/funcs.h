#pragma once

#include "matrix.h"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <fstream>

namespace test_funcs
{
	double get_result (const std::string& filename)
    {
        std::ifstream file(filename);
        if (file.bad()) throw std::runtime_error("file stream error");

        int size = 0;
        file >> size;

        Linear::Matrix<double> matrix {size};

        for (size_t row = 0; row < size; ++row)
        {
            for (size_t col = 0; col < size; ++col)
            {
                file >> matrix.at(row, col);
            }
        }

        file.close();

        return matrix.determinant(Linear::Determinant::Type::GAUSS);
    }

    double get_answer(const std::string& filename)
    {
        std::ifstream answer_file(filename);

        double ans = 0.0;

        if (answer_file.good())
            answer_file >> ans;
        else
            throw std::runtime_error("file stream error");

        answer_file.close();

        return ans;
    }

	int run_test (const std::string& test_name)
	{
        try
        {
            std::string test_directory = "/tests";

            std::string test_path = std::string(TEST_DATA_DIR) + test_directory + test_name;

            double res = get_result(test_path + ".dat");
            double ans = get_answer(test_path + ".ans");

            if (std::fabs(res - ans) <= ACR) res = ans;

            EXPECT_EQ(res, ans);

            return 0;
        }
        catch (std::exception& expt)
        {
            std::cerr << expt.what() << std::endl;
            return 1;
        }
	}
}
