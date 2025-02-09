#include <matrix.h>

int main()
try
{
    if (std::cin.bad()) throw std::runtime_error("stream error");

    int size = 0;
    std::cin >> size;

    Linear::Matrix<double> matrix {size};

    for (size_t row = 0; row < size; ++row)
    {
        for (size_t col = 0; col < size; ++col)
        {
            std::cin >> matrix.at(row, col);
        }
    }

    std::cout << matrix.determinant(Linear::Determinant::Type::GAUSS) << std::endl;
}
catch (std::exception& expt)
{
    std::cerr << expt.what() << std::endl;
    return 1;
}
