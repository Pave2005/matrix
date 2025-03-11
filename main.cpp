#include <matrix.h>

int main()
try
{
    size_t size = 0;
    std::cin >> size;
    if (!std::cin.good())
        throw std::runtime_error("stream error: invalid size");

    Linear::Matrix<double> matrix{size};

    for (size_t row = 0; row < size; ++row)
    {
        for (size_t col = 0; col < size; ++col)
        {
            std::cin >> matrix.at(row, col);
            if (!std::cin.good())
                throw std::runtime_error("stream error: invalid element");
        }
    }

    std::cout << matrix.determinant(Linear::Determinant::Type::GAUSS) << std::endl;
}
catch (std::exception &expt)
{
    std::cerr << expt.what() << std::endl;
    return 1;
}
