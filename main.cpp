#include <matrix.h>

int main()
{
    size_t dim = 0;
    std::cin >> dim;

    Linear::Matrix<double> matrix {dim};

    for (size_t row = 0; row < dim; ++row)
        for (size_t col = 0; col < dim; ++col)
            std::cin >> matrix.at(row, col);

    std::cout << matrix.determinant(Linear::Determinant::Type::FULL) << std::endl;

    matrix.dump(std::cout);
}
