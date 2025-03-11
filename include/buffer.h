#pragma once

#include <algorithm>
#include <memory>

namespace detail
{
template <typename T> class MatrixBuffer
{
  protected:
    size_t size_ = 0;
    T *data_ = nullptr;

    MatrixBuffer(size_t size = 0) : size_(size)
    {
        size_t memorySize = size * sizeof(T);
        data_ = (size == 0) ? nullptr : static_cast<T *>(::operator new[](memorySize));
    }

    void swap(MatrixBuffer &rhs) noexcept
    {
        std::swap(size_, rhs.size_);
        std::swap(data_, rhs.data_);
    }

    ~MatrixBuffer()
    {
        for (size_t i = 0; i < size_; ++i)
            std::destroy_at(data_ + i);
        operator delete[](data_);
    }

  public:
    MatrixBuffer(const MatrixBuffer &rhs) : MatrixBuffer<T>(rhs.size_)
    {
        std::copy(rhs.data_, rhs.data_ + rhs.size_, data_);
    }

    MatrixBuffer &operator=(const MatrixBuffer &rhs)
    {
        MatrixBuffer tmp{rhs};
        this->swap(tmp);

        return *this;
    }

    MatrixBuffer(MatrixBuffer &&rhs) noexcept
    {
        swap(rhs);
    }

    MatrixBuffer &operator=(MatrixBuffer &&rhs) noexcept
    {
        this->swap(rhs);
        return *this;
    }
};
} // namespace detail
