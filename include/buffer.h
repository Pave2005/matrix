#pragma once

#include <algorithm>

template <typename T>
class MatrixBuffer
{
    protected:
        size_t size_ = 0;
        size_t used_ = 0;
        T*  data_ = nullptr;

        MatrixBuffer (size_t size = 0) : size_(size), used_(0), data_(nullptr)
        {
            int memorySize = size * sizeof(T);
            data_ = (size == 0) ? nullptr : static_cast<T*>(::operator new [] (memorySize));
        }

        void swap (MatrixBuffer& rhs)
        {
            std::swap (this->size_, rhs.size_);
            std::swap (this->used_, rhs.used_);
            std::swap (this->data_, rhs.data_);
        }

        ~MatrixBuffer ()
        {
            for (int i = 0; i < used_; ++i) data_[i].~T();
            operator delete [] (data_);
        }

    public:
        MatrixBuffer (const MatrixBuffer& rhs) : MatrixBuffer<T>(rhs.size_)
        {
            std::copy(rhs.data_, rhs.data_ + rhs.size_, data_);

            this->used_ = rhs.used_;
        }

        MatrixBuffer& operator= (const MatrixBuffer& rhs)
        {
            MatrixBuffer lhs {rhs};
            size_ = rhs.size_;
            used_ = rhs.used_;

            data_ = std::exchange(lhs.data_, nullptr );

            return *this;
        }


        MatrixBuffer (MatrixBuffer&& rhs) noexcept : size_(rhs.size_), used_(rhs.used_),
                                                     data_(std::exchange(rhs.data_, nullptr)) {}

        MatrixBuffer& operator= (MatrixBuffer&& rhs) noexcept
        {
            std::swap (size_, rhs.size_);
            std::swap (used_, rhs.used_);
            std::swap (data_, rhs.data_);
            return *this;
        }
};
