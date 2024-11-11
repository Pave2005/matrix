#pragma once

#include <algorithm>

template <typename T>
class MatrixBuffer
{
    protected:
        int size_ = 0;
        int used_ = 0;
        T*  data_ = nullptr;

        MatrixBuffer (int size = 0) : size_(size), used_(0), data_(nullptr)
        {
            int memorySize = size * sizeof (T);
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
};
