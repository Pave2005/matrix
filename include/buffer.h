#pragma once

#include <algorithm>

template <typename T>
class MatrixBuffer
{
    protected:
        int capacity_ = 0;
        int size_ = 0;
        T*  data_ = nullptr;

        MatrixBuffer (int capacity = 0) : capacity_(capacity), size_(0), data_(nullptr)
        {
            int memorySize = capacity * sizeof (T);
            data_ = capacity == 0 ? nullptr : static_cast<T*>(::operator new [] (memorySize));
        }

        void swap (MatrixBuffer& rhs)
        {
            std::swap (this->capacity_, rhs.capacity_);
            std::swap (this->size_, rhs.size_);
            std::swap (this->data_, rhs.data_);
        }

        ~MatrixBuffer ()
        {
            for (int i = 0; i < size_; ++i) data_[i].~T ();
            operator delete [] (data_);
        }

    public:
        MatrixBuffer (const MatrixBuffer& rhs) = delete;

        MatrixBuffer& operator= (const MatrixBuffer& rhs) = delete;
};
