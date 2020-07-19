#pragma once
// Lightweight ghetto 2D array library. Sorry, all the ones with suitable licenses I could find were
// much heavier weight than we need for this.
//
// The design philosophy of this class is that it never silently allocates memory that the caller has
// to free. So all operations work on the object instance in place. This means that there is no operator
// overloading, as this would require allocation of intermediate instances that would need to be cleaned
// up, and I don't want to deal with that mess.
//
// The in-place operations are intended to be chained, for brevity.
//
// Tom Joseph <thomas.joseph@uphs.upenn.edu>
#include <iostream>
#include <csignal>
#include <cstring>

class FloatArray2D {
    float *data_ = nullptr;
    int num_rows_, num_cols_;

    // Sanity bounds check. I guess this could be an assert, because no user data should cause this
    // type of error unless something else is flagrantly wrong
    void check_out_of_bounds(int r, int c) const {
        if(r > num_rows_ || c > num_cols_ || r < 0 || c < 0) {
            std::cerr << "FloatArray2D: out of bounds: trying to access (" << r << ", " << c << ") but bounds are ("
                << num_rows_ << ", " << num_cols_ << ")" << std::endl;
            std::raise(SIGABRT);
        }
    }

public:
    // Your standard constructor. Allocates and zeroes space as specified.
    FloatArray2D(int num_rows, int num_cols) : num_rows_(num_rows), num_cols_(num_cols) {
        reset(num_rows, num_cols);
    }

    FloatArray2D() : num_rows_(0), num_cols_(0) {
        // OK to have an empty array. We will FloatArray2D::reset() it later
    }

    ~FloatArray2D() {
        delete[] data_;
    }

    // Resets to specified dimensions, freeing any existing data, and initializing the data to zeros
    FloatArray2D& reset(int num_rows, int num_cols) {
        delete[] data_;

        if(num_rows == 0 || num_cols == 0) {
            std::cerr << "FloatArray2D::reset: Asked to reset to nothing. That's weird" << std::endl;
            exit(EXIT_FAILURE);
        }

        data_ = new float[num_rows * num_cols];
        num_rows_ = num_rows;
        num_cols_ = num_cols;
        zero();
        return *this;
    }

    // Frees any existing data, initializes data to zero, but keeps existing dimensions
    FloatArray2D& reset() {
        reset(num_rows_, num_cols_);
        return *this;
    }

    // Initialize by copying from another instance. Frees any existing data.
    // Optionally only copies a chunk of rows from the source array
    FloatArray2D& init_from(const FloatArray2D &src, int first_row = 0, int row_count = 0) {
        // Check that there's enough data in the source to copy in the first place
        if(src.num_rows() < (first_row + row_count)) {
            std::cerr << "FloatArray2D::init_from: Not enough data in source to copy into this new instance" << std::endl;
            std::raise(SIGABRT);
        }

        if(row_count == 0) {
            row_count = src.num_rows();
        }
        reset(row_count, src.num_cols());

        // Depends heavily on row-major ordering
        std::memcpy((void *) data_, (void *) (src.data_ + first_row*src.num_cols()),
                    sizeof(float)*row_count*src.num_cols());
        return *this;
    }

    // Gets an element at the specified position.
    // Does a bounds check.
    float get(int r, int c) const {
        check_out_of_bounds(r, c);
        return data_[r*num_cols_ + c];
    }

    // Sets an element to the specified value.
    void set(int r, int c, float val) {
        check_out_of_bounds(r, c);
        data_[r*num_cols_ + c] = val;
    }

    // Returns the sum of all elements.
    // This is a double because perhaps a float won't have enough bits.
    // Yes, the caller could easily do this by calling get() a lot but if we do it this way
    // it's faster.
    double sum() const {
        auto total = 0.0;
        for(int i = 0; i < (num_rows_ * num_cols_); i++) {
            total += data_[i];
        }
        return total;
    }

    // Adds in place to a single cell
    FloatArray2D& add(int r, int c, float val) {
        check_out_of_bounds(r, c);
        data_[r*num_cols_ + c] += val;
        return *this;
    }

    // Pairwise add in place with another FloatArray2D
    FloatArray2D& add(const FloatArray2D &other) {
        if(other.num_rows() != num_rows_ || other.num_cols() != num_cols_) {
            std::cerr << "FloatArray2D::add: Operand dimensions are different from mine" << std::endl;
            std::raise(SIGABRT);
        }
        for(int i = 0; i < (num_rows_*num_cols_); i++) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    // Pairwise subtract another FloatArray2D, in place. That is to say: "this -= other"
    FloatArray2D& subtract(const FloatArray2D &other) {
        if(other.num_rows() != num_rows_ || other.num_cols() != num_cols_) {
            std::cerr << "FloatArray2D::subtract: Operand dimensions are different from mine" << std::endl;
            std::raise(SIGABRT);
        }
        for(int i = 0; i < (num_rows_*num_cols_); i++) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    // Multiply everything in place by specified value
    FloatArray2D& multiply(float val) {
        for(int i = 0; i < (num_rows_*num_cols_); i++) {
            data_[i] *= val;
        }
        return *this;
    }

    // Pairwise multiply with another FloatArray2D, in place. This is NOT a matrix multiplication.
    FloatArray2D& multiply(const FloatArray2D &other) {
        if(other.num_rows() != num_rows_ || other.num_cols() != num_cols_) {
            std::cerr << "FloatArray2D::multiply: Operand dimensions are different from mine" << std::endl;
            std::raise(SIGABRT);
        }
        for(int i = 0; i < (num_rows_*num_cols_); i++) {
            data_[i] *= other.data_[i];
        }
        return *this;
    }

    // Divides a single cell, in place, by a value
    FloatArray2D& divide(int r, int c, float denom) {
        check_out_of_bounds(r, c);
        data_[r*num_cols_ + c] /= denom;
        return *this;
    }

    // Pairwise divide with another FloatArray2D in place.
    FloatArray2D& divide(const FloatArray2D &denom) {
        if(denom.num_rows() != num_rows_ || denom.num_cols() != num_cols_) {
            std::cerr << "FloatArray2D::divide: Operand denom dimensions are different from mine" << std::endl;
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < (num_rows_*num_cols_); i++) {
            data_[i] /= denom.data_[i];
        }
        return *this;
    }

    // Pairwise divide unless denominator is zero
    FloatArray2D& divide_ignore_zero_denom(const FloatArray2D &denom) {
        if(denom.num_rows() != num_rows_ || denom.num_cols() != num_cols_) {
            std::cerr << "FloatArray2D::divide_ignore_zero_denom: Operand denom dimensions are different from mine" << std::endl;
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < (num_rows_ * num_cols_); i++) {
            if(denom.data_[i] != 0) {
                data_[i] /= denom.data_[i];
            }
        }
        return *this;
    }

    // Pairwise divide but assume answer is zero if denominator is zero
    FloatArray2D& divide_nan_is_zero(const FloatArray2D &denom) {
        if(denom.num_rows() != num_rows_ || denom.num_cols() != num_cols_) {
            std::cerr << "FloatArray2D::divide_ignore_zero_denom: Operand denom dimensions are different from mine" << std::endl;
            exit(EXIT_FAILURE);
        }
        for(int i = 0; i < (num_rows_ * num_cols_); i++) {
            if(denom.data_[i] != 0) {
                data_[i] /= denom.data_[i];
            } else {
                data_[i] = 0.0f;
            }
        }
        return *this;
    }

    // Copies a single row from a source array into this one, in place.
    // This is meant to copy things like individual xyz tuples that comprise a row.
    FloatArray2D& copy_row(const FloatArray2D &src, int src_row, int dest_row) {
        if(num_cols_ != src.num_cols()) {
            std::cerr << "FloatArray2D::copy_row: Uh oh, attempt to copy from source with different row size" << std::endl;
            exit(EXIT_FAILURE);
        }

        float *dest_data = data_ + (dest_row*num_cols_);
        float *src_data = src.data_ + (src_row*num_cols_);
        std::memcpy((void *) dest_data, (void *) src_data, sizeof(float) * num_cols_);
        return *this;
    }

    // Copy all data into a raw output array, in row-major order
    void copy_raw_to(float *out) {
        std::memcpy((void *) out, (void *) data_, sizeof(float)*num_rows_*num_cols_);
    }

    // Zeroes us out
    FloatArray2D& zero() {
        memset(data_, 0, num_rows_*num_cols_*sizeof(float));
        return *this;
    }

    // Ostensibly these getters will be optimized into simple memory accesses
    int num_rows() const {
        return num_rows_;
    }

    int num_cols() const {
        return num_cols_;
    }
};

// Stream insertion operator overload
std::ostream& operator<<(std::ostream &os, const FloatArray2D &arr);