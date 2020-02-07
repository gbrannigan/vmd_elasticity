// FloatArray2D-related functions. The implementation itself is actually just in the header file.
// Here we only have an operator<< to aid in printing a FloatArray2D.
#include <iostream>
#include "floatarray2d.h"

std::ostream& operator<<(std::ostream &os, const FloatArray2D &arr) {
    os << "FloatArray2D " << arr.num_rows() << "x" << arr.num_cols() << " [" << std::endl;
    for(int r = 0; r < arr.num_rows(); r++) {
        os << "  ";
        for(int c = 0; c < arr.num_cols(); c++) {
            os << arr.get(r, c) << " ";
        }
        os << std::endl;
    }
    os << "] ";
    return os;
}
