/*
This file is part of C++lex, a project by Tommaso Urli.

C++lex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

C++lex is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with C++lex.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CPPLEX_MATRIX_H
#define CPPLEX_MATRIX_H

#include <cstdio>
#include <Eigen/Dense>

namespace optimization {

    using Matrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using RowVector = Eigen::Matrix<long double, 1, Eigen::Dynamic>;
    using ColumnVector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;

    /** Auxiliary function, number comparison with tolerance. */
    template <typename Scalar>
    bool tol_equal(Scalar n, Scalar m, Scalar tol=1e-16) {
        return std::fabs( n - m ) < tol;
    }

    /** Compare two values (with tolerance). */
    template <typename Matrix, typename Scalar>
    bool more_equal_than (const Matrix& matrix, Scalar value, Scalar tol = 1e-16) {
        for (int i = 0; i < matrix.rows(); ++i)
            for (int j = 0; j < matrix.cols(); ++j)
                if ( matrix(i,j) + tol < value ) return false;
        return true;    
    }
	
    /** Compare two values (with tolerance). */
    template <typename Matrix, typename Scalar>
    bool less_equal_than (const Matrix& matrix, Scalar value, Scalar tol = 1e-16) {
        for (int i = 0; i < matrix.rows(); ++i)
            for (int j = 0; j < matrix.cols(); ++j)
                if ( matrix(i,j) - tol > value ) return false;
        return true;    
    }

    /** Prints the matrix with a name for debug. */
    template <typename Scalar,
              int RowsAtCompileTime,
              int ColsAtCompileTime,
              int Options,
              int MaxRowsAtCompileTime,
              int MaxColsAtCompileTime>
    void log_matrix(const Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime,
                    Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> matrix,
                    std::string name) {
        // Printing
        printf("-- %s\n", name.c_str());
        for (int i = 0; i < matrix.rows(); ++i) {
            printf(" ");
            for (int j = 0; j < matrix.cols(); ++j) {
                printf("%10.5f ", (double)matrix(i,j));
            }
            printf("\n");
        }
        printf("--\n");     
    }

    /** Updates inverse of the matrix after a column has changed. */
    template <typename Scalar,
              int RowsAtCompileTime,
              int ColsAtCompileTime,
              int Options,
              int MaxRowsAtCompileTime,
              int MaxColsAtCompileTime>
    auto update_inverse_after_column_change(Eigen::Matrix<Scalar,
                                            RowsAtCompileTime, ColsAtCompileTime,
                                            Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>
                                            const& old_inverse,
                                            ColumnVector const& new_column, int q) {

        // Prepare result
        using Matrix = Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime,
                                     Options, MaxRowsAtCompileTime, MaxColsAtCompileTime>;
        Matrix new_inverse(old_inverse.rows(), old_inverse.cols());
        ColumnVector a_tilde = old_inverse * new_column;
                	    
        for (int i = 0; i < old_inverse.rows(); ++i)
            for (int j = 0; j < old_inverse.cols(); ++j)
                if ( i != q )
                    new_inverse(i,j) = old_inverse(i,j) - (( old_inverse(q,j) * a_tilde(i) ) / a_tilde(q));
                else
                    new_inverse(i,j) = old_inverse(q,j) / a_tilde(q);
                
        return new_inverse;
    }
}

#endif
