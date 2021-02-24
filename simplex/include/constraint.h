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

#ifndef CPPLEX_CONSTRAINT_H
#define CPPLEX_CONSTRAINT_H

#include "matrix.h"
#include "datamismatchexception.h"

// From the STL
#include <iostream>

namespace optimization {

    enum ConstraintType {

        CT_LESS_EQUAL,
        CT_MORE_EQUAL,
        CT_EQUAL,
        CT_NON_NEGATIVE,
        CT_BOUNDS

    };

    template <typename Scalar> class Simplex;

    /**
        Constraint
        ==========
        Class that represents a constraint: A_i * x_j = b_i.

    */
    template <typename Scalar>
    class Constraint {

        friend class Simplex<Scalar>;

    public:

        Constraint( RowVector<Scalar> const & coefficients, ConstraintType type, Scalar value ) {

            this->coefficients = coefficients;
            this->type = type;
            this->value = value;
        }

        Constraint( RowVector<Scalar> const & coefficients, ConstraintType type, Scalar lower, Scalar upper ) {
            if ( type != CT_BOUNDS )
                throw(DataMismatchException("Invalid constraint type for provided data"));

            this->coefficients = coefficients;
            this->type = type;
            this->lower = lower;
            this->upper = upper;
        }

        // Debug
        void log() const {
            for (int i = 0; i < coefficients.size(); ++i)
                std::cout << coefficients(i) << "\t";

            switch(type) {

            case CT_EQUAL:
                std::cout << "=\t";
                break;

            case CT_LESS_EQUAL:
                std::cout << "<=\t";
                break;

            case CT_MORE_EQUAL:
                std::cout << ">=\t";
                break;

            case CT_BOUNDS:
                std::cout << "bounded to ";
                break;

            case CT_NON_NEGATIVE:
                std::cout << "non-negative ";
            }

            if (type == CT_NON_NEGATIVE)
                std::cout << std::endl;
            else if ( type == CT_BOUNDS )
                std::cout << lower << " <= " << "value" << " <= " << upper << std::endl;
            else
                std::cout << value << std::endl;
        }

        void add_column(Scalar value) {
            coefficients.conservativeResize(1, coefficients.size()+1);
            coefficients(coefficients.size()-1) = value;
        }

        int size() const {
            return coefficients.size();
        }

    private:

        ConstraintType type;
        RowVector<Scalar> coefficients;
        Scalar value;
        Scalar upper;
        Scalar lower;

    };

}

#endif
