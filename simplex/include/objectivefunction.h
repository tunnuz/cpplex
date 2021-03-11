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

#ifndef CPPLEX_OBJECTIVE_FUNCTION_H
#define CPPLEX_OBJECTIVE_FUNCTION_H

#include "matrix.h"

// From the STL
#include <iostream>

namespace cpplex {

    enum ObjectiveFunctionType {

        OFT_MAXIMIZE,
        OFT_MINIMIZE

    };

    template <typename Scalar> class Simplex;

    template <typename Scalar>
    class ObjectiveFunction {

        friend class Simplex<Scalar>;

    public:

        ObjectiveFunction() = default;
        ObjectiveFunction( ObjectiveFunctionType type, RowVector<Scalar> const & costs ) :
            type(type),
            costs(costs)
        { }

        ObjectiveFunction<Scalar>& operator=( ObjectiveFunction<Scalar> const & objective_function ) = default;

        // Solution value
        auto get_value( RowVector<Scalar> const & x ) const {
            return costs * x;
        }

        // Manipulation
        void add_column(Scalar value) {
            costs.conservativeResize(costs.size()+1);
            costs(costs.size()-1) = value;
        }

        // Debug
        void log() const {

            if (type == OFT_MINIMIZE)
                std::cout << "min ( ";
            else
                std::cout << "max ( ";

            for (int i = 0; i < costs.cols(); ++i)
                std::cout << costs(i) << "  ";
            std::cout << ") * x " << std::endl;

        }

    private:

        ObjectiveFunctionType type;
        RowVector<Scalar> costs;

    };
}

#endif
