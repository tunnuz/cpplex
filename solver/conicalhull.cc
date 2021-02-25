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

#include <iostream>
#include "simplex.h"

using namespace cpplex;

using Scalar = long double;

int main( int, char** ) {
    Matrix<Scalar> conical_vectors(5, 3);
    conical_vectors <<
        0.5,    -0.866025   , -1.71266e-15,
        0.5,     0.737725   , -0.453609,
        0.5,    -0.288675   ,  0.816497,
        0.5,     0.673575   ,  0.544331,
        1.0,    -2.02968e-15, -1.71266e-15;

    std::cout << "Normal vectors forming cone:\n"
              << conical_vectors << "\n" << std::endl;
    std::cout << conical_vectors.transpose() << "\n" << std::endl;

    const int d = conical_vectors.cols();
    const int nvars = conical_vectors.rows()-1;

    std::vector<bool> redundant(conical_vectors.rows());

    for (int row = 0; row < conical_vectors.rows(); ++row)
    {
        Simplex<Scalar> problem("conical hull", NORMAL);

        for (int var = 0; var < nvars; ++var)
        {
            std::string variable_name = "x" + std::to_string(var);
            problem.add_variable(new Variable<Scalar>(&problem, variable_name.c_str()));

            RowVector<Scalar> eye = RowVector<Scalar>::Zero(nvars);
            eye(var) = 1;
            problem.add_constraint( Constraint<Scalar>( eye, CT_NON_NEGATIVE, 0 ) );
        }

        // Remove the current vector from the set of constraints. If there are still solutions,
        // then this vector is not a member of the conical hull.
        Matrix<Scalar> constraints = conical_vectors.transpose();
        constraints.block(0, row, d, nvars-row) = constraints.block(0, row+1, d, nvars-row);
        constraints.conservativeResize(d, nvars);

        RowVector<Scalar> x = conical_vectors.row(row);

        for (int c = 0; c < d; ++c)
        {
            problem.add_constraint( Constraint<Scalar>( constraints.row(c), CT_EQUAL, x(c) ) );
        }

        auto costs = RowVector<Scalar>::Zero(nvars);
        problem.set_objective_function( ObjectiveFunction<Scalar>( OFT_MAXIMIZE, costs) );

        problem.solve();
        std::cout << std::endl;
        assert(not problem.must_be_fixed());

        if ( problem.has_solutions() ) {
            if ( !problem.is_unlimited() )
                problem.print_solution();
            else
                std::cout << "Problem is unlimited." << std::endl;

        } else {
            std::cout << "Problem is overconstrained." << std::endl;
        }

        std::cout << std::endl;
    }

    return 0;
}
