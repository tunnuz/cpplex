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

#ifndef CPPLEX_CONICALHULL_H
#define CPPLEX_CONICALHULL_H

#include "simplex.h"

namespace cpplex
{
    template <typename Scalar>
    class ConicalHull
    {
        using LogicalVector = Eigen::Matrix<bool, Eigen::Dynamic, 1>;

    public:
        ConicalHull(Matrix<Scalar> normal_vectors,
                    Verbosity verbosity=SILENT)
            : normal_vectors(normal_vectors),
              verbosity(verbosity)
        {
            this->solve_linear_programs();
        }


        LogicalVector is_dependent() const
        {
            return this->dependents;
        }

        LogicalVector is_independent() const
        {
            return this->dependents.array() != true;
        }

        Matrix<Scalar> hull() const
        {
            LogicalVector independent = this->is_independent();
            Matrix<Scalar> hull(independent.count(), this->normal_vectors.cols());

            int i = 0;
            for (int row = 0; row < independent.size(); ++row)
            {
                if (independent(row))
                    hull.row(i++) = this->normal_vectors.row(row);
            }

            return hull;
        }

    private:
        Matrix<Scalar> normal_vectors;
        Verbosity verbosity;

        LogicalVector dependents;

        void solve_linear_programs()
        {
            this->dependents = LogicalVector::Zero(this->normal_vectors.rows());

            const int d = this->normal_vectors.cols();

            for (int row = 0; row < this->normal_vectors.rows(); ++row)
            {
                Simplex<Scalar> problem((std::string("row ") + std::to_string(1+row)).c_str(),
                                        this->verbosity);

                this->dependents(row) = true;
                const int nvars = this->normal_vectors.rows() - this->dependents.count();

                for (int var = 0; var < nvars; ++var)
                {
                    std::string variable_name = "x" + std::to_string(var);
                    problem.add_variable(new Variable<Scalar>(&problem, variable_name.c_str()));

                    RowVector<Scalar> eye = RowVector<Scalar>::Zero(nvars);
                    eye(var) = 1;
                    problem.add_constraint( Constraint<Scalar>( eye, CT_NON_NEGATIVE, 0 ) );
                }

                // Remove redundant vectors (including the current vector) from the set of
                // constraints. If there are still solutions, then this vector is not a member of
                // the conical hull.
                Matrix<Scalar> constraints(d, nvars);
                int var = 0;
                for (int row = 0; row < this->normal_vectors.rows(); ++row)
                {
                    if (not this->dependents(row))
                        constraints.col(var++) = this->normal_vectors.row(row);
                }

                RowVector<Scalar> x = this->normal_vectors.row(row);

                for (int c = 0; c < d; ++c)
                {
                    problem.add_constraint( Constraint<Scalar>( constraints.row(c), CT_EQUAL, x(c) ) );
                }

                auto costs = RowVector<Scalar>::Zero(nvars);
                problem.set_objective_function( ObjectiveFunction<Scalar>( OFT_MAXIMIZE, costs) );

                problem.solve();
                assert(not problem.must_be_fixed());

                this->dependents(row) = problem.has_solutions();
            }
        }
    };
}

#endif
