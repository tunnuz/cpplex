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
#include "conicalhull.h"

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

    ConicalHull<Scalar> hull(conical_vectors);
    std::cout << "The following vectors are dependent on the others:\n"
              << hull.is_dependent().transpose() << "\n" << std::endl;

    std::cout << "The conical hull spans these vectors:\n"
              << hull.hull() << std::endl;

    return 0;
}
