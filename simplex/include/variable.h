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

#ifndef CPPLEX_VARIABLE_H
#define CPPLEX_VARIABLE_H

#include "constraint.h"

namespace cpplex {
    template <typename Scalar> class Simplex;
    template <typename Scalar> class AuxiliaryVariable;

    /**
        Variable
        ========
        Base variable class.

    */
    template <typename Scalar>
    class Variable {

        friend class Simplex<Scalar>;

    public:

        Variable( Simplex<Scalar> * creator, char const * name) :
            name(name),
            creator(creator)
        { }

        virtual ~Variable() { }

        virtual void process(Matrix<Scalar>& calculated_solution, Matrix<Scalar>& solution, int index) {
            solution(index) = calculated_solution(index);
        }

    protected:
        std::string name;
        Simplex<Scalar> * creator;

    };

    /**
        SplittedVariable
        ================
        Variables that are splitted in two (one) AuxiliaryVariables during
        the translation in standard form.

    */
    template <typename Scalar>
    class SplittedVariable : public Variable<Scalar> {

        friend class Simplex<Scalar>;

    public:
        SplittedVariable( Simplex<Scalar>* creator,
                          char const * name,
                          AuxiliaryVariable<Scalar>* aux) :
            Variable<Scalar>(creator, name),
            aux(aux)
        { }

        ~SplittedVariable() { }

        void process(Matrix<Scalar>& calculated_solution, Matrix<Scalar>& solution, int index) {
            solution(index) = calculated_solution(index) - calculated_solution(aux->index);
        }

    private:
        AuxiliaryVariable<Scalar>* aux;

    };

    /**
        SlackVariable
        =============
        Type of variable added when transforming a <= or >= constraint
        into a = constraint.

    */
    template <typename Scalar>
    class SlackVariable : public Variable<Scalar> {

        friend class Simplex<Scalar>;

    public:
        SlackVariable(Simplex<Scalar>* creator, char const * name) :
            Variable<Scalar>(creator, name)
        { }

        ~SlackVariable() { }

        void process(Matrix<Scalar>& calculated_solution, Matrix<Scalar>& solution, int index) { }

    };

    /**
        AuxiliaryVariable
        =================
        Variable created when transforming a variable in a splitted
        variable. The relation:

            x = x+ - x-

        holds between the original variable, the SplittedVariable and
        the AuxiliaryVariable.

    */
    template <typename Scalar>
    class AuxiliaryVariable : public Variable<Scalar> {

        friend class Simplex<Scalar>;
        friend class SplittedVariable<Scalar>;

    public:
        AuxiliaryVariable(Simplex<Scalar>* creator, char const * name, int index) :
            Variable<Scalar>(creator, name),
            index(index)
        { }

        ~AuxiliaryVariable() { }

        void process(Matrix<Scalar>& calculated_solution, Matrix<Scalar>& solution, int index) { }

    private:
        int index;

    };
}

#endif
