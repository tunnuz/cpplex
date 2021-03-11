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

#ifndef CPPLEX_COLUMNSET_H
#define CPPLEX_COLUMNSET_H

// From the STL
#include <iostream>
#include <vector>
#include <algorithm>

namespace cpplex {
    template <typename Scalar> class Simplex;

    /**
        ColumnSet
        =========
        Class that represents a set of columns.

    */
    template <typename Scalar>
    class ColumnSet {

        friend class Simplex<Scalar>;

    public:

        void insert ( int column ) {
            columns.push_back(column);
        }

        void remove ( int column ) {
            if ( std::find( columns.begin(), columns.end(), column) != columns.end() )
                columns.erase(std::find( columns.begin(), columns.end(), column));
        }

        void substitute ( int old_column, int new_column ) {
            if ( std::find( columns.begin(), columns.end(), old_column) != columns.end() )
                *(std::find( columns.begin(), columns.end(), old_column)) = new_column;
        }

        void log(char const * prelude) const  {

            std::cout << prelude;

            for (   std::vector<int>::const_iterator it = columns.begin();
                    it != columns.end();
                    ++it ) {
                std::cout << *it << " ";
            }

            std::cout << std::endl;
        }

        bool contains(int column) const {
            return std::find( columns.begin(), columns.end(), column) != columns.end();
        }

        int& column(int idx) {
            return columns.at(idx);
        }

        int index_of(int column) {
            int pos = 0;
            for ( std::vector<int>::iterator it = columns.begin(); it != columns.end(); ++it)
                if ( *it != column )
                    ++pos;
                else
                    return pos;
            return -1;
        }

        unsigned int size() const {
            return columns.size();
        }


    private:

        // A int vector stores the indices of the columns in the set.
        std::vector<int> columns;
    };

}


#endif
