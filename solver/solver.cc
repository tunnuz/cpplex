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
#include "pilal.h"
#include "simplex.h"

using namespace pilal;
using namespace optimization;

int main( int argc, char* argv[]) {

    if ( argc == 2 ) {
        Simplex problem("Simplex Instance");

        try {                                         
            problem.load_problem( argv[1] );
            
			
            // Solve
            problem.solve();         
            std::cout << std::endl;
            
            if (problem.must_be_fixed()) {
                std::cout << "Problem formulation is incorrect." << std::endl;
                return 1;
            }
            
            if ( problem.has_solutions() ) {
                if ( !problem.is_unlimited() ) 
                    problem.print_solution();
                else
                    std::cout << "Problem is unlimited." << std::endl;
                
            } else {
                std::cout << "Problem is overconstrained." << std::endl;
            }                                                           
            
        } catch ( DataMismatchException c ) {
            std::cout << "Error: " << c.error << std::endl;
        } 
        
        return 0;
    } else {  
        std::cout << "Error: omitted problem file." << std::endl;
        return 1;
    }    
    
	std::cout << "Quitting ..." << std::endl;
    
}


