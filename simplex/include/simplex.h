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

#ifndef CPPLEX_SIMPLEX_H
#define CPPLEX_SIMPLEX_H

// Simplex classes
#include "simplexexceptions.h"
#include "constraint.h"
#include "objectivefunction.h"
#include "columnset.h"


// From the STL
#include <iostream>
#include <vector>

namespace optimization {
    class ColumnSet;
    class Constraint;
    class ObjectiveFunction;
    class Variable;

    class Simplex {
        
        public:  
        
            // Constructor
            Simplex(char const * name);
            ~Simplex();
            
            // Settings                                
            void load_problem(char const * problem_name);
            void add_variable(Variable* variable);
            void add_constraint(Constraint const & constraint);
            void set_objective_function(ObjectiveFunction const & objective_function);
            
            // Solving procedures
            void solve();          

            // Print
            void print_solution() const;
            void log() const; 
            
            bool is_unlimited() const;
            bool has_solutions() const;
            bool must_be_fixed() const;
            RowVector const& get_dual_variables() const;
                                  
            
        protected:
                                                             
            std::string name;
            
            // Preprocessing
            void process_to_standard_form();
            void process_to_artificial_problem();
         
            // Solving
            void solve_with_base( ColumnSet const& base );
            
            // Column sets
            ColumnSet suggested_base;
            ColumnSet current_base;
            ColumnSet current_out_of_base;
            
            // Data
            ObjectiveFunction objective_function;
            std::vector< Constraint> constraints;
            std::vector< Constraint> nn_constraints;            
            std::vector< Variable*> variables;
            
            // Processed data
            RowVector costs;
            Matrix coefficients_matrix;
            ColumnVector constraints_vector;
            Matrix base_inverse;   
            RowVector dual_variables;
            Matrix column_p;
            int solution_dimension, old_column;
            
            // Results
            ColumnVector base_solution;
            Matrix solution;
            Matrix reduced_cost;
            long double solution_value;
            
            bool    optimal, 
                    unlimited, 
                    overconstrained, 
                    has_to_be_fixed, 
                    changed_sign;
            
            int inverse_recalculation_rate;
          
    };      

}

#endif
