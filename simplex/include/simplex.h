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
#include "variable.h"


// From the STL
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <vector>
#include <sys/stat.h>


namespace cpplex {
    // Auxiliary function
    inline int file_exists( const char *filename ) {
        struct stat buffer ;
        if ( stat( filename, &buffer ) ) return 1 ;
        return 0 ;
    }

    // Enum to individuate parsing context
    enum ParsingContext {
         PB_METADATA,
         PB_VARS,
         PB_CONSTRAINTS,
         PB_OBJECTIVE
    };


    /**
        Simplex
        =======
        Implementation of the class that will incapsulate the
        simplex behavior.
    */
    template <typename Scalar>
    class Simplex {
        static constexpr Scalar TOL = std::numeric_limits<Scalar>::epsilon();
        static constexpr bool VERBOSE = false;

    public:

        // Constructor
        Simplex( char const * name ) :
            name(name),
            solution_dimension(0),
            changed_sign(false),
            inverse_recalculation_rate(10) {
        }

        ~Simplex() {
            // Cleanup variables
            for (auto it = variables.begin(); it != variables.end(); ++it)
                if ( (*it)->creator == this )
                    delete *it;
        }

        // Settings

        void load_problem( char const * problem_name ) {

            if ( file_exists ( problem_name )  )
                throw ( DataMismatchException( "file not found." ) );

            /*
              File parsing
            */

            std::ifstream file(problem_name);

            ParsingContext current_parsing_block;
            int current_var = 0, solution_dimension = 0;

            if (file.is_open()) {

                while (!file.eof()) {
                    std::string buffer_init, token;
                    std::getline(file, buffer_init);
                    std::stringstream buffer(buffer_init);

                    // Extract token
                    buffer >> token;

                    if ( token.length() ) {
                        if ( token == "[METADATA]" )
                            current_parsing_block = PB_METADATA;
                        else if ( token == "[VARIABLES]" )
                            current_parsing_block = PB_VARS;
                        else if ( token == "[CONSTRAINTS]")
                            current_parsing_block = PB_CONSTRAINTS;
                        else if ( token == "[OBJECTIVE]")
                            current_parsing_block = PB_OBJECTIVE;
                        else {
                            switch ( current_parsing_block ) {

                            case PB_METADATA:
                            {
                                if ( token == "name" ) {

                                    buffer_init.erase(0,5);
                                    name = buffer_init;

                                } else if ( token == "vars" ) {
                                    buffer >> solution_dimension;
                                }
                            }
                            break;

                            case PB_VARS:
                            {

                                Matrix<Scalar> eye = Matrix<Scalar>::Zero(1,solution_dimension);
                                eye(current_var) = 1;
                                std::string variable_name, lower_bound, upper_bound;

                                lower_bound = token;
                                buffer >> variable_name;
                                buffer >> upper_bound;

                                // Add variable for tracking
                                add_variable(new Variable<Scalar>(this, variable_name.c_str()));


                                if ( lower_bound != "inf" ) {
                                    if ( atof(lower_bound.c_str()) == 0 )
                                        add_constraint( Constraint<Scalar>( eye, CT_NON_NEGATIVE, 0 ) );
                                    else
                                        add_constraint( Constraint<Scalar>( eye, CT_MORE_EQUAL, atof(lower_bound.c_str()) ) );
                                }
                                if ( upper_bound != "inf")
                                    add_constraint( Constraint<Scalar>( eye, CT_LESS_EQUAL, atof(upper_bound.c_str()) ) );

                                current_var++;
                            }
                            break;

                            case PB_CONSTRAINTS:
                            {

                                if ( current_var != solution_dimension )
                                    throw (DataMismatchException("Mismatch between declared and defined variables."));

                                Matrix<Scalar> coefficients(1, solution_dimension);
                                coefficients(0) = atof(token.c_str());

                                for ( int i = 1; i < solution_dimension; ++i)
                                    buffer >> coefficients(i);

                                std::string ct;
                                Scalar bound;
                                buffer >> ct;
                                buffer >> bound;

                                if ( ct == ">" )
                                    add_constraint( Constraint<Scalar>( coefficients, CT_MORE_EQUAL, bound ) );
                                else if ( ct== "<" )
                                    add_constraint( Constraint<Scalar>( coefficients, CT_LESS_EQUAL, bound ) );
                                else if ( ct== "=" )
                                    add_constraint( Constraint<Scalar>( coefficients, CT_EQUAL, bound ) );
                                else {
                                    std::stringstream parse_error;
                                    parse_error << " near ";
                                    parse_error << buffer_init.c_str();
                                    throw( DataMismatchException( parse_error.str().c_str() ) );
                                }

                            }
                            break;

                            case PB_OBJECTIVE:
                            {
                                std::string oft = token;
                                Matrix<Scalar> costs(1, solution_dimension);
                                for (int i = 0; i < solution_dimension; ++i)
                                    buffer >> costs(i);

                                if ( oft == "maximize" )
                                    set_objective_function( ObjectiveFunction<Scalar>( OFT_MAXIMIZE, costs) );
                                else if ( oft == "minimize" )
                                    set_objective_function( ObjectiveFunction<Scalar>( OFT_MINIMIZE, costs) );
                                else {
                                    throw( DataMismatchException( "Unknown objective function kind." ) );
                                }
                            }

                            break;
                            }
                        }
                    }
                }
            }


            file.close();
            return;
        }

        void add_variable(Variable<Scalar>* variable) {
            variables.push_back(variable);
        }

        void add_constraint ( Constraint<Scalar> const & constraint ) {

            if ( constraints.size() != 0 ) {
                if ( solution_dimension != constraint.coefficients.size() )
                    throw(DataMismatchException("Constraints must have the same size"));
            } else {
                solution_dimension = constraint.size();
            }

            if ( constraint.type == CT_NON_NEGATIVE ) {
                nn_constraints.push_back(constraint);
            } else {
                constraints.push_back(constraint);
            }
        }

        void set_objective_function ( ObjectiveFunction<Scalar> const& objective_function ) {

            if ( solution_dimension != objective_function.costs.size() )
                throw(DataMismatchException("Objective function must have same size as solution"));

            this->objective_function = objective_function;
        }

        // Solving procedures
        void solve() {

            ColumnSet<Scalar> initial_base;

            // Create alias to *this
            Simplex& original_problem = *this;

            // Create problem to work on
            Simplex standard_form_problem = original_problem;

            has_to_be_fixed = false;

            log();

            // Preprocessing
            std::cout << "Generating problem in standard form ...";
            standard_form_problem.process_to_standard_form();
            std::cout << " done." << std::endl;

            if (VERBOSE) standard_form_problem.log();

            // Generate and solve artificial problem
            {
                // Create copy of standard form problem to create artificial problem
                Simplex artificial_problem = standard_form_problem;

                std::cout << "Generating artificial problem ...";
                artificial_problem.process_to_artificial_problem();
                std::cout << " done." << std::endl;

                if (VERBOSE) artificial_problem.log();

                // Use artificial problem suggested base to solve it
                std::cout << "Solving artificial problem ..." << std::endl;
                artificial_problem.solve_with_base( artificial_problem.suggested_base );
                std::cout << "Done." << std::endl;

                if (artificial_problem.solution_value != 0) {

                    std::cout << "Problem has no solution." << std::endl;
                    overconstrained = true;
                    return;

                } else {

                    overconstrained = false;

                    if (VERBOSE) std::cout << "Suggested initial base for original problem:";
                    if (VERBOSE) artificial_problem.current_base.log(" ");

                    // If initial base doesn't contain artificial variables
                    // I can just use it, otherwise it may contain an artificial
                    // variable.

                    // Check for existence of a column index related  to an artificial
                    // variable by reading costs vector
                    int artificial_variable = -1;

                    for (int i = 0; i < artificial_problem.solution_dimension; ++i)
                        if ( artificial_problem.objective_function.costs(i) == 1 && artificial_problem.current_base.contains(i) )
                            artificial_variable = i;

                    // If index is still -1 (no artificial variables)
                    if (artificial_variable == -1) {

                        std::cout << "Base is clear about artificial variables, proceed ..." << std::endl;
                        standard_form_problem.suggested_base = artificial_problem.current_base;

                    } else {

                        /*
                          If an artificial variable exists ... I can change the i (artificial)
                          column with a j column in current_out_of base so that:

                          *   j is not an auxiliary variable
                          *   (B^-1)_q * A^j != 0
                          */

                        if (VERBOSE) std::cout << "Artificial variable detected in base: " << artificial_variable << std::endl;
                        int q = artificial_problem.current_base.index_of(artificial_variable);
                        RowVector<Scalar> bi_row_q(artificial_problem.current_base.size());

                        for (unsigned int k = 0; k < artificial_problem.current_base.size(); ++k )
                            bi_row_q(k) = artificial_problem.base_inverse(q,k);

                        // Find j
                        int j = -1;
                        for (unsigned int i = 0; i < standard_form_problem.current_out_of_base.size() && j == -1; ++i) {

                            // Pick the ones that doesn't refer to an artificial variable
                            if ( artificial_problem.costs(i) == 0 ) {
                                ColumnVector<Scalar> column_j(standard_form_problem.current_base.size());

                                for (unsigned int k = 0; k < standard_form_problem.current_base.size(); ++k )
                                    column_j(k) = artificial_problem.coefficients_matrix(k,i);

                                if ((double)( bi_row_q * column_j ) != 0)
                                    j = i;
                            }
                        }

                        if ( j != -1 ) {

                            // Found a j, substitute artificial_value with j
                            standard_form_problem.suggested_base = artificial_problem.current_base;
                            standard_form_problem.suggested_base.substitute(artificial_variable, j);
                            if (VERBOSE) standard_form_problem.suggested_base.log("Now initial base is");

                        } else {

                            /*
                              I didn't find a j which respected the requirements.
                              It may happen that for each j we have (B^-1)_q * A^j = 0,
                              this means that the rows of A are linearly dependent and
                              we can eliminate one of them. Let d be

                              d = e_q * B^-1

                              We have to eliminate a row for which d is non-zero.
                            */

                            std::cout << "Constraints are linearly dependent!" << std::endl;

                            // Find a constraint to eliminate (change)
                            int change = -1;
                            for (unsigned int i = 0; i < standard_form_problem.constraints.size() && change == -1; ++i)
                                if ( bi_row_q(i) != 0 )
                                    change = i;

                            std::cout << "Constraint #" << change << " must be eliminated." << std::endl;
                            has_to_be_fixed = true;
                            return;
                        }
                    }
                }
            }

            std::cout << "Solving problem ..." << std::endl;
            standard_form_problem.solve_with_base(standard_form_problem.suggested_base);
            std::cout << "Done." << std::endl;

            /*
              The solution of the standard form problem must be transposed to
              the original problem.
            */

            std::cout << "Processing standard form solution ..." << std::endl;

            if ( standard_form_problem.unlimited ) {
                unlimited = true;
            } else {
                unlimited = false;
                solution.resize(solution_dimension,1);

                int index = 0;
                for(auto it = standard_form_problem.variables.begin(); it != standard_form_problem.variables.end(); ++it, ++index)
                    (*it)->process(standard_form_problem.solution, solution, index);

                solution_value = standard_form_problem.solution_value;
                dual_variables = standard_form_problem.dual_variables;
                constraints_vector = standard_form_problem.constraints_vector;
                changed_sign = standard_form_problem.changed_sign;
            }
            std::cout << "Done." << std::endl;
        }

        // Print
        void print_solution() const {

            std::cout << "Optimal solution is:" << std::endl;
            for (int i = 0; i < solution_dimension; ++i)
                std::cout << variables.at(i)->name << ":\t\t\t" << solution(i) << std::endl;

            std::cout << std::endl;
            std::cout << "Solution value/cost:\t\t" << solution_value << std::endl;

            Scalar dual_problem_value = (dual_variables * constraints_vector);
            if (changed_sign)
                dual_problem_value *= -1;

            std::cout << "Dual problem value:\t\t" << dual_problem_value << std::endl;

            return;
        }


        void log() const {

            // Title
            for (unsigned int i = 0; i < name.length(); ++i)
                std::cout << "=";
            std::cout << std::endl;
            std::cout << name << std::endl;
            for (unsigned int i = 0; i < name.length(); ++i)
                std::cout << "=";
            std::cout << std::endl;


            // Regular
            std::cout << std::endl;
            std::cout << "Number of regular constraints: " << constraints.size() << std::endl;
            for (auto it = constraints.begin(); it != constraints.end(); ++it)
                it->log();

            // Non negativity
            std::cout << std::endl;
            std::cout << "Number of non-negativity constraints: " << nn_constraints.size() << std::endl;
            for (auto it = nn_constraints.begin(); it != nn_constraints.end(); ++it)
                it->log();

            std::cout << std::endl;
            std::cout << "Objective function:" << std::endl;
            objective_function.log();
            std::cout << std::endl;
        }

        bool has_solutions() const {
            return !overconstrained;
        }

        bool is_unlimited() const {
            return unlimited;
        }

        bool must_be_fixed() const {
            return has_to_be_fixed;
        }

        RowVector<Scalar> const& get_dual_variables() const {
            return dual_variables;
        }


    protected:

        std::string name;

        // Preprocessing
        void process_to_standard_form() {

            // Process non-negative constraints
            int initial_solution_dimension = solution_dimension;
            for ( int i = 0; i < initial_solution_dimension; ++i ) {                        // For each component of x

                bool has_constraint = false;

                // Find an x that doesn't have a non-negativity constraint on it
                for ( auto it = nn_constraints.begin(); it != nn_constraints.end() && !has_constraint; ++it)
                    if ( it->coefficients(i) == 1 )
                        has_constraint = true;

                if ( !has_constraint ) {

                    // Add a non-negativity constraint
                    Matrix<Scalar> eye(1,solution_dimension);
                    eye(i) = 1;
                    this->add_constraint( Constraint<Scalar>( eye, CT_NON_NEGATIVE, 0) );

                    ++solution_dimension;

                    // Add a column to all constraints
                    for ( auto mit = nn_constraints.begin(); mit != nn_constraints.end(); ++mit)
                        mit->add_column(0);


                    // Add another non-negativity constraint
                    Matrix<Scalar> n_eye(1,solution_dimension);
                    n_eye(solution_dimension-1) = 1;

                    this->add_constraint( Constraint<Scalar>( n_eye, CT_NON_NEGATIVE, 0) );

                    // Add a regular constraint
                    for ( auto mit = constraints.begin(); mit != constraints.end(); ++mit)
                        mit->add_column(-mit->coefficients(i));

                    objective_function.add_column(-objective_function.costs(i));

                    // Update variables status
                    std::string aux_name(variables.at(i)->name);
                    Variable<Scalar>* auxiliary = new AuxiliaryVariable<Scalar>(this, (aux_name + "_minus").c_str(), variables.size());
                    Variable<Scalar>* splitted = new SplittedVariable<Scalar>(this, variables.at(i)->name.c_str(), (AuxiliaryVariable<Scalar>*)auxiliary);

                    // Modify variables
                    variables.at(i) = splitted;
                    variables.push_back(auxiliary);
                }
            }

            // Process regular constraints
            for ( auto it = constraints.begin(); it != constraints.end(); ++it ) {

                if ( it->type == CT_MORE_EQUAL )  {

                    // Add empty column to all regular constraints except the current
                    for ( auto mit = constraints.begin(); mit != constraints.end(); ++mit)
                        if ( mit != it )
                            mit->add_column(0);

                    for ( auto mit = nn_constraints.begin(); mit != nn_constraints.end(); ++mit)
                        mit->add_column(0);

                    // Add a 1 column to the current
                    it->add_column(-1);
                    it->type = CT_EQUAL;
                    objective_function.add_column(0);
                    ++solution_dimension;

                    // Add constraint
                    Matrix<Scalar> eye(1,solution_dimension);
                    eye(solution_dimension-1) = 1;
                    this->add_constraint( Constraint<Scalar>( eye, CT_NON_NEGATIVE, 0) );

                    // Update variables vector
                    std::stringstream variable_name;
                    variable_name << "slack_";
                    variable_name << solution_dimension-1;
                    variables.push_back(new SlackVariable<Scalar>(this, variable_name.str().c_str()));


                } else if ( it->type == CT_LESS_EQUAL ) {

                    // Add empty column to all regular constraints except the current
                    for ( auto mit = constraints.begin(); mit != constraints.end(); ++mit)
                        if ( mit != it )
                            mit->add_column(0);

                    for ( auto mit = nn_constraints.begin(); mit != nn_constraints.end(); ++mit)
                        mit->add_column(0);

                    // Add a 1 column to the current
                    it->add_column(1);
                    it->type = CT_EQUAL;
                    objective_function.add_column(0);
                    ++solution_dimension;

                    // Add constraint
                    Matrix<Scalar> eye(1,solution_dimension);
                    eye(solution_dimension-1) = 1;
                    this->add_constraint( Constraint<Scalar>( eye, CT_NON_NEGATIVE, 0) );

                    // Update variables vector
                    std::stringstream variable_name;
                    variable_name << "slack_";
                    variable_name << solution_dimension-1;
                    variables.push_back(new SlackVariable<Scalar>(this, variable_name.str().c_str()));

                }
            }

            // Manipulate objective function
            if ( objective_function.type == OFT_MAXIMIZE ) {
                objective_function.type = OFT_MINIMIZE;
                changed_sign = true;
                objective_function.costs = Matrix<Scalar>::Zero(1,solution_dimension) - objective_function.costs;
            }

            // Update name
            name +=  " (Standard Form)";
        }

        void process_to_artificial_problem() {

            ColumnSet<Scalar> identity;

            // Scans all the columns, when I find a column that is an eye for i
            // put it in the base at position i

            for (unsigned int i = 0; i < constraints.size(); ++i) {

                if (VERBOSE) std::cout << std::endl;
                if (VERBOSE) std::cout << "Checking for column " << i << " of identity." << std::endl;

                bool column_not_found = true;

                for ( int c = solution_dimension-1; c > -1 && column_not_found; --c) {

                    if (VERBOSE) std::cout << "Checking against column " << c << std::endl;

                    bool column_match = true;

                    for ( unsigned int j = 0; j < constraints.size() && column_match; ++j ) {

                        column_match = true;

                        if ( ( i == j && ( constraints.at(j).coefficients(c) != 1 ) ) ||
                             ( i != j && constraints.at(j).coefficients(c) != 0 )  )
                            column_match = false;

                    }

                    if (column_match) {
                        if (VERBOSE) std::cout << "Row match." << std::endl;
                        identity.insert(c);
                        column_not_found = false;
                    }
                }

                if (column_not_found) {
                    if (VERBOSE) std::cout << "Column not found, added artificial variable." << std::endl;
                    identity.insert(-1);
                } else {
                    if (VERBOSE) std::cout << "Column found and added to identity." << std::endl;
                }

            }

            // If artificial variables are needed
            //objective_function.costs.empty();
            objective_function.costs.setZero();

            if ( identity.contains(-1)) {

                for ( unsigned int i = 0; i < identity.size(); ++i) {

                    if ( identity.column(i) == -1 ) {

                        // Add column 1 to constraint i
                        for (unsigned int k = 0; k < constraints.size(); ++k)
                            if ( k == i )
                                constraints.at(k).add_column(1);
                            else
                                constraints.at(k).add_column(0);

                        // Solution vector is bigger
                        identity.column(i) = solution_dimension;
                        ++solution_dimension;

                        // Create non-negative constraint for new variable
                        Matrix<Scalar> eye(1,solution_dimension);
                        eye(solution_dimension-1) = 1;

                        for (unsigned int k = 0; k < nn_constraints.size(); ++k)
                            nn_constraints.at(k).add_column(0);

                        // Objective function costs updated
                        objective_function.add_column(1);

                        this->add_constraint( Constraint<Scalar>( eye, CT_NON_NEGATIVE, 0) );

                    }
                }
            }
            suggested_base = identity;

            name += " --> (Artificial problem)";
        }


        // Solving
        void solve_with_base( ColumnSet<Scalar> const& initial_base ) {
            
            // Preprocess constraints data to lead to matrices
            coefficients_matrix.resize(                                             // A
                constraints.size(),
                solution_dimension
                );

            constraints_vector.resize( constraints.size(), 1);

            for ( unsigned int i = 0; i < constraints.size(); ++i) {
                // Set b
                constraints_vector(i) = constraints.at(i).value;

                for ( int j = 0; j < solution_dimension; ++j)
                    coefficients_matrix(i,j) = constraints.at(i).coefficients(j);
            }

            // Copy costs
            costs = objective_function.costs;

            // Initialize algorithm
            int step = 0;

            if ( initial_base.size() == constraints.size() )
                current_base = initial_base;
            else
                throw(DataMismatchException("Wrong initial base size!"));

            // Exported variables
            optimal = false;
            unlimited = false;



            while ( !optimal && !unlimited ) {

                // Temporary matrices
                RowVector<Scalar>  u;                                                // c_b * B^-1
                RowVector<Scalar>  base_costs (current_base.size());                 // Costs of base

                // Populate current_out_of_base
                current_out_of_base.columns.clear();
                for (int i = 0; i < solution_dimension; ++i)
                    if (!current_base.contains(i))
                        current_out_of_base.insert(i);

                // Every inverse_recalculation steps recompute inverse from scratch
                if ( step % inverse_recalculation_rate == 0 ) {

                    Matrix<Scalar>  base_matrix(  current_base.size(), current_base.size() );

                    // Unpack current base and objective costs
                    for (unsigned int j = 0; j < current_base.size(); ++j) {
                        base_costs(j) = costs( current_base.column(j) );

                        for (unsigned int i = 0; i < current_base.size(); ++i)
                            base_matrix(i,j) = coefficients_matrix(i, current_base.column(j));
                    }

                    // Compute inverse
                    base_inverse = base_matrix.inverse();

                } else {

                    // Unpack objective costs
                    for (unsigned int j = 0; j < current_base.size(); ++j) {
                        base_costs(j) = costs( current_base.column(j) );
                    }

                    Matrix<Scalar> old_inverse = base_inverse;

                    // Compute inverse
                    base_inverse = update_inverse_after_column_change<Scalar>(old_inverse, column_p, old_column);

                }

                if (VERBOSE) std::cout << "Step: " << step << std::endl;
                ++step;



                if (VERBOSE) current_base.log("Columns in base: ");
                if (VERBOSE) current_out_of_base.log("Out of base: ");


                if (VERBOSE) log_matrix(base_inverse, "Base inverse is:");

                // Compute x_B = B^-1 * b
                base_solution = base_inverse * constraints_vector;

                // Compute u = c_B * A;
                u = base_costs * base_inverse;

                if (VERBOSE) log_matrix(u, "U");

                // Compute reduced cost
                reduced_cost = costs - (u * coefficients_matrix);

                if (VERBOSE) log_matrix(reduced_cost, "Current reduced cost is");

                optimal = more_equal_than(reduced_cost, 0.L, TOL);

                bool degenerate = false;
                for (unsigned int i = 0; i < current_base.size() && degenerate == false; ++i)
                    if ( tol_equal(base_solution(i), 0.L, TOL) )
                        degenerate = true;

                if (!optimal) {

                    if (VERBOSE) std::cout << "Base not optimal since reduced cost is negative." << std::endl;

                    // Column of reduced cost with min value (one of the policies)
                    int p = -1;
                    column_p.resize(constraints.size(), 1);
                    Matrix<Scalar> a_tilde;

                    // Bland's strategy
                    for (unsigned int i = 0; i < current_out_of_base.size() && p == -1; ++i)
                        if ( reduced_cost(current_out_of_base.column(i)) < 0 )
                            p = current_out_of_base.column(i);

                    for (unsigned int i = 0; i < constraints.size(); ++i)
                        column_p(i) = coefficients_matrix(i,p);

                    if (VERBOSE) std::cout << "The column to insert is " << p << std::endl;
                    if (VERBOSE) log_matrix(column_p, "That is ...");

                    // Compute a_tilde
                    a_tilde = base_inverse * column_p;

                    if (VERBOSE) log_matrix(a_tilde, "a_tilde");

                    unlimited = less_equal_than(a_tilde, 0.L, TOL);

                    if (!unlimited) {
                        if (VERBOSE) std::cout << "Problem not unlimited." << std::endl;

                        // Bland's strategy
                        int q_position = -1;
                        for (unsigned int i = 0; i < current_base.size(); ++i) {

                            Scalar value = base_solution(i) / a_tilde(i);

                            if ( a_tilde(i) > 0 &&
                                 ( q_position == -1 || value < ( base_solution(q_position) / a_tilde(q_position) ) ) )
                                q_position = i;

                        }


                        int q = current_base.column(q_position);
                        old_column = q_position;

                        if (VERBOSE) std::cout << "The column to take off is " << q <<  std::endl;

                        // Take off q, push in p
                        current_base.substitute(q,p);

                    } else {
                        std::cout << "Problem unlimited." << std::endl;
                    }

                } else {
                    std::cout << "Optimal found at step " << step << "." << std::endl;
                    RowVector<Scalar> objective_function_base(current_base.size());
                    ColumnVector<Scalar> full_solution(solution_dimension);

                    // Update dual variables
                    dual_variables = u;

                    for (unsigned int i = 0; i < constraints.size(); ++i)
                        objective_function_base(i) = costs(current_base.column(i));

                    for ( int i = 0; i < solution_dimension; ++i)
                        if ( current_base.contains(i) )
                            full_solution(i) = base_solution( current_base.index_of( i ) );

                    if (VERBOSE) log_matrix(full_solution, "Solution:");

                    // Saves some flops
                    solution_value = (objective_function_base * base_solution);

                    if ( changed_sign )
                        solution_value = -solution_value;

                    if (VERBOSE) std::cout << "Solution value: " << solution_value << std::endl;

                    solution = full_solution;

                }
            }
        }


        // Column sets
        ColumnSet<Scalar> suggested_base;
        ColumnSet<Scalar> current_base;
        ColumnSet<Scalar> current_out_of_base;

        // Data
        ObjectiveFunction<Scalar> objective_function;
        std::vector< Constraint<Scalar> > constraints;
        std::vector< Constraint<Scalar> > nn_constraints;
        std::vector< Variable<Scalar>* > variables;

        // Processed data
        RowVector<Scalar> costs;
        Matrix<Scalar> coefficients_matrix;
        ColumnVector<Scalar> constraints_vector;
        Matrix<Scalar> base_inverse;
        RowVector<Scalar> dual_variables;
        Matrix<Scalar> column_p;
        int solution_dimension, old_column;

        // Results
        ColumnVector<Scalar> base_solution;
        Matrix<Scalar> solution;
        Matrix<Scalar> reduced_cost;
        Scalar solution_value;

        bool    optimal,
            unlimited,
            overconstrained,
            has_to_be_fixed,
            changed_sign;

        int inverse_recalculation_rate;

    };

}

#endif
