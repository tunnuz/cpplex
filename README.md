Introduction
============

At the very current status the project is composed of two chunks:

    *   Simplex (The real simplex algorithm implementation)
    *   The solver

The solver depends on Simplex. The project can 
be built either under Windows (there is a MS Visual Studio 2008 solution which
can be automatically coverted to a 2010 solution without complications), Linux 
(make) or MacOS (make).

The project is currently working acceptably, however some optimizations and 
improvements could be added in the future:

    *   Optimization to exploit upper and lower limitations in variables
    *   Branch-and-bound for integer linear programming
	*	Multi-threading (for branch and bound?)
	*   Simplex compiled as libraries under MacOS and Linux
	*   Better language to define problems
	*   Doxygen documentation

Basically, anything that can improve performance and usability is a good
feature candidate.

How to use it
=============

All you have to do to compile the project is:

    make
    
If you want to clean up all the dependencies files and object files just type:

    make wipe
    
The "lib" directory was used in a previous Linux version to compile
Simplex as a library, but was discontinued because I wanted the Makefile to
work under MacOS as well and I had not enough time to learn how to write
portable Makefiles. "lib" is currently used by the Visual Studio solution to
compile Simplex as a static library. Once you have built the project,
to execute the solver with a problem file (look for some example problems in 
the "problems" directory) run under MacOS or Linux:

    ./solver <name_of_the_problem_file>
    
Under Windows:

    solver.exe <name_of_the_problem_file>
    
Problems have a very straightforward structure, therefore you can learn the
syntax from the examples in the "problems" or from the Google Code page.

License
=======

See LICENSE file.

Drop a line
===========

If you use this code for your software, please let me know with a mail
message at tunnuz@gmail.com, or not.
