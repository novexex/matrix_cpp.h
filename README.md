# s21_matrix+

Implementation of the s21_matrix_oop.h library.

## Implementation of the matrix_cpp.h library functions

- The program must be developed in C++ language of C++17 standard using gcc compiler
- The program code must be located in the src folder
- When writing code it is necessary to follow the Google style
- Implement the matrix as an `S21Matrix` class
- Use only the `matrix`, `rows` and `columns` fields as private.
- Implement the access to private fields `rows` and `columns` via accessor and mutator. If the matrix increases in size, it is filled with zeros. If it decreases in size, the excess is simply discarded.
- Make it as a static library (with s21_matrix_oop.h header file)
- Implement the operations described [above](#matrix-operations)
- Overload the operators according to the table in the chapter [above](#matrix-operations)
- Prepare full coverage of library functions code with unit-tests using the GTest library
- Provide a Makefile for building the library and tests (with targets all, clean, test, s21_matrix_oop.a)

