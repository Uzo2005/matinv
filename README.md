We are currently doing linear algebra as a course and so I do lot of boring-ass calculations involving matrices. In this repo, I have set out to automate all the boring calculations, so that I can spend less time on my assignments.
#### Usage
    The intended use-case for this project is quick verification of hand-computed results of linear algebra problem sets and theorems. To be able to do this and maintain accuracy, I entirely avoided using floating point types and instead implemented `Fraction` as a numeric type. Using fractions to compute arithmetic operations ensured that every result obtained maintained 100% accuracy and that no floating-point error could hinder exploration.

    For operations which involve row operations, calling `showSteps` on the result of these functions would print the necessary row operation steps needed to arrive at the final values.

    I would strongly advise readers not to use code from this repository in their own applications because I have sacrificed speed of execution for the fun of exploration.


#### So far I have implemented:

    1. Determinants Of Matrices
    2. Adjoint(Adjucate?) Of Matrices
    3. Inverse Of Matrices (will throw exception for singular matrices)
    4. Normal and Reduced Echelon Forms of Matrix, and hence verification of consistency for a linear system of equations
    5. LU Factorisation of matrices

#### I hope to also implement:

    1. Eigenvalues and Eigenvectors
    2. Diagonilization of Matrices

Also if I have time, I will make this a website so other people can use it.
