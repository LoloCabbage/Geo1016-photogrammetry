# Triangulation.cpp
Finish this function for reconstructing 3D geometry from corresponding image points.
# To Do List
- [x] 1. check if the input is valid
- [x] 2. Estimate fundamental matrix F
- [x] 3. Recover relative pose (R and t)
- [x] 4. Determine the 3D coordinates 
# Note
Don't forget to
- write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
- write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera, which can help you check if R and t are correct).
You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the viewer will be notified to visualize the 3D points and update the view).

There are a few cases you should return 'false' instead, for example:
- function not implemented yet;
- input not valid (e.g., not enough points, point numbers don't match);
- encountered failure in any step.
# Related functions and example code given by the teacher
**Functions for linear method**

*Triangulation/matrix.h*  Matrices of arbitrary dimensions and related functions

*Triangulation/vector.h*  Vectors of arbitrary dimensions and related functions

*Triangulation/matrix_algo.h*  Determinant, inverse, SVD, linear least-squares"

**Functions for non-linear method (optional task)**

*Tutorial_NonlinearLeastSquares/main.cpp* for an example and some explanations.

**Example codes**
    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4
