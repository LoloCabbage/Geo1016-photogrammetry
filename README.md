# calibration.cpp
  Calls the Calibration::calibration function of the calibtaion_method.cpp file. 
  In this file they insert the points saved in the text file to the calibration function. And it inserts the other parameters, which needs to be updates in the calibration::calibration function.
  --> I do not know what happens when the calibration is a succes??
  
# calibtaion_method.cpp
TO DO
1. check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
     --> They already check this in line 204 of calibration.cpp right?
3. construct the P matrix (so P * m = 0)
4. solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition. Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point should be very close to your input images points.
5. extract intrinsic parameters from M.
6. extract extrinsic parameters from M.

# Instruction (to be removed) inside the code
    std::cout << "\nTODO: I am going to implement the calibration() function in the following file:\n"
                 "\t    - calibration_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tCamera calibration requires computing the SVD and inverse of matrices.\n"
                 "\tIn this assignment, I provide you with a 'Matrix' and a 'Vector' data structures for storing and\n"
                 "\tmanipulating matrices and vectors of arbitrary sizes. I also wrote some code to show you how to:\n"
                 "\t    - compute the SVD of a matrix;\n"
                 "\t    - compute the inverse of a matrix;\n"
                 "\t    - compute the transpose of a matrix.\n\n"
                 "\tFeel free to use any of the provided data structures and functions. The commonly used linear algebra\n"
                 "\tfunctions are provided in the following files:\n"
                 "\t    - Calibration/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;
