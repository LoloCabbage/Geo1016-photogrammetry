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
