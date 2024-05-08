/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"
#include <cmath>
#include "vector.h"

using namespace easy3d;

/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by fx, fy, cx, cy, skew, R, and t).
 */

//// check if input is valid
bool check_input(const std::vector<Vector3D>& points_3d, const std::vector<Vector2D>& points_2d) {
    // check if the number of correspondences >= 6
    if (points_3d.size() < 6 || points_2d.size() < 6) {
        std::cerr << "Error: the number of correspondences must be at least 6." << std::endl;
        return false;
    }

    // check if the sizes of 2D/3D points match
    else if (points_3d.size() != points_2d.size()) {
        std::cerr << "Error: the sizes of 2D/3D points must match." << std::endl;
        return false;
    }

    std::cout << "Input data is valid." << std::endl;
    return true;
}

Matrix construct_P(const std::vector<Vector3D>& points_3d, const std::vector<Vector2D>& points_2d, int size_input_points){
    Matrix P(size_input_points * 2, 12);
    P.load_zero();

    for (int i = 0; i < points_3d.size(); ++i) {
        int i_in_P = i * 2;

        double Xi = points_3d[i][0];
        double Yi = points_3d[i][1];
        double Zi = points_3d[i][2];
        double ui = points_2d[i][0];
        double vi = points_2d[i][1];

        P[i_in_P][0] = Xi;
        P[i_in_P][1] = Yi;
        P[i_in_P][2] = Zi;
        P[i_in_P][3] = 1.0;

        P[i_in_P][8] = -ui * Xi;
        P[i_in_P][9] = -ui * Yi;
        P[i_in_P][10] = -ui * Zi;
        P[i_in_P][11] = -ui;

        P[i_in_P+1][4] = Xi;
        P[i_in_P+1][5] = Yi;
        P[i_in_P+1][6] = Zi;
        P[i_in_P+1][7] = 1.0;

        P[i_in_P+1][8] = -vi * Xi;
        P[i_in_P+1][9] = -vi * Yi;
        P[i_in_P+1][10] = -vi * Zi;
        P[i_in_P+1][11] = -vi;

    }
    return P;
}

Matrix construct_m(Matrix &P, int size_input_points){
    Matrix U(size_input_points*2, size_input_points*2);
    Matrix S(size_input_points*2, 12);
    Matrix V(12, 12);
    svd_decompose(P, U, S, V);
    int numCols = V.cols();
    Vector m = V.get_column(numCols - 1);
    Matrix M(3, 4, m.data());
    return M;
}

void proj_2D (const Matrix &M,const std::vector<Vector3D>& points_3d, const std::vector<Vector2D>& points_2d){
    Matrix homo_matrix(4,points_3d.size());
    for (int i = 0; i < points_3d.size(); ++i) {
        double x = points_3d[i][0];
        double y = points_3d[i][1];
        double z = points_3d[i][2];

        homo_matrix[0][i] = x;
        homo_matrix[1][i] = y;
        homo_matrix[2][i] = z;
        homo_matrix[3][i] = 1;
    }

    Matrix M_homo = M * homo_matrix;
    Matrix proj_2D(2,points_3d.size());
    double rmse = 0;
    for (int i = 0; i < points_3d.size(); ++i) {
        proj_2D[0][i] = M_homo[0][i]/M_homo[2][i];
        proj_2D[1][i] = M_homo[1][i]/M_homo[2][i];
        rmse = rmse + (proj_2D[0][i]-points_2d[i][0])*(proj_2D[0][i]-points_2d[i][0]);
        rmse = rmse + (proj_2D[1][i]-points_2d[i][1])*(proj_2D[1][i]-points_2d[i][1]);
    }
    rmse = sqrt(rmse/points_3d.size()/2);
    std::cout << "The projected 2D coordinates" << proj_2D << std::endl;
    std::cout << "The rmse of projected data is :"<< rmse << std::endl;
}

void extracting_parameters(Matrix &M, double& fx, double& fy, double& cx, double& cy, double& s, Matrix33& R, Vector3D& t){
    Vector3D a1 = {M[0][0], M[0][1], M[0][2]};
    Vector3D a2 = {M[1][0], M[1][1], M[1][2]};
    Vector3D a3 = {M[2][0], M[2][1], M[2][2]};
    Vector3D b = {M[0][3], M[1][3], M[2][3]};

    double rho = 1 / length(a3);
    double theta = acos(-
                                dot(cross(a1,a3), cross(a2,a3)) /
                        (length(cross(a1, a3)) * length(cross(a2, a3))));
    double alpha = rho * rho * length(cross(a1, a3)) * sin(theta);
    double beta = rho * rho * length(cross(a2, a3))* sin(theta);

    fx = alpha;
    fy = beta / sin(theta);
    s = - alpha / tan(theta);
    cx = rho * rho * dot (a1, a3);
    cy = rho * rho * dot (a2, a3);

    Vector3D r1 = cross(a2,a3) / length(cross(a2,a3));
    Vector3D r3 = rho * a3;
    Vector3D r2 = cross(r3, r1);
    R.set_row(0, r1);
    R.set_row(1,r2);
    R.set_row(2,r3);

    Matrix33 K;
    K.set(0,0,fx);
    K.set(0,1,s);
    K.set(0,2,cx);
    K.set(1,1,fy);
    K.set(1,2,cy);
    K.set(2,2,1);

    t =  rho * inverse(K) * b;
}

bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,  /// output: focal length (i.e., K[0][0]).
        double& fy,  /// output: focal length (i.e., K[1][1]).
        double& cx,  /// output: x component of the principal point (i.e., K[0][2]).
        double& cy,  /// output: y component of the principal point (i.e., K[1][2]).
        double& s,   /// output: skew factor (i.e., K[0][1]), which is s = -alpha * cot(theta).
        Matrix33& R, /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t) /// output：a 3D vector encoding camera translation.
{

  // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
  bool valid = check_input(points_3d, points_2d);
  if (!valid){
      return false;
  }
  // TODO: construct the P matrix (so P * m = 0).
  int size_input_points = points_3d.size();
  Matrix P = construct_P(points_3d, points_2d, size_input_points);

  // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
  Matrix M = construct_m(P, size_input_points);

  /// Optional: you can check if your M is correct by applying M on the 3D points.
  /// If correct, the projected point should be very close to your input images points.
  proj_2D(M, points_3d, points_2d);

  // TODO: extract intrinsic parameters from M.
  // TODO: extract extrinsic parameters from M.
  extracting_parameters(M, fx, fy,cx, cy, s, R, t);

  return true;
}


















