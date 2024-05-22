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

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */

//// struct for normalization information
struct normalization {
    Vector2D centroid; // image center
    double avg_dist; // average distance to the image center
    std::vector<Vector2D> normal_points; // normalized image points
    Matrix T; // transformation matrix
};

//// check if the input is valid
bool check_input(const std::vector<Vector2D>& points_0, const std::vector<Vector2D>& points_1) {
    // check if the number of correspondences >= 8
    if (points_0.size() < 8 || points_1.size() < 8) {
        std::cerr << "Error: the number of correspondences must be at least 6." << std::endl;
        return false;
    }

        // check if the sizes of 2D/3D points match
    else if (points_0.size() != points_1.size()) {
        std::cerr << "Error: the sizes of 2D/3D points must match." << std::endl;
        return false;
    }

    std::cout << "Input data is valid." << std::endl;
    return true;
}

//// normalize the points
normalization normalize(const std::vector<Vector2D>& points){
    std::vector<Vector2D> normal_points;
    // compute the centroid of the points
    Vector2D centroid(0, 0);
    for (const auto& point : points) {
        centroid += point;
    }
    centroid /= points.size();
    // use the centroid as the origin
    for (const auto& point : points) {
        normal_points.push_back(point - centroid);
    }
    // compute the average distance to the origin
    double avg_dist = 0;
    for (const auto& point : normal_points) {
        avg_dist += point.norm();
    }
    avg_dist /= normal_points.size();
    // scale the points so that the average distance to the origin is sqrt(2)
    for (auto& point : normal_points) {
        point *= sqrt(2) / avg_dist;
    }

    // compute the transformation matrix
    Matrix T(3, 3);
    T.set(0, 0, sqrt(2) / avg_dist);
    T.set(1, 1, sqrt(2) / avg_dist);
    T.set(2, 2, 1);
    T.set(0, 2, -centroid.x() * sqrt(2) / avg_dist);
    T.set(1, 2, -centroid.y() * sqrt(2) / avg_dist);

    std::cout << "T: " << T << std::endl;

    normalization norm_data;
    norm_data.centroid = centroid;
    norm_data.avg_dist = avg_dist;
    norm_data.normal_points = normal_points;
    norm_data.T = T;
    return norm_data;
}

//// generate the fundamental matrix
Matrix F (const std::vector<Vector2D>& normal_points_0, const std::vector<Vector2D>& normal_points_1){
    int n = normal_points_0.size();
    Matrix A(n, 9);
    for (int i = 0; i < n; i++) {
        double x0 = normal_points_0[i].x();
        double y0 = normal_points_0[i].y();
        double x1 = normal_points_1[i].x();
        double y1 = normal_points_1[i].y();

        A(i, 0) = x0 * x1;
        A(i, 1) = y0 * x1;
        A(i, 2) = x1;
        A(i, 3) = x0 * y1;
        A(i, 4) = y0 * y1;
        A(i, 5) = y1;
        A(i, 6) = x0;
        A(i, 7) = y0;
        A(i, 8) = 1;
    }

    // svd decomposition
    Matrix U(n, n);
    Matrix D(n, 9);
    Matrix V(9, 9);
    svd_decompose(A, U, D, V);
    int numCols = V.cols();
    Vector F = V.get_column(numCols - 1);
    Matrix F_matrix(3, 3, F.data());

    // enforce rank 2
    Matrix U_F(3, 3);
    Matrix D_F(3, 3);
    Matrix V_F(3, 3);
    svd_decompose(F_matrix, U_F, D_F, V_F);
    D_F.set(2, 2, 0);
    F_matrix = U_F * D_F * V_F.transpose(); //// TRANSPOSE OR NOT?? to CHECK
    return F_matrix;
}

//// denormalize the fundamental matrix
Matrix denormalize(const Matrix& F, const Matrix& T0, const Matrix& T1){
    return T1.transpose() * F * T0;
}

//// Compute essential matrix E
Matrix compute_matrix_E(double fx, double fy, double cx, double cy, double s, Matrix &F){
    Matrix33 K;
    K(0,0) = fx;
    K(0,1) = s;
    K(0,2) = cx;
    K(1,1) = fy;
    K(1,2) = cy;
    K(2,2) = 1;
    Matrix E = K.transpose() * F * K;
    return K;
}

//// Recover the translation vector and rotation matrix
void find_possible_R_and_t(Matrix &E, Matrix33 &R1, Matrix33 &R2, Vector3D t1, Vector3D t2){
    Matrix33 W;
    W(0,1) = -1;
    W(1,0) = 1;
    W(2,2) = 1;

    Matrix33 U;
    Matrix33 D;
    Matrix33 V_transpose;
    svd_decompose(E, U, D, V_transpose);

    R1 = determinant(U * W * V_transpose) * U * W * V_transpose;
    R2 = determinant(U * W.transpose() * V_transpose) * U * W.transpose() * V_transpose;
    t1[0] = U(0,2);
    t1[1] = U(1,2);
    t1[2] = U(2,2);
    t2[0] = U(0,2);
    t2[1] = U(1,2);
    t2[2] = U(2,2);
}


bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        double s,                 /// input: the skew factor (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    // TODO: check if the input is valid (always good because you never known how others will call your function).
    bool valid = check_input(points_0, points_1);
    if (!valid){
        return false;
    }

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    normalization norm_data_0 = normalize(points_0);
    normalization norm_data_1 = normalize(points_1);
    std::vector<Vector2D> normal_points_0 = norm_data_0.normal_points;
    std::vector<Vector2D> normal_points_1 = norm_data_1.normal_points;

    Matrix F_matrix = F(normal_points_0, normal_points_1);
    Matrix F_denormalized = denormalize(F_matrix, norm_data_0.T, norm_data_1.T);

    // TODO: - compute the essential matrix E;
    Matrix E = compute_matrix_E(fx, fy, cx, cy, s, F_denormalized);

    // TODO: - recover rotation R and t.
    Matrix33 R1, R2;
    Vector3D t1, t2;
    find_possible_R_and_t(E, R1, R2, t1, t2);


    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)


    return points_3d.size() > 0;
}