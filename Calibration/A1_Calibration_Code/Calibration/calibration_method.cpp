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
    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    bool valid = check_input(points_3d, points_2d);
    // TODO: construct the P matrix (so P * m = 0).

    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    // TODO: extract intrinsic parameters from M.

    // TODO: extract extrinsic parameters from M.

    return false;
}

















