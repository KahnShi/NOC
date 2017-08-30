// -*- mode: c++ -*-
/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2016, JSK Lab
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/o2r other materials provided
 *     with the distribution.
 *   * Neither the name of the JSK Lab nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

#ifndef HYDRUS_DYNAMICS_H
#define HYDRUS_DYNAMICS_H

#include <iostream>
#include <vector>
/* linear algebra */
#include <math.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Eigenvalues>

using namespace Eigen;

namespace hydrus_dynamics{
  class HydrusDynamics{
  public:
    HydrusDynamics(int n_links, double l_length, std::vector<double> *link_weight_vec_ptr, MatrixXd *I_ptr);
    ~HydrusDynamics();
    int n_links_;
    double link_length_;
    std::vector<double> link_weight_vec_;
    double weight_sum_;
    MatrixXd Inertial_;
    /* state */
    double px_;
    double py_;
    double pz_;
    double er_;
    double ep_;
    double ey_;
    double d_px_;
    double d_py_;
    double d_pz_;
    double d_er_;
    double d_ep_;
    double d_ey_;
    /* control */
    double q1_;
    double q2_;
    double q3_;
    double f1_;
    double f2_;
    double f3_;
    double f4_;
    /* matrix */
    VectorXd x_vec_;
    VectorXd q_vec_;
    VectorXd u_vec_;
    MatrixXd Ds_;
    MatrixXd Cs_;
    VectorXd Bs_;
    VectorXd Csds_;
    VectorXd gs_;
    MatrixXd Ds3_;
    MatrixXd Cs3_;
    VectorXd dds2_;
    // Ds_x Bs_x Bs_u Csds_x Csds_dx gs_x Ds3_x Cs3_x Cs3_dx
    MatrixXd Ds_x_;
    MatrixXd Bs_x_;
    MatrixXd Bs_u_;
    MatrixXd Csds_x_;
    MatrixXd Csds_dx_;
    MatrixXd gs_x_;
    MatrixXd Ds3_x_;
    MatrixXd Cs3_x_;
    MatrixXd Cs3_dx_;

    // mid result
    MatrixXd R_local_;
    std::vector<MatrixXd> R_local_d_vec_; // d er, ep, eq
    std::vector<MatrixXd> R_link_local_vec_;
    MatrixXd T_local_;
    std::vector<MatrixXd> T_local_d_vec_; // d er, ep, eq
    MatrixXd link_center_pos_local_;
    std::vector<MatrixXd> link_center_pos_local_d_vec_; // d q1, q2, q3
    std::vector<MatrixXd> Jacobian_P_vec_;
    std::vector<MatrixXd> Jacobian_W_vec_;
    MatrixXd S_operation_result_;
    MatrixXd D11_;
    MatrixXd D12_;
    MatrixXd D13_;
    MatrixXd D22_;
    MatrixXd D23_;
    MatrixXd D33_;

    void getCurrentState(VectorXd *s, VectorXd *q);
    MatrixXd vectorSkewToMatrix(Vector3d s);
    void calculateSkewOperation(std::vector<MatrixXd> *S_vec_ptr);
    void updateMatrixD();
  };
}


#endif
