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
  #define P_X 0
  #define P_Y 1
  #define P_Z 2
  #define E_R 3
  #define E_P 4
  #define E_Y 5
  #define V_X 6
  #define V_Y 7
  #define V_Z 8
  #define DE_R 9
  #define DE_P 10
  #define DE_Y 11
  #define Q_1 6
  #define Q_2 7
  #define Q_3 8
  #define U_1 0
  #define U_2 1
  #define U_3 2
  #define U_4 3
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
    VectorXd x_vec_; // size 12, s, ds
    VectorXd q_vec_; // size 9, q, dq, ddq
    VectorXd u_vec_;
    MatrixXd Ds_;
    MatrixXd Ds_inv_;
    MatrixXd C_;
    VectorXd Bs_;
    VectorXd gs_;
    MatrixXd Ds3_;
    VectorXd dds2_;
    // Ds_x Bs_x Bs_u Cs_x Cs_dx gs_x Ds3_x Cs3_x Cs3_dx
    std::vector<MatrixXd> D_dx_vec_;
    std::vector<MatrixXd> D_ddx_vec_;
    std::vector<MatrixXd> C_dx_vec_; // d x
    std::vector<MatrixXd> C_d_dx_vec_; // d (dx)
    std::vector<VectorXd> Bs_du_vec_; // d f1, f2, f3, f4
    std::vector<VectorXd> Bs_dx_vec_; // d er, ep, eq
    std::vector<VectorXd> gs_dx_vec_; // d er, ep, eq

    // mid result
    MatrixXd R_local_;
    std::vector<MatrixXd> R_local_dx_vec_; // d er, ep, eq
    std::vector<MatrixXd> R_local_ddx_vec_; // d er, ep, eq
    std::vector<MatrixXd> R_link_local_vec_;
    std::vector<MatrixXd> R_link_local_dx_vec_; // (d R_l0 q1, q2, q3), ..., (d R_l3 q1, q2, q3
    std::vector<MatrixXd> R_link_local_ddx_vec_; // (d R_l0 q1, q2, q3), ..., (d R_l3 q1, q2, q3)
    MatrixXd T_local_;
    std::vector<MatrixXd> T_local_dx_vec_; // d er, ep, eq
    std::vector<MatrixXd> T_local_ddx_vec_; // d er, ep, eq
    MatrixXd Q_local_;
    std::vector<MatrixXd> Q_local_dx_vec_; // d er, ep, eq
    std::vector<MatrixXd> Q_local_ddx_vec_; // d er, ep, eq
    MatrixXd link_center_pos_local_;
    std::vector<MatrixXd> link_center_pos_local_dx_vec_; // d q1, q2, q3
    std::vector<MatrixXd> link_center_pos_local_ddx_vec_; // d q1, q2, q3
    std::vector<MatrixXd> Jacobian_P_vec_;
    std::vector<MatrixXd> Jacobian_P_dq_vec_; // (JacoP[0] d q1, q2, q3), ..., (JacoP[3] d q1, q2, q3)
    std::vector<MatrixXd> Jacobian_P_ddq_vec_; // (JacoP[0] d q1, q2, q3), ..., (JacoP[3] d q1, q2, q3)
    std::vector<MatrixXd> Jacobian_W_vec_;
    MatrixXd S_operation_result_;
    std::vector<MatrixXd> S_operation_dx_vec_; // d er, ep, eq, q1, q2, q3
    std::vector<MatrixXd> S_operation_ddx_vec_; // d er, ep, eq, q1, q2, q3
    MatrixXd D11_;
    MatrixXd D12_;
    MatrixXd D13_;
    MatrixXd D22_;
    MatrixXd D23_;
    MatrixXd D33_;

    void updateHrydrusDynamicParamater(VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr);
    VectorXd getStateDerivative();
    void updateMiddleVariable(VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr);
    MatrixXd vectorToSkewMatrix(VectorXd s);
    void updateMainMatrix();
    void linaerizeState(MatrixXd *s_mat_ptr, MatrixXd *u_mat_ptr);
  };
}


#endif
