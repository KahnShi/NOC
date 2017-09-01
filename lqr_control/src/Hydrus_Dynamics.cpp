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
 *   * Redistributions in binary form must rep_roduce the above
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

#include <lqr_control/Hydrus_Dynamics.h>
namespace hydrus_dynamics{
  HydrusDynamics::HydrusDynamics(int n_links, double l_length, std::vector<double> *link_weight_vec_ptr, MatrixXd *I_ptr){
    n_links_ = n_links;
    link_length_ = l_length;
    link_weight_vec_ = (*link_weight_vec_ptr);
    weight_sum_ = 0.0;
    for (int i = 0; i < 4; ++i)
      weight_sum_ += link_weight_vec_[i];

    Inertial_ = *I_ptr;

    // init
    Ds_ = MatrixXd::Zero(6, 6);
    Ds3_ = MatrixXd::Zero(6, 3);
    C_ = MatrixXd::Zero(9, 9);

    for (int i = 0; i < 3; ++i){
      R_local_dx_vec_.push_back(MatrixXd::Zero(3, 3));
      T_local_dx_vec_.push_back(MatrixXd::Zero(3, 3));
      Q_local_dx_vec_.push_back(MatrixXd::Zero(3, 3));
      link_center_pos_local_dx_vec_.push_back(MatrixXd::Zero(3, 4));
      for (int j = 0; j < 3; ++j){
        R_local_ddx_vec_.push_back(MatrixXd::Zero(3, 3));
        T_local_ddx_vec_.push_back(MatrixXd::Zero(3, 3));
        Q_local_ddx_vec_.push_back(MatrixXd::Zero(3, 3));
        link_center_pos_local_ddx_vec_.push_back(MatrixXd::Zero(3, 4));
      }
      Bs_dx_vec_.push_back(VectorXd::Zero(6));
      gs_dx_vec_.push_back(VectorXd::Zero(6));
    }
    for (int i = 0; i < 4; ++i){
      R_link_local_vec_.push_back(MatrixXd::Zero(3, 3));
      Jacobian_P_vec_.push_back(MatrixXd::Zero(3, 3));
      Jacobian_W_vec_.push_back(MatrixXd::Zero(3, 3));
      for (int j = 0; j < 3; ++j){
        Jacobian_P_dq_vec_.push_back(MatrixXd::Zero(3, 3));
        R_link_local_dx_vec_.push_back(MatrixXd::Zero(3, 3));
        for (int k = 0; k < 3; ++k){
          Jacobian_P_ddq_vec_.push_back(MatrixXd::Zero(3, 3));
          R_link_local_ddx_vec_.push_back(MatrixXd::Zero(3, 3));
        }
      }
      Bs_du_vec_.push_back(VectorXd::Zero(6));
    }
    for (int i = 0; i < 6; ++i){
      S_operation_dx_vec_.push_back(MatrixXd::Zero(3, 4));
      for (int j = 0; j < 6; ++j){
        S_operation_ddx_vec_.push_back(MatrixXd::Zero(3, 4));
      }
      C_dx_vec_.push_back(MatrixXd::Zero(9, 9));
      C_d_dx_vec_.push_back(MatrixXd::Zero(9, 9));
    }
    for (int i = 0; i <= Q_3; ++i){
      D_dx_vec_.push_back(MatrixXd::Zero(9, 9));
      for (int j = 0; j <= Q_3; ++j){
        D_ddx_vec_.push_back(MatrixXd::Zero(9, 9));
      }
    }
  }

  HydrusDynamics::~HydrusDynamics(){
  }

  VectorXd HydrusDynamics::getStateDerivative(){
    VectorXd d_x(6);
    for (int i = 0; i < 6; ++i)
      d_x(i) = x_vec_(6 + i);
    VectorXd d_q(3), dd_q(3);
    for (int i = 0; i < 3; ++i){
      d_q(i) = q_vec_(3 + i);
      dd_q(i) = q_vec_(6 + i);
    }
    VectorXd dd_x(6);
    dd_x = Ds_inv_ * (Bs_ - C_.block<6, 6>(0, 0) * d_x - gs_ - Ds3_ * dd_q
                      - C_.block<6, 3>(0, 6) * d_q);
    VectorXd d_state(12);
    for (int i = 0; i < 6; ++i){
      d_state(i) = d_x(i);
      d_state(6 + i) = dd_x(i);
    }
    return d_x;
  }

  void HydrusDynamics::linaerizeState(VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr, MatrixXd *s_mat_ptr, MatrixXd *u_mat_ptr){
    updateMiddleVariable(x_ptr, u_ptr, joint_ptr);
    updateMainMatrix();

    VectorXd d_x(6);
    for (int i = 0; i < 6; ++i)
      d_x(i) = x_vec_(6 + i);
    VectorXd d_q(3), dd_q(3);
    for (int i = 0; i < 3; ++i){
      d_q(i) = q_vec_(3 + i);
      dd_q(i) = q_vec_(6 + i);
    }
    *s_mat_ptr = MatrixXd::Zero(12, 12);
    *u_mat_ptr = MatrixXd::Zero(12, 4);
    // s_mat
    s_mat_ptr->block<6, 6>(0, 6) = MatrixXd::Identity(6, 6);
    for (int i = E_R; i <= E_Y; ++i){
      s_mat_ptr->block<6, 1>(6, i) += -Ds_inv_ * D_dx_vec_[i].block<6, 6>(0, 0)
        * Ds_inv_ * (Bs_ + C_.block<6, 6>(0, 0) * d_x - gs_ - Ds3_ * dd_q
                     - C_.block<6, 3>(0, 6) * d_q)
        + Ds_inv_ * (Bs_dx_vec_[i-E_R] - C_dx_vec_[i].block<6, 6>(0, 0) * d_x
                     - gs_dx_vec_[i-E_R] - D_dx_vec_[i].block<6, 3>(0, 6) * dd_q
                     - C_dx_vec_[i].block<6, 3>(0, 6) * d_q);
    }
    for (int i = V_X; i <= DE_Y; ++i){
      VectorXd d_x_d_x = VectorXd::Zero(6); d_x_d_x(i-V_X) = 1.0;
      s_mat_ptr->block<6, 1>(6, i) += Ds_inv_
        * (- C_.block<6, 6>(0, 0) * d_x_d_x
           - C_d_dx_vec_[i-V_X].block<6, 6>(0, 0) * d_x
           - C_dx_vec_[i-V_X].block<6, 3>(0, 6) * d_q);
    }
    // u_mat
    for (int i = U_1; i <= U_3; ++i){
      u_mat_ptr->block<6, 1>(6, i) += Ds_inv_ * Bs_du_vec_[i];
    }
  }

  void HydrusDynamics::updateMiddleVariable(VectorXd *x_ptr, VectorXd *u_ptr,VectorXd *joint_ptr){
    x_vec_ = *x_ptr;
    u_vec_ = *u_ptr;
    q_vec_ = *joint_ptr;
    px_ = (*x_ptr)[0];
    py_ = (*x_ptr)[1];
    pz_ = (*x_ptr)[2];
    er_ = (*x_ptr)[3];
    ep_ = (*x_ptr)[4];
    ey_ = (*x_ptr)[5];
    d_px_ = (*x_ptr)[6];
    d_py_ = (*x_ptr)[7];
    d_pz_ = (*x_ptr)[8];
    d_er_ = (*x_ptr)[9];
    d_ep_ = (*x_ptr)[10];
    d_ey_ = (*x_ptr)[11];

    q1_ = (*joint_ptr)[0];
    q2_ = (*joint_ptr)[1];
    q3_ = (*joint_ptr)[2];

    // mid result
    // R_local_ and its derivative
    R_local_ << cos(ep_)*cos(ey_), cos(ey_)*sin(ep_)*sin(er_) - cos(er_)*sin(ey_), sin(er_)*sin(ey_) + cos(er_)*cos(ey_)*sin(ep_),
      cos(ep_)*sin(ey_), cos(er_)*cos(ey_) + sin(ep_)*sin(er_)*sin(ey_), cos(er_)*sin(ep_)*sin(ey_) - cos(ey_)*sin(er_),
      -sin(ep_), cos(ey_)*sin(er_), cos(er_)*cos(ey_);
    // d er, ep_, ey
    R_local_dx_vec_[0] << 0, sin(er_)*sin(ey_) + cos(er_)*cos(ey_)*sin(ep_),   cos(er_)*sin(ey_) - cos(ey_)*sin(ep_)*sin(er_),
      0, cos(er_)*sin(ep_)*sin(ey_) - cos(ey_)*sin(er_), - cos(er_)*cos(ey_) - sin(ep_)*sin(er_)*sin(ey_),
      0, cos(er_)*cos(ey_), -cos(ey_)*sin(er_);

    R_local_dx_vec_[1] << -cos(ey_)*sin(ep_), cos(ep_)*cos(ey_)*sin(er_), cos(ep_)*cos(er_)*cos(ey_),
      -sin(ep_)*sin(ey_), cos(ep_)*sin(er_)*sin(ey_), cos(ep_)*cos(er_)*sin(ey_),
      -cos(ep_), 0, 0;

    R_local_dx_vec_[2] << -cos(ep_)*sin(ey_), - cos(er_)*cos(ey_) - sin(ep_)*sin(er_)*sin(ey_), cos(ey_)*sin(er_) - cos(er_)*sin(ep_)*sin(ey_),
      cos(ep_)*cos(ey_),   cos(ey_)*sin(ep_)*sin(er_) - cos(er_)*sin(ey_), sin(er_)*sin(ey_) + cos(er_)*cos(ey_)*sin(ep_),
      0, -sin(er_)*sin(ey_), -cos(er_)*sin(ey_);

    // dd er, ep, ey
    R_local_ddx_vec_[0] << 0,   cos(er_)*sin(ey_) - cos(ey_)*sin(ep_)*sin(er_), - sin(er_)*sin(ey_) - cos(er_)*cos(ey_)*sin(ep_),
      0, - cos(er_)*cos(ey_) - sin(ep_)*sin(er_)*sin(ey_),   cos(ey_)*sin(er_) - cos(er_)*sin(ep_)*sin(ey_),
      0,                            -cos(ey_)*sin(er_),                            -cos(er_)*cos(ey_);
    R_local_ddx_vec_[1] << 0, cos(ep_)*cos(er_)*cos(ey_), -cos(ep_)*cos(ey_)*sin(er_),
      0, cos(ep_)*cos(er_)*sin(ey_), -cos(ep_)*sin(er_)*sin(ey_),
      0,                       0,                        0;
    R_local_ddx_vec_[2] << 0, cos(ey_)*sin(er_) - cos(er_)*sin(ep_)*sin(ey_), cos(er_)*cos(ey_) + sin(ep_)*sin(er_)*sin(ey_),
      0, sin(er_)*sin(ey_) + cos(er_)*cos(ey_)*sin(ep_), cos(er_)*sin(ey_) - cos(ey_)*sin(ep_)*sin(er_),
      0,                          -cos(er_)*sin(ey_),                           sin(er_)*sin(ey_);
    R_local_ddx_vec_[3] << 0, cos(ep_)*cos(er_)*cos(ey_), -cos(ep_)*cos(ey_)*sin(er_),
      0, cos(ep_)*cos(er_)*sin(ey_), -cos(ep_)*sin(er_)*sin(ey_),
      0,                       0,                        0;
    R_local_ddx_vec_[4] << -cos(ep_)*cos(ey_), -cos(ey_)*sin(ep_)*sin(er_), -cos(er_)*cos(ey_)*sin(ep_),
      -cos(ep_)*sin(ey_), -sin(ep_)*sin(er_)*sin(ey_), -cos(er_)*sin(ep_)*sin(ey_),
      sin(ep_),                        0,                        0;
    R_local_ddx_vec_[5] << sin(ep_)*sin(ey_), -cos(ep_)*sin(er_)*sin(ey_), -cos(ep_)*cos(er_)*sin(ey_),
      -cos(ey_)*sin(ep_),  cos(ep_)*cos(ey_)*sin(er_),  cos(ep_)*cos(er_)*cos(ey_),
      0,                        0,                        0;
    R_local_ddx_vec_[6] << 0, cos(ey_)*sin(er_) - cos(er_)*sin(ep_)*sin(ey_), cos(er_)*cos(ey_) + sin(ep_)*sin(er_)*sin(ey_),
      0, sin(er_)*sin(ey_) + cos(er_)*cos(ey_)*sin(ep_), cos(er_)*sin(ey_) - cos(ey_)*sin(ep_)*sin(er_),
      0,                          -cos(er_)*sin(ey_),                           sin(er_)*sin(ey_);
    R_local_ddx_vec_[7] << sin(ep_)*sin(ey_), -cos(ep_)*sin(er_)*sin(ey_), -cos(ep_)*cos(er_)*sin(ey_),
      -cos(ey_)*sin(ep_),  cos(ep_)*cos(ey_)*sin(er_),  cos(ep_)*cos(er_)*cos(ey_),
      0,                        0,                        0;
    R_local_ddx_vec_[8] << -cos(ep_)*cos(ey_),   cos(er_)*sin(ey_) - cos(ey_)*sin(ep_)*sin(er_), - sin(er_)*sin(ey_) - cos(er_)*cos(ey_)*sin(ep_),
      -cos(ep_)*sin(ey_), - cos(er_)*cos(ey_) - sin(ep_)*sin(er_)*sin(ey_),   cos(ey_)*sin(er_) - cos(er_)*sin(ep_)*sin(ey_),
      0,                            -cos(ey_)*sin(er_),                            -cos(er_)*cos(ey_);
    // T_local_ddx
    T_local_ddx_vec_[0] << 0,        0,                0,
      0, -cos(er_), -cos(ep_)*sin(er_),
      0,  sin(er_), -cos(ep_)*cos(er_);
    T_local_ddx_vec_[1] << 0, 0,                0,
      0, 0, -cos(er_)*sin(ep_),
      0, 0,  sin(ep_)*sin(er_);
    T_local_ddx_vec_[3] << 0, 0,                0,
      0, 0, -cos(er_)*sin(ep_),
      0, 0,  sin(ep_)*sin(er_);
    T_local_ddx_vec_[4] << 0, 0,          sin(ep_),
      0, 0, -cos(ep_)*sin(er_),
      0, 0, -cos(ep_)*cos(er_);
    // T_local_ddx_vec: all 0 index: 2, 5, 6-8

    // R_link_local_vec_
    R_link_local_vec_[0] = MatrixXd::Identity(3, 3);
    for (int i = 1; i < 4; ++i){
      double cur_q = 0.0;
      for (int j = 1; j <= i; ++j) cur_q += q_vec_(j-1);
      R_link_local_vec_[i] << cos(cur_q), -sin(cur_q), 0,
        sin(cur_q),  cos(cur_q), 0,
        0,        0, 1;
    }
    // d R_link_local_vec_
      // d q1
    R_link_local_dx_vec_[3*1+0] << -sin(q1_), -cos(q1_), 0,
      cos(q1_), -sin(q1_), 0,
      0,        0, 0;
    R_link_local_dx_vec_[3*2+0] << -sin(q1_ + q2_), -cos(q1_ + q2_), 0,
      cos(q1_ + q2_), -sin(q1_ + q2_), 0,
      0,             0, 0;
    R_link_local_dx_vec_[3*3+0] << -sin(q1_ + q2_ + q3_), -cos(q1_ + q2_ + q3_), 0,
      cos(q1_ + q2_ + q3_), -sin(q1_ + q2_ + q3_), 0,
      0,                  0, 0;
      // d q2_
    R_link_local_dx_vec_[3*2+1] = R_link_local_dx_vec_[3*2+0];
    R_link_local_dx_vec_[3*3+1] = R_link_local_dx_vec_[3*3+0];
      // d q3
    R_link_local_dx_vec_[3*3+2] = R_link_local_dx_vec_[3*3+0];
    // dd R_link_local_vec_
      // R_link_local_vec_[1] dd
    R_link_local_ddx_vec_[9*1+3*0+0] = -R_link_local_vec_[1];
      // R_link_local_vec_[2] dd
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        R_link_local_ddx_vec_[9*2+3*i+j] = -R_link_local_vec_[2];
      // R_link_local_vec_[3] dd
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        R_link_local_ddx_vec_[9*3+3*i+j] = -R_link_local_vec_[3];


    // T_local_ and its derivative
    T_local_ << 1, 0, -sin(ep_),
      0, cos(er_), cos(ep_)*sin(er_),
      0, -sin(er_), cos(ep_)*cos(er_);
    MatrixXd T_d = MatrixXd::Zero(3, 3); // d er, ep_, ey
    T_d << 0, 0, 0,
      0, -sin(er_),  cos(ep_)*cos(er_),
      0, -cos(er_), -cos(ep_)*sin(er_);
    T_local_dx_vec_[0] = T_d;
    T_d << 0, 0, -cos(ep_),
      0, 0, -sin(ep_)*sin(er_),
      0, 0, -cos(er_)*sin(ep_);
    T_local_dx_vec_[1] = T_d;
    T_d = MatrixXd::Zero(3, 3);
    T_local_dx_vec_[2] = T_d;

    // Q_local and its derivative, d derivative
    Q_local_ = R_local_.transpose() * T_local_;
    for (int i = 0; i < 3; ++i){
      Q_local_dx_vec_[i] = R_local_dx_vec_[i].transpose() * T_local_
        + R_local_.transpose() * T_local_dx_vec_[i];
      // Q_local_ddx
      for (int j = 0; j < 3; ++j){
        Q_local_ddx_vec_[3*i+j] = R_local_ddx_vec_[3*i+j].transpose() * T_local_
          + R_local_dx_vec_[i].transpose() * T_local_dx_vec_[j]
          + R_local_dx_vec_[j].transpose() * T_local_dx_vec_[i]
          + R_local_.transpose() * T_local_ddx_vec_[3*i+j];
      }
    }

    // link_center and its derivative
    link_center_pos_local_ << link_length_/2, link_length_ + (link_length_*cos(q1_))/2, link_length_ + (link_length_*cos(q1_ + q2_))/2 + link_length_*cos(q1_), link_length_ + link_length_*cos(q1_ + q2_) + link_length_*cos(q1_) + (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, (link_length_*sin(q1_))/2, (link_length_*sin(q1_ + q2_))/2 + link_length_*sin(q1_), link_length_*sin(q1_ + q2_) + link_length_*sin(q1_) + (link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    // d q1_, q2_, q3_
    link_center_pos_local_dx_vec_[0] << 0, -(link_length_*sin(q1_))/2, - (link_length_*sin(q1_ + q2_))/2 - link_length_*sin(q1_), - link_length_*sin(q1_ + q2_) - link_length_*sin(q1_) - (link_length_*sin(q1_ + q2_ + q3_))/2,
      0,  (link_length_*cos(q1_))/2,   (link_length_*cos(q1_ + q2_))/2 + link_length_*cos(q1_),   link_length_*cos(q1_ + q2_) + link_length_*cos(q1_) + (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_dx_vec_[1] << 0, 0, -(link_length_*sin(q1_ + q2_))/2, - link_length_*sin(q1_ + q2_) - (link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0,  (link_length_*cos(q1_ + q2_))/2,   link_length_*cos(q1_ + q2_) + (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_dx_vec_[2] << 0, 0, 0, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0,  (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    // dd q1, q2, q3
    link_center_pos_local_ddx_vec_[0] << 0, -(link_length_*cos(q1_))/2, - (link_length_*cos(q1_ + q2_))/2 - link_length_*cos(q1_), - link_length_*cos(q1_ + q2_) - link_length_*cos(q1_) - (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, -(link_length_*sin(q1_))/2, - (link_length_*sin(q1_ + q2_))/2 - link_length_*sin(q1_), - link_length_*sin(q1_ + q2_) - link_length_*sin(q1_) - (link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_ddx_vec_[1] << 0, 0, -(link_length_*cos(q1_ + q2_))/2, - link_length_*cos(q1_ + q2_) - (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, -(link_length_*sin(q1_ + q2_))/2, - link_length_*sin(q1_ + q2_) - (link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_ddx_vec_[2] << 0, 0, 0, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_ddx_vec_[3] << 0, 0, -(link_length_*cos(q1_ + q2_))/2, - link_length_*cos(q1_ + q2_) - (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, -(link_length_*sin(q1_ + q2_))/2, - link_length_*sin(q1_ + q2_) - (link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_ddx_vec_[4] << 0, 0, -(link_length_*cos(q1_ + q2_))/2, - link_length_*cos(q1_ + q2_) - (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, -(link_length_*sin(q1_ + q2_))/2, - link_length_*sin(q1_ + q2_) - (link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_ddx_vec_[5] << 0, 0, 0, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_ddx_vec_[6] << 0, 0, 0, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0,  (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_ddx_vec_[7] = link_center_pos_local_ddx_vec_[6];
    link_center_pos_local_ddx_vec_[8] = link_center_pos_local_ddx_vec_[6];

    // Jacobian vec
    Jacobian_P_vec_[1] << -(link_length_*sin(q1_))/2, 0, 0,
      (link_length_*cos(q1_))/2, 0, 0,
      0, 0, 0;
    Jacobian_P_vec_[2] << -(link_length_*(sin(q1_ + q2_) + 2*sin(q1_)))/2, -(link_length_*sin(q1_ + q2_))/2, 0,
      (link_length_*(cos(q1_ + q2_) + 2*cos(q1_)))/2,  (link_length_*cos(q1_ + q2_))/2, 0,
      0,                             0, 0;
    Jacobian_P_vec_[3] << -(link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_) + 2*sin(q1_)))/2, -(link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      (link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_) + 2*cos(q1_)))/2,  (link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2,  (link_length_*cos(q1_ + q2_ + q3_))/2,
      0,                                                     0,                                  0;

    Jacobian_W_vec_[1](2, 0) = 1.0;
    Jacobian_W_vec_[2] = Jacobian_W_vec_[1]; Jacobian_W_vec_[2](2, 1) = 1.0;
    Jacobian_W_vec_[3] = Jacobian_W_vec_[2]; Jacobian_W_vec_[3](2, 2) = 1.0;
    // d Jacobian_P
      // d q1
    Jacobian_P_dq_vec_[1*3+0] << -(link_length_*cos(q1_))/2, 0, 0,
      -(link_length_*sin(q1_))/2, 0, 0,
      0, 0, 0;
    Jacobian_P_dq_vec_[2*3+0] << -(link_length_*(cos(q1_ + q2_) + 2*cos(q1_)))/2, -(link_length_*cos(q1_ + q2_))/2, 0,
      -(link_length_*(sin(q1_ + q2_) + 2*sin(q1_)))/2, -(link_length_*sin(q1_ + q2_))/2, 0,
      0,                             0, 0;
    Jacobian_P_dq_vec_[3*3+0] << -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_) + 2*cos(q1_)))/2, -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      -(link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_) + 2*sin(q1_)))/2, -(link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0;
      // d q2_
    Jacobian_P_dq_vec_[2*3+1] << -(link_length_*cos(q1_ + q2_))/2, -(link_length_*cos(q1_ + q2_))/2, 0,
      -(link_length_*sin(q1_ + q2_))/2, -(link_length_*sin(q1_ + q2_))/2, 0,
      0,                             0, 0;
    Jacobian_P_dq_vec_[3*3+1] << -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      -(link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, -(link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0;
      // d q3_
    Jacobian_P_dq_vec_[3*3+2] << -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      -(link_length_*sin(q1_ + q2_ + q3_))/2, -(link_length_*sin(q1_ + q2_ + q3_))/2, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0;

    // Jacobian_P_ddq_vec_
    Jacobian_P_ddq_vec_[9] << (link_length_*sin(q1_))/2, 0, 0,
      -(link_length_*cos(q1_))/2, 0, 0,
      0, 0, 0;
    Jacobian_P_ddq_vec_[18] << (link_length_*(sin(q1_ + q2_) + 2*sin(q1_)))/2, (link_length_*sin(q1_ + q2_))/2, 0,
      -(link_length_*(cos(q1_ + q2_) + 2*cos(q1_)))/2, -(link_length_*cos(q1_ + q2_))/2, 0,
      0, 0, 0;
    Jacobian_P_ddq_vec_[19] << (link_length_*sin(q1_ + q2_))/2, (link_length_*sin(q1_ + q2_))/2, 0,
      -(link_length_*cos(q1_ + q2_))/2, -(link_length_*cos(q1_ + q2_))/2, 0,
      0, 0, 0;
    Jacobian_P_ddq_vec_[21] << (link_length_*sin(q1_ + q2_))/2, (link_length_*sin(q1_ + q2_))/2, 0,
      -(link_length_*cos(q1_ + q2_))/2, -(link_length_*cos(q1_ + q2_))/2, 0,
      0, 0, 0;
    Jacobian_P_ddq_vec_[22] << (link_length_*sin(q1_ + q2_))/2, (link_length_*sin(q1_ + q2_))/2, 0,
      -(link_length_*cos(q1_ + q2_))/2, -(link_length_*cos(q1_ + q2_))/2, 0,
      0, 0, 0;
    Jacobian_P_ddq_vec_[27] << (link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_) + 2*sin(q1_)))/2, (link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_) + 2*cos(q1_)))/2, -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;
    Jacobian_P_ddq_vec_[28] << (link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, (link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;
    Jacobian_P_ddq_vec_[29] << (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;
    Jacobian_P_ddq_vec_[30] << (link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, (link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;
    Jacobian_P_ddq_vec_[31] << (link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, (link_length_*(sin(q1_ + q2_ + q3_) + 2*sin(q1_ + q2_)))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*(cos(q1_ + q2_ + q3_) + 2*cos(q1_ + q2_)))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;
    Jacobian_P_ddq_vec_[32] << (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;
    Jacobian_P_ddq_vec_[33] << (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;
    Jacobian_P_ddq_vec_[34] << (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;
    Jacobian_P_ddq_vec_[35] << (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2, (link_length_*sin(q1_ + q2_ + q3_))/2,
      -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2, -(link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0;

    // S(R_b * P_bli_b)
    S_operation_result_ = MatrixXd::Zero(3, 4);
    S_operation_result_ = R_local_ * link_center_pos_local_;
    // d S(R_b * P_bli_b)
    for (int i = 0; i < 3; ++i){ // d er,ep,ey
      S_operation_dx_vec_[i] = R_local_dx_vec_[i] * link_center_pos_local_;
      for (int j = 0; j < 3; ++j) // d e(i) d e(j)
        S_operation_ddx_vec_[i*6+j] = R_local_ddx_vec_[i*3+j] * link_center_pos_local_;
      for (int j = 0; j < 3; ++j) // d e(i) d q(j)
        S_operation_ddx_vec_[i*6+j+3] = R_local_dx_vec_[i] * link_center_pos_local_dx_vec_[j];
    }
    for (int i = 3; i < 6; ++i){ // d q1_,q2_,q3_
      S_operation_dx_vec_[i] = R_local_ * link_center_pos_local_dx_vec_[i-3];
      for (int j = 0; j < 3; ++j) // d q(i) d e(j)
        S_operation_ddx_vec_[i*6+j] = R_local_dx_vec_[j] * link_center_pos_local_dx_vec_[i-3];
      for (int j = 0; j < 3; ++j) // d q(i) d q(j)
        S_operation_ddx_vec_[i*6+j+3] = R_local_ * link_center_pos_local_ddx_vec_[(i-3)*3+j];
    }
  }

  MatrixXd HydrusDynamics::vectorToSkewMatrix(VectorXd s){
    MatrixXd res = MatrixXd::Zero(3, 3);
    res << 0.0, -s(2), s(1),
      s(2), 0.0, -s(0),
      -s(1), s(0), 0.0;
    return res;
  }

  void HydrusDynamics::updateMainMatrix(){
    Ds_ = MatrixXd::Zero(6, 6);
    for (int i = 0; i < 3; ++i) // D11
      Ds_(i, i) = weight_sum_;
    D12_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n_links_; ++i)
      D12_ = D12_ - link_weight_vec_[i] * vectorToSkewMatrix(S_operation_result_.col(i));
    D12_ = D12_ * T_local_;
    D13_ = MatrixXd::Zero(3, n_links_-1);
    for (int i = 0; i < n_links_; ++i)
      D13_ = D13_ + link_weight_vec_[i] * Jacobian_P_vec_[i];
    D13_ = R_local_ * D13_;
    D22_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n_links_; ++i)
      D22_ = D22_ + link_weight_vec_[i] * T_local_.transpose() *
        vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
        vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_
        + Q_local_.transpose() * R_link_local_vec_[i] * Inertial_ *
        R_link_local_vec_[i].transpose() * Q_local_;
    D23_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n_links_; ++i)
      D23_ = D23_ + Q_local_.transpose() * R_link_local_vec_[i] *
        Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i]
        - link_weight_vec_[i] * T_local_.transpose() *
        vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_ * Jacobian_P_vec_[i];
    D33_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n_links_; ++i)
      D33_ = D33_ + link_weight_vec_[i] * Jacobian_P_vec_[i].transpose() * Jacobian_P_vec_[i]
        + Jacobian_W_vec_[i].transpose() * R_link_local_vec_[i] * Inertial_ *
        R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];

    Ds_.block<3, 3>(0, 3) = D12_;
    Ds_.block<3, 3>(3, 0) = D12_.transpose();
    Ds_.block<3, 3>(3, 3) = D22_;
    Ds3_.block<3, 3>(0, 0) = D13_;
    Ds3_.block<3, 3>(3, 0) = D23_;
    Ds_inv_ = Ds_.inverse();

    // D_d
    for (int i = 0; i <= Q_3; ++i){ // init
      D_dx_vec_[i] = MatrixXd::Zero(9, 9);
      for (int j = 0; j <= Q_3; ++j)
        D_ddx_vec_[i*9+j] = MatrixXd::Zero(9, 9);
    }
    // D12_d
    for (int i = 0; i < n_links_; ++i){
      // d T_local
      for (int j = E_R; j <= E_Y; ++j){
        MatrixXd D12_d = -link_weight_vec_[i] * vectorToSkewMatrix(S_operation_result_.col(i))
          * T_local_dx_vec_[j - E_R];
        D_dx_vec_[j].block<3, 3>(0, 3) += D12_d;
        D_dx_vec_[j].block<3, 3>(3, 0) += D12_d.transpose();
        for (int k = E_R; k <= E_Y; ++k){ // d T_local d T_local
          MatrixXd D12_dd = -link_weight_vec_[i] * vectorToSkewMatrix(S_operation_result_.col(i))
            * T_local_ddx_vec_[3*(j-E_R) + (k-E_R)];
          D_ddx_vec_[j*9+k].block<3, 3>(0, 3) += D12_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(3, 0) += D12_dd.transpose();
        }
        for (int k = E_R; k <= Q_3; ++k){ // d T_local d S(R_b * P_bli_b)
          MatrixXd D12_dd = -link_weight_vec_[i] *
            vectorToSkewMatrix(S_operation_dx_vec_[k-E_R].col(i)) * T_local_dx_vec_[j-E_R];
          D_ddx_vec_[j*9+k].block<3, 3>(0, 3) += D12_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(3, 0) += D12_dd.transpose();
        }
      }
      // d S(R_b * P_bli_b)
      for (int j = E_R; j <= Q_3; ++j){
        MatrixXd D12_d = -link_weight_vec_[i] *
          vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)) * T_local_;
        D_dx_vec_[j].block<3, 3>(0, 3) += D12_d;
        D_dx_vec_[j].block<3, 3>(3, 0) += D12_d.transpose();
        for (int k = E_R; k <= E_Y; ++k){ // d S(R_b * P_bli_b) d T_local
          MatrixXd D12_dd = -link_weight_vec_[i] *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)) * T_local_dx_vec_[k-E_R];
          D_ddx_vec_[j*9+k].block<3, 3>(0, 3) += D12_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(3, 0) += D12_dd.transpose();
        }
        for (int k = E_R; k <= Q_3; ++k){ // d S(R_b * P_bli_b) d S(R_b * P_bli_b)
          MatrixXd D12_dd = -link_weight_vec_[i] *
            vectorToSkewMatrix(S_operation_ddx_vec_[(j-E_R)*6+(k-E_R)].col(i)) * T_local_;
          D_ddx_vec_[j*9+k].block<3, 3>(0, 3) += D12_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(3, 0) += D12_dd.transpose();
        }
      }
    }
    // D13_d
    for (int i = 0; i < n_links_; ++i){
      // d R_b
      for (int j = E_R; j <= E_Y; ++j){
        MatrixXd D13_d = link_weight_vec_[i] * R_local_dx_vec_[j-E_R] * Jacobian_P_vec_[i];
        D_dx_vec_[j].block<3, 3>(0, 6) += D13_d;
        D_dx_vec_[j].block<3, 3>(6, 0) += D13_d.transpose();
        for (int k = E_R; k <= E_Y; ++k){ // d R_b d R_b
          MatrixXd D13_dd = link_weight_vec_[i] * R_local_ddx_vec_[3*(j-E_R)+(k-E_R)]
            * Jacobian_P_vec_[i];
          D_ddx_vec_[j*9+k].block<3, 3>(0, 6) += D13_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 0) += D13_dd.transpose();
        }
        for (int k = Q_1; k <= Q_3; ++k){ // d R_b d J_p
          MatrixXd D13_dd = link_weight_vec_[i] * R_local_dx_vec_[j-E_R]
            * Jacobian_P_dq_vec_[i*3+(k-Q_1)];
          D_ddx_vec_[j*9+k].block<3, 3>(0, 6) += D13_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 0) += D13_dd.transpose();
        }
      }
      // d J_p
      for (int j = Q_1; j <= Q_3; ++j){
        MatrixXd D13_d = link_weight_vec_[i] * R_local_ * Jacobian_P_dq_vec_[i*3+j-Q_1];
        D_dx_vec_[j].block<3, 3>(0, 6) = D_dx_vec_[j].block<3, 3>(0, 6) + D13_d;
        D_dx_vec_[j].block<3, 3>(6, 0) = D_dx_vec_[j].block<3, 3>(6, 0) + D13_d.transpose();
        for (int k = E_R; k <= E_Y; ++k){ // d J_p d R_b
          MatrixXd D13_dd = link_weight_vec_[i] * R_local_dx_vec_[k-E_R]
            * Jacobian_P_dq_vec_[i*3+j-Q_1];
          D_ddx_vec_[j*9+k].block<3, 3>(0, 6) += D13_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 0) += D13_dd.transpose();
        }
        for (int k = Q_1; k <= Q_3; ++k){ // d J_p d J_p
          MatrixXd D13_dd = link_weight_vec_[i] * R_local_
            * Jacobian_P_ddq_vec_[i*9+(j-Q_1)*3+(k-Q_1)];
          D_ddx_vec_[j*9+k].block<3, 3>(0, 6) += D13_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 0) += D13_dd.transpose();
        }
      }
    }
    // D22_d
    for (int i = 0; i < n_links_; ++i){
      for (int j = E_R; j <= E_Y; ++j){
        // left part: d T_local & d T_local'
        MatrixXd D22_d = link_weight_vec_[i] * T_local_dx_vec_[j-E_R].transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_;
        D22_d += link_weight_vec_[i] * T_local_.transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)) *
          T_local_dx_vec_[j-E_R];
        for (int k = E_R; k <= E_Y; ++k){ // d T_local d T_local
          MatrixXd D22_dd = link_weight_vec_[i] * T_local_ddx_vec_[3*(j-E_R)+(k-E_R)].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_;
          D22_dd += link_weight_vec_[i] * T_local_dx_vec_[j-E_R].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_dx_vec_[k-E_R];
          D22_dd += link_weight_vec_[i] * T_local_dx_vec_[k-E_R].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)) *
            T_local_dx_vec_[j-E_R];
          D22_dd += link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)) *
            T_local_ddx_vec_[3*(j-E_R)+(k-E_R)];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 3) += D22_dd;
        }
        for (int k = E_R; k <= Q_3; ++k){ // d T_local d S(R_b * P_bli_b)
          MatrixXd D22_dd = link_weight_vec_[i] * T_local_dx_vec_[j-E_R].transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[k-E_R].col(i)).transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_;
          D22_dd += link_weight_vec_[i] * T_local_dx_vec_[j-E_R].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[k-E_R].col(i)) * T_local_;
          D_ddx_vec_[j*9+k].block<3, 3>(3, 3) += D22_dd;
        }
        // right part: d Q_local and their transpose
        D22_d += Q_local_dx_vec_[j-E_R].transpose() * R_link_local_vec_[i] *
          Inertial_ * R_link_local_vec_[i].transpose() * Q_local_;
        D22_d += Q_local_.transpose() * R_link_local_vec_[i] * Inertial_ *
          R_link_local_vec_[i].transpose() * Q_local_dx_vec_[j-E_R];
        D_dx_vec_[j].block<3, 3>(3, 3) = D_dx_vec_[j].block<3, 3>(3, 3) + D22_d;
        for (int k = E_R; k <= E_Y; ++k){ // d Q_local d Q_local
          MatrixXd D22_dd = Q_local_ddx_vec_[3*(j-E_R)+(k-E_R)].transpose() * R_link_local_vec_[i] *
            Inertial_ * R_link_local_vec_[i].transpose() * Q_local_;
          D22_dd += Q_local_dx_vec_[j-E_R].transpose() * R_link_local_vec_[i] *
            Inertial_ * R_link_local_vec_[i].transpose() * Q_local_dx_vec_[k-E_R];
          D22_dd += Q_local_dx_vec_[k-E_R].transpose() * R_link_local_vec_[i] * Inertial_ *
            R_link_local_vec_[i].transpose() * Q_local_dx_vec_[j-E_R];
          D22_dd += Q_local_.transpose() * R_link_local_vec_[i] * Inertial_ *
            R_link_local_vec_[i].transpose() * Q_local_ddx_vec_[3*(j-E_R)+(k-E_R)];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 3) += D22_dd;
        }
        for (int k = E_R; k <= E_Y; ++k){ // d Q_local d R_link_local
          MatrixXd D22_dd = Q_local_dx_vec_[j-E_R].transpose() * R_link_local_dx_vec_[3*i+(k-E_R)]
            * Inertial_ * R_link_local_vec_[i].transpose() * Q_local_;
          D22_dd += Q_local_dx_vec_[j-E_R].transpose() * R_link_local_vec_[i] *
            Inertial_ * R_link_local_dx_vec_[3*i+(k-E_R)].transpose() * Q_local_;
          D22_dd += Q_local_.transpose() * R_link_local_dx_vec_[3*i+(k-E_R)] * Inertial_ *
            R_link_local_vec_[i].transpose() * Q_local_dx_vec_[j-E_R];
          D22_dd += Q_local_.transpose() * R_link_local_vec_[i] * Inertial_ *
            R_link_local_dx_vec_[3*i+(k-E_R)].transpose() * Q_local_dx_vec_[j-E_R];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 3) += D22_dd;
        }
      }
      for (int j = E_R; j <= Q_3; ++j){
        // left part: d S(R_b * P_bli_b) & d S(R_b * P_bli_b)'
        MatrixXd D22_d = link_weight_vec_[i] * T_local_.transpose() *
          vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)).transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_;
        D22_d += link_weight_vec_[i] * T_local_.transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
          vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)) * T_local_;
        D_dx_vec_[j].block<3, 3>(3, 3) = D_dx_vec_[j].block<3, 3>(3, 3) + D22_d;
        for (int k = E_R; k <= E_Y; ++k){ // d S(R_b * P_bli_b) d T_local_
          MatrixXd D22_dd = link_weight_vec_[i] * T_local_dx_vec_[k-E_R].transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)).transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_;
          D22_dd += link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)).transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_dx_vec_[k-E_R];
          D22_dd += link_weight_vec_[i] * T_local_dx_vec_[k-E_R].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)) * T_local_;
          D22_dd += link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)) * T_local_dx_vec_[k-E_R];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 3) += D22_dd;
        }
        for (int k = E_R; k <= Q_3; ++k){ // d S(R_b * P_bli_b) d S(R_b * P_bli_b)
          MatrixXd D22_dd = link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_ddx_vec_[6*(j-E_R)+(k-E_R)].col(i)).transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)) * T_local_;
          D22_dd += link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)).transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[k-E_R].col(i)) * T_local_;
          D22_dd += link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[k-E_R].col(i)).transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)) * T_local_;
          D22_dd += link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            vectorToSkewMatrix(S_operation_ddx_vec_[6*(j-E_R)+(k-E_R)].col(i)) * T_local_;
          D_ddx_vec_[j*9+k].block<3, 3>(3, 3) += D22_dd;
        }
      }
      for (int j = Q_1; j <= Q_3; ++j){
        // right part: R_link_local & R_link_local'
        MatrixXd D22_d = Q_local_.transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
          Inertial_ * R_link_local_vec_[i].transpose() * Q_local_;
        D22_d += Q_local_.transpose() * R_link_local_vec_[i] * Inertial_ *
          R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Q_local_;
        D_dx_vec_[j].block<3, 3>(3, 3) = D_dx_vec_[j].block<3, 3>(3, 3) + D22_d;
        for (int k = E_R; k <= E_Y; ++k){ // d R_link_local d Q_local
          MatrixXd D22_dd = Q_local_dx_vec_[k-E_R].transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
            Inertial_ * R_link_local_vec_[i].transpose() * Q_local_;
          D22_dd += Q_local_.transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
            Inertial_ * R_link_local_vec_[i].transpose() * Q_local_dx_vec_[k-E_R];
          D22_dd += Q_local_dx_vec_[k-E_R].transpose() * R_link_local_vec_[i] * Inertial_ *
            R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Q_local_;
          D22_dd += Q_local_.transpose() * R_link_local_vec_[i] * Inertial_ *
            R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Q_local_dx_vec_[k-E_R];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 3) += D22_dd;
        }
        for (int k = Q_1; k <= Q_3; ++k){ // d R_link_local d R_link_local
          MatrixXd D22_dd = Q_local_.transpose() * R_link_local_ddx_vec_[9*i+3*(j-Q_1)+(k-Q_1)] *
            Inertial_ * R_link_local_vec_[i].transpose() * Q_local_;
          D22_dd += Q_local_.transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
            Inertial_ * R_link_local_dx_vec_[3*i+k-Q_1].transpose() * Q_local_;
          D22_dd += Q_local_.transpose() * R_link_local_dx_vec_[3*i+k-Q_1] *
            Inertial_ * R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Q_local_;
          D22_dd += Q_local_.transpose() * R_link_local_vec_[i] * Inertial_ *
            R_link_local_ddx_vec_[9*i+3*(j-Q_1)+(k-Q_1)].transpose() * Q_local_;
          D_ddx_vec_[j*9+k].block<3, 3>(3, 3) += D22_dd;
        }
      }
    }
    // D23_d
    for (int i = 0; i < n_links_; ++i){
      for (int j = E_R; j <= E_Y; ++j){
        // left part: d Q_local
        MatrixXd D23_d = Q_local_dx_vec_[j-E_R].transpose() * R_link_local_vec_[i] *
          Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];
        // right part: d T_local R_local_
        D23_d = D23_d - link_weight_vec_[i] * T_local_dx_vec_[j-E_R].transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_
          * Jacobian_P_vec_[i];
        D23_d = D23_d - link_weight_vec_[i] * T_local_.transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_dx_vec_[j-E_R]
          * Jacobian_P_vec_[i];
        D_dx_vec_[j].block<3, 3>(3, 6) = D_dx_vec_[j].block<3, 3>(3, 6) + D23_d;
        D_dx_vec_[j].block<3, 3>(6, 3) = D_dx_vec_[j].block<3, 3>(6, 3) + D23_d.transpose();
        for (int k = E_R; k <= E_Y; ++k){ // d  d (Q, T, R)
          MatrixXd D23_dd = Q_local_ddx_vec_[3*(j-E_R)+(k-E_R)].transpose() * R_link_local_vec_[i] *
            Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];
          D23_dd -= link_weight_vec_[i] * T_local_ddx_vec_[3*(j-E_R)+(k-E_R)].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_
            * Jacobian_P_vec_[i];
          D23_dd -= link_weight_vec_[i] * T_local_dx_vec_[j-E_R].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_dx_vec_[k-E_R]
            * Jacobian_P_vec_[i];
          D23_dd -= link_weight_vec_[i] * T_local_dx_vec_[k-E_R].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_dx_vec_[j-E_R]
            * Jacobian_P_vec_[i];
          D23_dd -= link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() *
            R_local_ddx_vec_[3*(j-E_R)+(k-E_R)] * Jacobian_P_vec_[i];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
        for (int k = Q_1; k <= Q_3; ++k){ // d  d J_p
          MatrixXd D23_dd = -link_weight_vec_[i] * T_local_dx_vec_[j-E_R].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_
            * Jacobian_P_dq_vec_[i*3+k-Q_1];
          D23_dd -= link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_dx_vec_[j-E_R]
            * Jacobian_P_dq_vec_[i*3+k-Q_1];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
        for (int k = Q_1; k <= Q_3; ++k){ // d  d R_link_local
          MatrixXd D23_dd = Q_local_dx_vec_[j-E_R].transpose() * R_link_local_dx_vec_[i*3+(k-Q_1)]
            * Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];
          D23_dd += Q_local_dx_vec_[j-E_R].transpose() * R_link_local_vec_[i] *
            Inertial_ * R_link_local_dx_vec_[i*3+(k-Q_1)].transpose() * Jacobian_W_vec_[i];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
        for (int k = E_R; k <= Q_3; ++k){ // d d S(R_b * P_bli_b)
          MatrixXd D23_dd = -link_weight_vec_[i] * T_local_dx_vec_[j-E_R].transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[k-E_R].col(i)).transpose() * R_local_
            * Jacobian_P_vec_[i];
          D23_dd -= link_weight_vec_[i] * T_local_.transpose() *
          vectorToSkewMatrix(S_operation_dx_vec_[k-E_R].col(i)).transpose() * R_local_dx_vec_[j-E_R]
          * Jacobian_P_vec_[i];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
      }
      for (int j = Q_1; j <= Q_3; ++j){
        // left part: R_link_local & R_link_local'
        MatrixXd D23_d = Q_local_.transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
          Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];
        D23_d = D23_d + Q_local_.transpose() * R_link_local_vec_[i] *
          Inertial_ * R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Jacobian_W_vec_[i];
        // right part: Jacobian_P
        D23_d = D23_d - link_weight_vec_[i] * T_local_.transpose() *
          vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_ *
          Jacobian_P_dq_vec_[i*3+j-Q_1];
        D_dx_vec_[j].block<3, 3>(3, 6) = D_dx_vec_[j].block<3, 3>(3, 6) + D23_d;
        D_dx_vec_[j].block<3, 3>(6, 3) = D_dx_vec_[j].block<3, 3>(6, 3) + D23_d.transpose();
        for (int k = E_R; k <= E_Y; ++k){ // d R_link_local d (Q, R, T)
          MatrixXd D23_dd = Q_local_dx_vec_[k-E_R].transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
            Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];
          D23_dd += Q_local_dx_vec_[k-E_R].transpose() * R_link_local_vec_[i] *
            Inertial_ * R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Jacobian_W_vec_[i];
          D23_dd -= link_weight_vec_[i] * T_local_dx_vec_[k-E_R].transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_ *
            Jacobian_P_dq_vec_[i*3+j-Q_1];
          D23_dd -= link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_dx_vec_[k-E_R] *
            Jacobian_P_dq_vec_[i*3+j-Q_1];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
        for (int k = Q_1; k <= Q_3; ++k){ // d R_link_local d (R_link, Jacobian_P)
          MatrixXd D23_dd = Q_local_.transpose() * R_link_local_ddx_vec_[9*i+3*(j-Q_1)+(k-Q_1)] *
            Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];
          D23_dd += Q_local_.transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
            Inertial_ * R_link_local_dx_vec_[3*i+k-Q_1].transpose() * Jacobian_W_vec_[i];
          D23_dd += Q_local_.transpose() * R_link_local_dx_vec_[3*i+k-Q_1] *
            Inertial_ * R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Jacobian_W_vec_[i];
          D23_dd += D23_d + Q_local_.transpose() * R_link_local_vec_[i] *
            Inertial_ * R_link_local_ddx_vec_[9*i+3*(j-Q_1)+(k-Q_1)].transpose() * Jacobian_W_vec_[i];
          D23_d -= link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_result_.col(i)).transpose() * R_local_ *
            Jacobian_P_ddq_vec_[i*9+(j-Q_1)*3+(k-Q_1)];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
        for (int k = E_R; k <= Q_3; ++k){ // d R_link_local d S(R_b * P_bli_b)
          MatrixXd D23_dd = - link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[k-E_R].col(i)).transpose() * R_local_ *
            Jacobian_P_dq_vec_[i*3+j-Q_1];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
      }
      for (int j = E_R; j <= Q_3; ++j){
        // right part: d S(R_b * P_bli_b)
        MatrixXd D23_d = - link_weight_vec_[i] * T_local_.transpose() *
          vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)).transpose()
          * R_local_ * Jacobian_P_vec_[i];
        D_dx_vec_[j].block<3, 3>(3, 6) = D_dx_vec_[j].block<3, 3>(3, 6) + D23_d;
        D_dx_vec_[j].block<3, 3>(6, 3) = D_dx_vec_[j].block<3, 3>(6, 3) + D23_d.transpose();
        for (int k = E_R; k <= E_Y; ++k){ // d S(R_b * P_bli_b) d (T, R)
          MatrixXd D23_dd = -link_weight_vec_[i] * T_local_dx_vec_[k-E_R].transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)).transpose()
            * R_local_ * Jacobian_P_vec_[i];
          D23_dd -= link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)).transpose()
            * R_local_dx_vec_[k-E_R] * Jacobian_P_vec_[i];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
        for (int k = Q_1; k <= Q_3; ++k){ // d S(R_b * P_bli_b) d Jacobian_P_vec
          MatrixXd D23_dd = -link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_dx_vec_[j-E_R].col(i)).transpose()
            * R_local_ * Jacobian_P_dq_vec_[i*3+(k-Q_1)];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
        for (int k = E_R; k <= Q_3; ++k){ // d S(R_b * P_bli_b) d S(R_b * P_bli_b)
          MatrixXd D23_dd = -link_weight_vec_[i] * T_local_.transpose() *
            vectorToSkewMatrix(S_operation_ddx_vec_[6*(j-E_R)+(k-E_R)].col(i)).transpose()
            * R_local_ * Jacobian_P_vec_[i];
          D_ddx_vec_[j*9+k].block<3, 3>(3, 6) += D23_dd;
          D_ddx_vec_[j*9+k].block<3, 3>(6, 3) += D23_dd.transpose();
        }
      }
    }
    // D33_d
    for (int i = 0; i < n_links_; ++i){
      for (int j = E_R; j <= E_Y; ++j){
        // right part: d R_link_local
        MatrixXd D33_d = Jacobian_W_vec_[i].transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
          Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];
        D33_d = D33_d + Jacobian_W_vec_[i].transpose() * R_link_local_vec_[i] * Inertial_ *
          R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Jacobian_W_vec_[i];
        D_dx_vec_[j].block<3, 3>(6, 6) = D_dx_vec_[j].block<3, 3>(6, 6) + D33_d;
        for (int k = Q_1; k <= Q_3; ++k){ // d R_link_local d d R_link_local
          MatrixXd D33_dd = Jacobian_W_vec_[i].transpose() * R_link_local_ddx_vec_[9*i+3*(j-Q_1)+(k-Q_1)] *
            Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i];
          D33_dd += Jacobian_W_vec_[i].transpose() * R_link_local_dx_vec_[3*i+j-Q_1] *
            Inertial_ * R_link_local_dx_vec_[3*i+k-Q_1] * Jacobian_W_vec_[i];
          D33_dd += Jacobian_W_vec_[i].transpose() * R_link_local_dx_vec_[3*i+k-Q_1] * Inertial_ *
            R_link_local_dx_vec_[3*i+j-Q_1].transpose() * Jacobian_W_vec_[i];
          D33_dd += Jacobian_W_vec_[i].transpose() * R_link_local_vec_[i] * Inertial_ *
            R_link_local_ddx_vec_[9*i+3*(j-Q_1)+(k-Q_1)].transpose() * Jacobian_W_vec_[i];
          D_ddx_vec_[j*9+k].block<3, 3>(6, 6) += D33_dd;
        }
      }
      for (int j = Q_1; j <= Q_3; ++j){
        // left part: d Jacobian_P
        MatrixXd D33_d = link_weight_vec_[i] * Jacobian_P_dq_vec_[i*3+j-Q_1].transpose() *
          Jacobian_P_vec_[i];
        D33_d = D33_d + link_weight_vec_[i] * Jacobian_P_vec_[i].transpose() *
          Jacobian_P_dq_vec_[i*3+j-Q_1];
        D_dx_vec_[j].block<3, 3>(6, 6) = D_dx_vec_[j].block<3, 3>(6, 6) + D33_d;
        for (int k = Q_1; k <= Q_3; ++k){ // d Jacobian_P d Jacobian_P
          MatrixXd D33_dd = link_weight_vec_[i] * Jacobian_P_ddq_vec_[i*9+3*(j-Q_1)+(k-Q_1)].transpose() *
            Jacobian_P_vec_[i];
          D33_dd += link_weight_vec_[i] * Jacobian_P_dq_vec_[i*3+j-Q_1].transpose() *
            Jacobian_P_dq_vec_[i*3+k-Q_1];
          D33_dd += link_weight_vec_[i] * Jacobian_P_dq_vec_[i*3+k-Q_1].transpose() *
            Jacobian_P_dq_vec_[i*3+j-Q_1];
          D33_dd += link_weight_vec_[i] * Jacobian_P_vec_[i].transpose() *
            Jacobian_P_ddq_vec_[i*9+3*(j-Q_1)+(k-Q_1)];
          D_ddx_vec_[j*9+k].block<3, 3>(6, 6) += D33_dd;
        }
      }
    }

    // C
    C_ = MatrixXd::Zero(9, 9);
    for (int i = 0; i < 9; ++i){
      C_dx_vec_[i] = MatrixXd::Zero(9, 9);
      C_d_dx_vec_[i] = MatrixXd::Zero(9, 9);
    }
    for (int k = 0; k <= Q_3; ++k)
      for (int j = 0; j <= Q_3; ++j)
        for (int i = 0; i <= Q_3; ++i){
          double x_d;
          if (i < 6) x_d = x_vec_[6+i]; // d s
          else x_d = q_vec_[3+(i-6)]; // d q
          double C_param = 0.5 * (D_dx_vec_[i](k, j) + D_dx_vec_[j](k, i) + D_dx_vec_[k](i, j));
          C_(k, j) += C_param * x_d;
          C_d_dx_vec_[i](k, j) += C_param;
          for (int i2 = 0; i2 <= Q_3; ++i2){
            double C_dx_param = 0.5 * (D_ddx_vec_[i*9+i2](k, j) + D_ddx_vec_[j*9+i2](k, i)
                                        + D_dx_vec_[k*9+i2](i, j));
            C_dx_vec_[i2](k, j) += C_dx_param * x_d;
          }
        }
    // gs
    gs_ = VectorXd::Zero(6);
    for (int i = 0; i < 3; ++i)
      gs_dx_vec_[i] = VectorXd::Zero(6);
    for (int i = 0; i < n_links_; ++i){
      Vector3d mid_result = link_weight_vec_[i] * 9.78 * Vector3d(0, 0, 1.0);
      // P_X P_Y P_Z
      gs_(P_X) = gs_(P_X) + mid_result.dot(Vector3d(1.0, 0, 0));
      gs_(P_Y) = gs_(P_Y) + mid_result.dot(Vector3d(0, 1.0, 0));
      gs_(P_Z) = gs_(P_Z) + mid_result.dot(Vector3d(0, 0, 1.0));
      // E_R E_P E_Y
      for (int j = E_R; j <= E_Y; ++j){
        gs_(j) = gs_(j) + mid_result.dot(R_local_dx_vec_[j-E_R]
                                         * link_center_pos_local_.col(i));
        for (int k = E_R; k <= E_Y; ++k)
          gs_dx_vec_[j](k) += mid_result.dot(R_local_ddx_vec_[3*(j-E_R)+k-E_R]
                                             * link_center_pos_local_.col(i));

      }
      // Q_1 Q_2 Q_3. do not have for simplified state
      //   for (int j = Q_1; j <= Q_3; ++j)
      //     gs_(j) = gs_(j) + mid_result.dot(R_local_
      //                                      * link_center_pos_local_dx_vec_[j-Q_1].col(i));
    }

    // Bs
    Bs_ = VectorXd::Zero(6);
    double f_sum = 0;
    for (int i = 0; i < 4; ++i)
      f_sum += u_vec_(i);
    Bs_.head(3) = R_local_ * Vector3d(0, 0, f_sum); // force
    VectorXd Bs_tau = VectorXd::Zero(3); // torque
    for (int i = 0; i < n_links_; ++i){
      Vector3d p_lci_b(link_center_pos_local_.col(i)(0), link_center_pos_local_.col(i)(1),
                       link_center_pos_local_.col(i)(2));
      Bs_.tail(3) = Bs_.tail(3) + p_lci_b.cross(Vector3d(0, 0, u_vec_(i)));
    }
    Bs_.tail(3) = Q_local_.transpose() * Bs_.tail(3);
    // Bs_du
    VectorXd Bs_du = R_local_ * Vector3d(0, 0, 1.0);
    for (int i = 0; i < 4; ++i)
      Bs_du_vec_[i].head(3) = Bs_du;
    for (int i = 0; i < n_links_; ++i){
      Vector3d p_lci_b(link_center_pos_local_.col(i)(0), link_center_pos_local_.col(i)(1),
                       link_center_pos_local_.col(i)(2));
      Bs_du_vec_[i].tail(3) = Q_local_.transpose() * p_lci_b.cross(Vector3d(0, 0, 1));
    }
    // Bs_dx
    for (int i = 0; i < 3; ++i)
      Bs_dx_vec_[i] = VectorXd::Zero(6);
    for (int i = E_R; i <= E_Y; ++i)
      Bs_dx_vec_[i-E_R].head(3) = R_local_dx_vec_[i-E_R] * Vector3d(0, 0, f_sum);
    for (int i = 0; i < n_links_; ++i){
      Vector3d p_lci_b(link_center_pos_local_.col(i)(0), link_center_pos_local_.col(i)(1),
                       link_center_pos_local_.col(i)(2));
      for (int j = E_R; j <= E_Y; ++j)
        Bs_dx_vec_[j-E_R].tail(3) = Bs_dx_vec_[j-E_R].tail(3) +
          Q_local_dx_vec_[j-E_R].transpose() * p_lci_b.cross(Vector3d(0, 0, u_vec_(i)));
    }
  }
}


