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
    Cs_ = MatrixXd::Zero(6, 6);
    Cs3_ = MatrixXd::Zero(6, 3);

    for (int i = 0; i < 3; ++i){
      R_local_d_vec_.push_back(MatrixXd::Zero(3, 3));
      T_local_d_vec_.push_back(MatrixXd::Zero(3, 3));
      link_center_pos_local_d_vec_.push_back(MatrixXd::Zero(3, 4));
    }
    for (int i = 0; i < 4; ++i){
      R_link_local_vec_.push_back(MatrixXd::Zero(3, 3));
      Jacobian_P_vec_.push_back(MatrixXd::Zero(3, 3));
      Jacobian_W_vec_.push_back(MatrixXd::Zero(3, 3));
    }
  }

  HydrusDynamics::~HydrusDynamics(){
  }

  void HydrusDynamics::getCurrentState(VectorXd *s_ptr, VectorXd *q_ptr){
    x_vec_ = *s_ptr;
    q_vec_ = *q_ptr;
    px_ = (*s_ptr)[0];
    py_ = (*s_ptr)[1];
    pz_ = (*s_ptr)[2];
    er_ = (*s_ptr)[3];
    ep_ = (*s_ptr)[4];
    ey_ = (*s_ptr)[5];
    d_px_ = (*s_ptr)[6];
    d_py_ = (*s_ptr)[7];
    d_pz_ = (*s_ptr)[8];
    d_er_ = (*s_ptr)[9];
    d_ep_ = (*s_ptr)[10];
    d_ey_ = (*s_ptr)[11];

    q1_ = (*q_ptr)[0];
    q2_ = (*q_ptr)[1];
    q3_ = (*q_ptr)[2];

    // mid result
    // R_local_ and its derivative
    R_local_ << cos(ep_)*cos(ey_), cos(ey_)*sin(ep_)*sin(er_) - cos(er_)*sin(ey_), sin(er_)*sin(ey_) + cos(er_)*cos(ey_)*sin(ep_),
      cos(ep_)*sin(ey_), cos(er_)*cos(ey_) + sin(ep_)*sin(er_)*sin(ey_), cos(er_)*sin(ep_)*sin(ey_) - cos(ey_)*sin(er_),
      -sin(ep_), cos(ey_)*sin(er_), cos(er_)*cos(ey_);
    MatrixXd rot_d = MatrixXd::Zero(3, 3); // d er, ep_, ey
    rot_d << 0, sin(er_)*sin(ey_) + cos(er_)*cos(ey_)*sin(ep_),   cos(er_)*sin(ey_) - cos(ey_)*sin(ep_)*sin(er_),
      0, cos(er_)*sin(ep_)*sin(ey_) - cos(ey_)*sin(er_), - cos(er_)*cos(ey_) - sin(ep_)*sin(er_)*sin(ey_),
      0, cos(er_)*cos(ey_), -cos(ey_)*sin(er_);
    R_local_d_vec_[0] = rot_d;

    rot_d << -cos(ey_)*sin(ep_), cos(ep_)*cos(ey_)*sin(er_), cos(ep_)*cos(er_)*cos(ey_),
      -sin(ep_)*sin(ey_), cos(ep_)*sin(er_)*sin(ey_), cos(ep_)*cos(er_)*sin(ey_),
      -cos(ep_), 0, 0;
    R_local_d_vec_[1] = rot_d;

    rot_d << -cos(ep_)*sin(ey_), - cos(er_)*cos(ey_) - sin(ep_)*sin(er_)*sin(ey_), cos(ey_)*sin(er_) - cos(er_)*sin(ep_)*sin(ey_),
      cos(ep_)*cos(ey_),   cos(ey_)*sin(ep_)*sin(er_) - cos(er_)*sin(ey_), sin(er_)*sin(ey_) + cos(er_)*cos(ey_)*sin(ep_),
      0, -sin(er_)*sin(ey_), -cos(er_)*sin(ey_);
    R_local_d_vec_[2] = rot_d;

    // R_link_local_vec_
    R_link_local_vec_[0] = MatrixXd::Identity(3, 3);
    for (int i = 1; i < 4; ++i){
      double cur_q = 0.0;
      for (int j = 1; j <= i; ++j) cur_q += q_vec_(j-1);
      R_link_local_vec_[i] << cos(cur_q), -sin(cur_q), 0,
        sin(cur_q),  cos(cur_q), 0,
        0,        0, 1;
    }

    // T_local_ and its derivative
    T_local_ << 1, 0, -sin(ep_),
      0, cos(er_), cos(ep_)*sin(er_),
      0, -sin(er_), cos(ep_)*cos(er_);
    MatrixXd T_d = MatrixXd::Zero(3, 3); // d er, ep_, ey
    T_d << 0, 0, 0,
      0, -sin(er_),  cos(ep_)*cos(er_),
      0, -cos(er_), -cos(ep_)*sin(er_);
    T_local_d_vec_[0] = T_d;
    T_d << 0, 0, -cos(ep_),
      0, 0, -sin(ep_)*sin(er_),
      0, 0, -cos(er_)*sin(ep_);
    T_local_d_vec_[1] = T_d;
    T_d = MatrixXd::Zero(3, 3);
    T_local_d_vec_[2] = T_d;

    // link_center and its derivative
    link_center_pos_local_ << link_length_/2, link_length_ + (link_length_*cos(q1_))/2, link_length_ + (link_length_*cos(q1_ + q2_))/2 + link_length_*cos(q1_), link_length_ + link_length_*cos(q1_ + q2_) + link_length_*cos(q1_) + (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, (link_length_*sin(q1_))/2, (link_length_*sin(q1_ + q2_))/2 + link_length_*sin(q1_), link_length_*sin(q1_ + q2_) + link_length_*sin(q1_) + (link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    MatrixXd link_pos_d = MatrixXd::Zero(3, 3); // d q1_, q2_, q3_
    link_pos_d << 0, -(link_length_*sin(q1_))/2, - (link_length_*sin(q1_ + q2_))/2 - link_length_*sin(q1_), - link_length_*sin(q1_ + q2_) - link_length_*sin(q1_) - (link_length_*sin(q1_ + q2_ + q3_))/2,
      0,  (link_length_*cos(q1_))/2,   (link_length_*cos(q1_ + q2_))/2 + link_length_*cos(q1_),   link_length_*cos(q1_ + q2_) + link_length_*cos(q1_) + (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_d_vec_[0] = link_pos_d;
    link_pos_d << 0, 0, -(link_length_*sin(q1_ + q2_))/2, - link_length_*sin(q1_ + q2_) - (link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0,  (link_length_*cos(q1_ + q2_))/2,   link_length_*cos(q1_ + q2_) + (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_d_vec_[1] = link_pos_d;
    link_pos_d << 0, 0, 0, -(link_length_*sin(q1_ + q2_ + q3_))/2,
      0, 0, 0,  (link_length_*cos(q1_ + q2_ + q3_))/2,
      0, 0, 0, 0;
    link_center_pos_local_d_vec_[2] = link_pos_d;

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

    // S(R_b * P_bli_b)
    S_operation_result_ = MatrixXd::Zero(3, 4);
    S_operation_result_ = R_local_ * link_center_pos_local_;
  }

  MatrixXd HydrusDynamics::vectorSkewToMatrix(Vector3d s){
    MatrixXd res = MatrixXd::Zero(3, 3);
    res << 0.0, -s(3), s(2),
      s(3), 0.0, -s(1),
      -s(2), s(1), 0.0;
  }

  void HydrusDynamics::updateMatrixD(){
    Ds_ = MatrixXd::Zero(6, 6);
    for (int i = 0; i < 3; ++i) // D11
      Ds_(i, i) = weight_sum_;
    D12_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n_links_; ++i)
      D12_ = D12_ - link_weight_vec_[i] * S_operation_result_.col(i);
    D12_ = D12_ * T_local_;
    D13_ = MatrixXd::Zero(3, n_links_-1);
    for (int i = 0; i < n_links_; ++i)
      D13_ = D13_ + link_weight_vec_[i] * Jacobian_P_vec_[i];
    D13_ = R_local_ * D13_;
    D22_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n_links_; ++i)
      D22_ = D22_ + link_weight_vec_[i] * T_local_.transpose() *
        S_operation_result_.col(i).transpose() * S_operation_result_.col(i) * T_local_
        + R_local_ * R_link_local_vec_[i] * Inertial_ *
        R_link_local_vec_[i].transpose() * R_local_.transpose();
    D23_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n_links_; ++i)
      D23_ = D23_ + T_local_.transpose() * R_local_ * R_link_local_vec_[i] *
        Inertial_ * R_link_local_vec_[i].transpose() * Jacobian_W_vec_[i]
        - link_weight_vec_[i] * T_local_.transpose() * S_operation_result_.col(i).transpose() *
        R_local_ * Jacobian_W_vec_[i];
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
  }

  void HydrusDynamics::calculateSkewOperation(std::vector<MatrixXd> *S_vec_ptr){
    for (int i = 0; i < 3; ++i) // d er,ep,ey
      S_vec_ptr->push_back(R_local_d_vec_[i] * link_center_pos_local_);
    for (int i = 0; i < 3; ++i) // d q1_,q2_,q3_
      S_vec_ptr->push_back(R_local_ * link_center_pos_local_d_vec_[i]);
  }
}


