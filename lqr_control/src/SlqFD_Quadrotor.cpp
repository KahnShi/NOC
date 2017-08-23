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

#include <lqr_control/SlqFD_Quadrotor.h>
namespace lqr_discrete{
  void SlqFiniteDiscreteControlQuadrotor::initSLQ(double freq, double period, VectorXd *x0, VectorXd *xn){
    /* Ros service */
    dare_client_ = nh_.serviceClient<lqr_control::Dare>("dare_solver");

    control_freq_ = freq;
    if (floor(freq * period) < freq * period){
      end_time_ = (floor(freq * period) + 1.0) / freq;
      iteration_times_ = floor(freq * period) + 1;
    }
    else{
      end_time_ = period;
      iteration_times_ = floor(freq * period);
    }
    std::cout << "[SLQ] Trajectory period: " << end_time_
              << ", Itetation times: " << iteration_times_ << "\n";
    quaternion_mode_ = false;
    if (quaternion_mode_)
      x_size_ = 13;
    else
      x_size_ = 12;
    u_size_ = 4;
    A_ptr_ = new MatrixXd(x_size_, x_size_);
    B_ptr_ = new MatrixXd(x_size_, u_size_);
    Q_ptr_ = new MatrixXd(x_size_, x_size_);
    R_ptr_ = new MatrixXd(u_size_, u_size_);
    x0_ptr_ = new VectorXd(x_size_);
    xn_ptr_ = new VectorXd(x_size_);
    x_ptr_ = new VectorXd(x_size_);
    u_ptr_ = new VectorXd(u_size_);
    M_para_ptr_ = new MatrixXd(3, 4);
    Riccati_P_ptr_ = new MatrixXd(x_size_, x_size_);
    P_ptr_ = new MatrixXd(x_size_, x_size_);
    p_ptr_ = new VectorXd(x_size_);
    H_ptr_ = new MatrixXd(u_size_, u_size_);
    G_ptr_ = new MatrixXd(u_size_, x_size_);
    K_ptr_ = new MatrixXd(u_size_, x_size_);
    g_ptr_ = new VectorXd(u_size_);
    l_ptr_ = new VectorXd(u_size_);
    r_ptr_ = new VectorXd(u_size_);
    q_ptr_ = new VectorXd(x_size_);

    *x0_ptr_ = (*x0);
    *xn_ptr_ = (*xn);

    /* init Q and R matrice */
    *Q_ptr_ = MatrixXd::Zero(x_size_, x_size_);
    for (int i = 0; i <= P_Z; ++i)
      (*Q_ptr_)(i, i) = 10.0;
    for (int i = V_X; i <= V_Z; ++i)
      (*Q_ptr_)(i, i) = 10.0;
    for (int i = W_X; i < x_size_; ++i)
      (*Q_ptr_)(i, i) = 1.0;

    *R_ptr_ = 50*MatrixXd::Identity(u_size_, u_size_);

    /* uav property from paper eth15-slq-window */
    I_ptr_ = new MatrixXd(3, 3);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        (*I_ptr_)(i, j) = 0.0;
    (*I_ptr_)(0, 0) = 0.03;
    (*I_ptr_)(1, 1) = 0.03;
    (*I_ptr_)(2, 2) = 0.05;

    uav_mass_ = 1.0;

    *M_para_ptr_ = MatrixXd::Zero(3, 4);
    double r_uav = 0.3, c_rf = 0.016;
    (*M_para_ptr_)(0, 1) = -r_uav;
    (*M_para_ptr_)(0, 3) = r_uav;
    (*M_para_ptr_)(1, 0) = r_uav;
    (*M_para_ptr_)(1, 2) = -r_uav;
    (*M_para_ptr_)(2, 0) = (*M_para_ptr_)(2, 2) = -c_rf;
    (*M_para_ptr_)(2, 1) = (*M_para_ptr_)(2, 3) = c_rf;

    uav_rotor_thrust_min_ = 0.0;
    uav_rotor_thrust_max_ = (uav_mass_ * 9.78 / 4.0) * 2.0;

    /* Assume initial and final state is still, namely dx = [v, a] = 0 */
    u0_ptr_ = new VectorXd(u_size_);
    un_ptr_ = new VectorXd(u_size_);
    for (int i = 0; i < 4; ++i){
      (*u0_ptr_)(i) = uav_mass_ * 9.78 / 4.0;
      (*un_ptr_)(i) = uav_mass_ * 9.78 / 4.0;
    }

    /* SLQ special initialization */
    // todo: assume start point the quadrotor is hovering
    VectorXd x_init(x_size_), u_init(u_size_);
    x_init = getRelativeState(x0_ptr_);
    u_init = VectorXd::Zero(u_size_);
    // test: real u
    // u_init = (*un_ptr_);

    for (int i = 0; i <= iteration_times_; ++i){
      x_vec_.push_back(x_init);
      u_vec_.push_back(u_init);
      u_fw_vec_.push_back(u_init);
      u_fb_vec_.push_back(u_init);
      K_vec_.push_back(MatrixXd::Zero(u_size_, x_size_));
    }

    debug_ = true;
    FDLQR();
    getRiccatiH();
    ROS_INFO("[SLQ] Initialization finished.");
  }

  void SlqFiniteDiscreteControlQuadrotor::getRiccatiH(){
    *x_ptr_ = x_vec_[iteration_times_];
    *u_ptr_ = u_vec_[iteration_times_];
    // *x_ptr_ = VectorXd::Zero(x_size_);
    // *u_ptr_ = VectorXd::Zero(u_size_);
    if (debug_){
      printStateInfo(x_ptr_, iteration_times_);
      printControlInfo(u_ptr_, iteration_times_);
    }
    if (quaternion_mode_)
      updateMatrixAB(x_ptr_, u_ptr_);
    else
      updateEulerMatrixAB(x_ptr_, u_ptr_);
    if (debug_){
      printMatrixAB();
    }

    lqr_control::float64Array dat_A, dat_B, dat_P;

    // fill out message:
    dat_A.array.layout.dim.push_back(std_msgs::MultiArrayDimension());
    dat_A.array.layout.dim.push_back(std_msgs::MultiArrayDimension());
    dat_A.array.layout.dim[0].label = "height";
    dat_A.array.layout.dim[1].label = "width";
    dat_A.array.layout.dim[0].size = x_size_;
    dat_A.array.layout.dim[1].size = x_size_;
    dat_A.array.layout.dim[0].stride = x_size_ * x_size_;
    dat_A.array.layout.dim[1].stride = x_size_;
    dat_A.array.layout.data_offset = 0;
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < x_size_; ++j)
        dat_A.array.data.push_back((*A_ptr_)(i, j));
    dat_B.array.layout.dim.push_back(std_msgs::MultiArrayDimension());
    dat_B.array.layout.dim.push_back(std_msgs::MultiArrayDimension());
    dat_B.array.layout.dim[0].label = "height";
    dat_B.array.layout.dim[1].label = "width";
    dat_B.array.layout.dim[0].size = x_size_;
    dat_B.array.layout.dim[1].size = u_size_;
    dat_B.array.layout.dim[0].stride = x_size_ * u_size_;
    dat_B.array.layout.dim[1].stride = u_size_;
    dat_B.array.layout.data_offset = 0;
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < u_size_; ++j)
        dat_B.array.data.push_back((*B_ptr_)(i, j));

    lqr_control::Dare dare_srv;
    dare_srv.request.A = dat_A;
    dare_srv.request.B = dat_B;
    ROS_INFO("[SLQ] Matrix AB is sent to Riccati solver.");
    if (dare_client_.call(dare_srv))
      dat_P = dare_srv.response.P;
    else
      ROS_ERROR("[SLQ] No response from dare sever.");
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < x_size_; ++j)
        (*Riccati_P_ptr_)(i, j) = dat_P.array.data[i * x_size_ + j];
    ROS_INFO("[SLQ] Matrix P is received from Riccati solver.");

    if (debug_){
      std::cout << "[Debug] print matrix P initial value:\n";
      for (int i = 0; i < x_size_; ++i){
        for (int j = 0; j < x_size_; ++j)
          std::cout << (*Riccati_P_ptr_)(i, j) << ", ";
        std::cout << "\n";
      }
    }

    // test P vector
    // if (debug_){
      // MatrixXd F = ((*R_ptr_) + B_ptr_->transpose() * (*Riccati_P_ptr_) * (*B_ptr_)).inverse() * (B_ptr_->transpose() * (*Riccati_P_ptr_) * (*A_ptr_));
      // VectorXd u = -F * (*x_ptr_);
      // std::cout << "[Debug] Test P's performance:\n";
      // VectorXd new_x(x_size_);
      // updateEulerNewState(&new_x, x_ptr_, &u);
      // VectorXd new_absolute_x;
      // new_absolute_x = getAbsoluteState(&new_x);
      // std::cout << "\n\n[debug] print current state:\n";
      // for (int j = 0; j < x_size_; ++j)
      //   std::cout << new_absolute_x(j) << ", ";
      // std::cout << "\n[debug] print current u:\n";
      // for (int j = 0; j < u_size_; ++j)
      //   std::cout << u(j) + (*un_ptr_)(j) << ", ";
      // std::cout << "\n";
      // }

    ROS_INFO("[SLQ] Get P matrice initial value from Ricatti function.");
  }

  void SlqFiniteDiscreteControlQuadrotor::updateAll(){
  }

  void SlqFiniteDiscreteControlQuadrotor::iterativeOptimization(){
    // *p_ptr_ = VectorXd::Zero(x_size_);
    // *r_ptr_ = (*R_ptr_) * (*un_ptr_);
    *P_ptr_ = *Riccati_P_ptr_;
    *p_ptr_ = (*P_ptr_) * (x_vec_[iteration_times_]);
    // test: real u
    // *r_ptr_ = VectorXd::Zero(u_size_);

    for (int i = iteration_times_ - 1; i >= 0; --i){
      // test: add weight for waypoint
      //updateQWeight(i * end_time_ / iteration_times_);

      *x_ptr_ = x_vec_[i];
      *u_ptr_ = u_vec_[i];
      if (quaternion_mode_)
        updateMatrixAB(x_ptr_, u_ptr_);
      else
        updateEulerMatrixAB(x_ptr_, u_ptr_);

      *q_ptr_ = (*Q_ptr_) * (x_vec_[i]);
      *r_ptr_ = (*R_ptr_) * (u_vec_[i]);
      // *r_ptr_ = (*R_ptr_) * ((u_vec_[i]) + (*un_ptr_));
      updateSLQEquations();

      Vector4d u_fb = (*K_ptr_) * (*x_ptr_);
      u_fb_vec_[i] = u_fb;
      u_fw_vec_[i] = (*l_ptr_);
      K_vec_[i] = (*K_ptr_);

      if ((i % 50 == 0 || i == iteration_times_ - 1) && debug_){
        printStateInfo(x_ptr_, i);
        printControlInfo(u_ptr_, i);
        std::cout << "\n[Debug] id[" << i << "] u feedback: ";
        for (int j = 0; j < u_size_; ++j)
          std::cout << u_fb_vec_[i](j) << ", ";
        std::cout << "\n[Debug] id[" << i << "] u feedforward: ";
        for (int j = 0; j < u_size_; ++j)
          std::cout << u_fw_vec_[i](j) << ", ";
        std::cout << "\n\n";
      }
    }

    /* Update control by finding the best alpha */
    alpha_ = 1.0;
    double alpha_candidate = 1.0, energy_min = -1.0;
    // if (feedforwardConverged() && debug_)
    //   std::cout << "[SLQ] feedforward converge.";
    // else{
    //   while (1){
    //     bool u_dynamic_flag = true;
    //     for (int i = 0; i < iteration_times_; ++i){
    //       VectorXd new_u = u_vec_[i] + (*un_ptr_)
    //         + alpha_ * u_fw_vec_[i] + u_fb_vec_[i];
    //       // test: real u
    //       //VectorXd new_u = u_vec_[i]
    //       // + alpha_ * u_fw_vec_[i] + u_fb_vec_[i];

    //       for (int j = 0; j < u_size_; ++j){
    //         if (new_u(j) < uav_rotor_thrust_min_
    //             || new_u(j) > uav_rotor_thrust_max_){
    //           u_dynamic_flag = false;
    //           break;
    //         }
    //       }
    //       if (!u_dynamic_flag){
    //         std::cout << "alpha: " << alpha_ << ", out of dynamic limitation\n";
    //         break;
    //       }
    //     }
    //     // maximum line search steps reached
    //     if (alpha_ < 1.0/32.0){
    //       alpha_ = 0.0;
    //       break;
    //     }
    //     else if (u_dynamic_flag){
    //       double energy = getSystemEnergy();
    //       if (energy_min < 0 || energy < energy_min){
    //         energy_min = energy;
    //         alpha_candidate = alpha_;
    //       }
    //     }
    //     alpha_ /= 2.0;
    //   }
    // }

    if (debug_)
      std::cout << "\nAlpha selected: " << alpha_candidate << "\n\n";

    // for (int i = 0; i < iteration_times_; ++i){
    //   VectorXd u = u_vec_[i] + alpha_candidate * u_fw_vec_[i]
    //     + u_fb_vec_[i];
    //   // Guarantee control is in bound
    //   checkControlInputFeasible(&u);
    //   u_vec_[i] = u;
    // }

    // *u_ptr_ = u_vec_[0];
    // *x_ptr_ = x_vec_[0];
    // for (int i = 0; i < iteration_times_; ++i){
    //   if (quaternion_mode_)
    //     updateMatrixAB(x_ptr_, u_ptr_);
    //   else
    //     updateEulerMatrixAB(x_ptr_, u_ptr_);
    //   VectorXd new_x(x_size_);
    //   if (quaternion_mode_)
    //     updateNewState(&new_x, x_ptr_, u_ptr_);
    //   else
    //     updateEulerNewState(&new_x, x_ptr_, u_ptr_);

    //   if ((i % 10 == 0 || i == iteration_times_ - 1) && debug_){
    //     printStateInfo(&new_x, i);
    //     printControlInfo(u_ptr_, i);
    //   }

    //   *u_ptr_ = u_vec_[i + 1];
    //   x_vec_[i + 1] = new_x;
    //   *x_ptr_ = new_x;
    // }


    // test K with every new state
    alpha_candidate = 0.5; // test
    VectorXd cur_x(x_size_);
    cur_x = x_vec_[0];
    for (int i = 0; i < iteration_times_; ++i){
      VectorXd cur_u(u_size_);
      cur_u = u_vec_[i] + alpha_candidate * u_fw_vec_[i]
        + K_vec_[i] * (cur_x - x_vec_[i]);
      checkControlInputFeasible(&cur_u);
      VectorXd new_x(x_size_);
      updateEulerNewState(&new_x, &cur_x, &cur_u);
      if ((i % 10 == 0 || i == iteration_times_ - 1) && debug_){
        printStateInfo(&cur_x, i);
        printControlInfo(&cur_u, i);
      }
      x_vec_[i] = cur_x;
      u_vec_[i] = cur_u;
      cur_x = new_x;
      if (i == iteration_times_ - 1)
        x_vec_[iteration_times_] = new_x;
    }
  }

  void SlqFiniteDiscreteControlQuadrotor::updateMatrixAB(VectorXd *x_ptr, VectorXd *u_ptr){
    updateMatrixA(x_ptr, u_ptr);
    updateMatrixB(x_ptr, u_ptr);
  }

  void SlqFiniteDiscreteControlQuadrotor::updateMatrixA(VectorXd *x_ptr, VectorXd *u_ptr){
    *A_ptr_ = MatrixXd::Zero(x_size_, x_size_);

    /* x, y, z */
    (*A_ptr_)(P_X, V_X) = 1;
    (*A_ptr_)(P_Y, V_Y) = 1;
    (*A_ptr_)(P_Z, V_Z) = 1;

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += ((*u_ptr)[i] + (*un_ptr_)[i]);

    /* u' = u / m */
    u = u / uav_mass_;
    /* d v_x = (2 * q_w * q_y + 2 * q_x * q_z) * u' */
    (*A_ptr_)(V_X, Q_W) = 2 * (*x_ptr)[Q_Y] * u;
    (*A_ptr_)(V_X, Q_X) = 2 * (*x_ptr)[Q_Z] * u;
    (*A_ptr_)(V_X, Q_Y) = 2 * (*x_ptr)[Q_W] * u;
    (*A_ptr_)(V_X, Q_Z) = 2 * (*x_ptr)[Q_X] * u;
    /* d v_y = (2 * q_y * q_z - 2 * q_w * q_x) * u' */
    (*A_ptr_)(V_Y, Q_W) = -2 * (*x_ptr)[Q_X] * u;
    (*A_ptr_)(V_Y, Q_X) = -2 * (*x_ptr)[Q_W] * u;
    (*A_ptr_)(V_Y, Q_Y) = 2 * (*x_ptr)[Q_Z] * u;
    (*A_ptr_)(V_Y, Q_Z) = 2 * (*x_ptr)[Q_Y] * u;
    /* d v_z = (q_w ^2 - q_x ^2 - q_y ^2 + q_z ^2) * u' */
    (*A_ptr_)(V_Z, Q_W) = 2 * (*x_ptr)[Q_W] * u;
    (*A_ptr_)(V_Z, Q_X) = -2 * (*x_ptr)[Q_X] * u;
    (*A_ptr_)(V_Z, Q_Y) = -2 * (*x_ptr)[Q_Y] * u;
    (*A_ptr_)(V_Z, Q_Z) = 2 * (*x_ptr)[Q_Z] * u;

    /* q_w, q_x, q_y, q_z */
    /* d q = 1/2 * q * [0, w_w]^T (the multiply opeation is under quaternion multiply) */
    /* d q = 1/2 * q * [0 R] * [0, w_b]^T */
    /* d q w = 1/2 * q * [0 R] * [0, d w_b]^T */
    Quaterniond qw((*x_ptr)[Q_W], (*x_ptr)[Q_X], (*x_ptr)[Q_Y], (*x_ptr)[Q_Z]);
    MatrixXd rot = qw.normalized().toRotationMatrix();
    MatrixXd rot4 = MatrixXd::Zero(4, 4);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        rot4(i+1, j+1) = rot(i, j);
    MatrixXd q_m = MatrixXd::Zero(4, 4);
    q_m << qw.w(), -qw.x(), -qw.y(), -qw.z(),
      qw.x(), qw.w(), qw.z(), -qw.y(),
      qw.y(), -qw.z(), qw.w(), qw.x(),
      qw.z(), qw.y(), -qw.x(), qw.w();
    Vector4d d_q_w_x = 0.5 * q_m * rot4 * Vector4d(0, 1, 0, 0);
    Vector4d d_q_w_y = 0.5 * q_m * rot4 * Vector4d(0, 0, 1, 0);
    Vector4d d_q_w_z = 0.5 * q_m * rot4 * Vector4d(0, 0, 0, 1);

    /* d q q = 1/2 * (d q * [0 R] + q * [0 dR]) * [0, w_b]^T */
    Vector4d w_b4(0, (*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    // d qw
    MatrixXd q_m_w = MatrixXd::Zero(4, 4);
    q_m_w(0, 0) = q_m_w(1, 1) = q_m_w(2, 2) = q_m_w(3, 3) = 1;
    MatrixXd rot4_w = MatrixXd::Zero(4, 4);
    rot4_w << 0, 0, 0, 0,
      0, qw.w(), -qw.z(), qw.y(),
      0, qw.z(), qw.w(), -qw.x(),
      0, -qw.y(), qw.x(), qw.w();
    rot4_w = 2 * rot4_w;
    Vector4d d_q_q_w = 0.5 * (q_m_w * rot4 + q_m * rot4_w) * w_b4;
    // d qx
    MatrixXd q_m_x = MatrixXd::Zero(4, 4);
    q_m_x(1, 0) = q_m_x(2, 3) = 1; q_m_x(0, 1) = q_m_x(3, 2) = -1;
    MatrixXd rot4_x = MatrixXd::Zero(4, 4);
    rot4_x << 0, 0, 0, 0,
      0, qw.x(), qw.y(), qw.z(),
      0, qw.y(), -qw.x(), -qw.w(),
      0, qw.z(), qw.w(), -qw.x();
    rot4_x = 2 * rot4_x;
    Vector4d d_q_q_x = 0.5 * (q_m_x * rot4 + q_m * rot4_x) * w_b4;
    // d qy
    MatrixXd q_m_y = MatrixXd::Zero(4, 4);
    q_m_y(2, 0) = q_m_y(3, 1) = 1; q_m_y(0, 2) = q_m_y(1, 3) = -1;
    MatrixXd rot4_y = MatrixXd::Zero(4, 4);
    rot4_y << 0, 0, 0, 0,
      0, -qw.y(), qw.x(), qw.w(),
      0, qw.x(), qw.y(), qw.z(),
      0, -qw.w(), qw.z(), qw.y();
    rot4_y = 2 * rot4_y;
    Vector4d d_q_q_y = 0.5 * (q_m_y * rot4 + q_m * rot4_y) * w_b4;
    // d qz
    MatrixXd q_m_z = MatrixXd::Zero(4, 4);
    q_m_z(1, 2) = q_m_z(3, 0) = 1; q_m_z(2, 1) = q_m_z(0, 3) = -1;
    MatrixXd rot4_z = MatrixXd::Zero(4, 4);
    rot4_z << 0, 0, 0, 0,
      0, -qw.z(), -qw.w(), qw.x(),
      0, qw.w(), -qw.z(), qw.y(),
      0, qw.x(), qw.y(), qw.z();
    rot4_z = 2 * rot4_z;
    Vector4d d_q_q_z = 0.5 * (q_m_z * rot4 + q_m * rot4_z) * w_b4;

    (*A_ptr_)(Q_W, Q_W) = d_q_q_w(0);
    (*A_ptr_)(Q_W, Q_X) = d_q_q_x(0);
    (*A_ptr_)(Q_W, Q_Y) = d_q_q_y(0);
    (*A_ptr_)(Q_W, Q_Z) = d_q_q_z(0);
    (*A_ptr_)(Q_W, W_X) = d_q_w_x(0);
    (*A_ptr_)(Q_W, W_Y) = d_q_w_y(0);
    (*A_ptr_)(Q_W, W_Z) = d_q_w_z(0);

    (*A_ptr_)(Q_X, Q_W) = d_q_q_w(1);
    (*A_ptr_)(Q_X, Q_X) = d_q_q_x(1);
    (*A_ptr_)(Q_X, Q_Y) = d_q_q_y(1);
    (*A_ptr_)(Q_X, Q_Z) = d_q_q_z(1);
    (*A_ptr_)(Q_X, W_X) = d_q_w_x(1);
    (*A_ptr_)(Q_X, W_Y) = d_q_w_y(1);
    (*A_ptr_)(Q_X, W_Z) = d_q_w_z(1);

    (*A_ptr_)(Q_Y, Q_W) = d_q_q_w(2);
    (*A_ptr_)(Q_Y, Q_X) = d_q_q_x(2);
    (*A_ptr_)(Q_Y, Q_Y) = d_q_q_y(2);
    (*A_ptr_)(Q_Y, Q_Z) = d_q_q_z(2);
    (*A_ptr_)(Q_Y, W_X) = d_q_w_x(2);
    (*A_ptr_)(Q_Y, W_Y) = d_q_w_y(2);
    (*A_ptr_)(Q_Y, W_Z) = d_q_w_z(2);

    (*A_ptr_)(Q_Z, Q_W) = d_q_q_w(3);
    (*A_ptr_)(Q_Z, Q_X) = d_q_q_x(3);
    (*A_ptr_)(Q_Z, Q_Y) = d_q_q_y(3);
    (*A_ptr_)(Q_Z, Q_Z) = d_q_q_z(3);
    (*A_ptr_)(Q_Z, W_X) = d_q_w_x(3);
    (*A_ptr_)(Q_Z, W_Y) = d_q_w_y(3);
    (*A_ptr_)(Q_Z, W_Z) = d_q_w_z(3);

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + tau), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    /* d w_w = I^-1 * (- d(w^) * (Iw) - (w^) * (I * d(w))) */
    Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd w_m = MatrixXd::Zero(3, 3);
    w_m << 0, -(*x_ptr)[W_Z], (*x_ptr)[W_Y],
      (*x_ptr)[W_Z], 0, -(*x_ptr)[W_X],
      -(*x_ptr)[W_Y], (*x_ptr)[W_X], 0;
    MatrixXd dw_m = MatrixXd::Zero(3, 3);
    dw_m(1, 2) = -1;
    dw_m(2, 1) = 1;
    Vector3d dw_x = I_ptr_->inverse() *
      ((-dw_m * ((*I_ptr_) * w))
       - w_m * ((*I_ptr_) * Vector3d(1.0, 0.0, 0.0)));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_X) = dw_x(i);

    dw_m = MatrixXd::Zero(3, 3);
    dw_m(0, 2) = 1;
    dw_m(2, 0) = -1;
    Vector3d dw_y = I_ptr_->inverse() *
      ((-dw_m * ((*I_ptr_) * w))
       - w_m * ((*I_ptr_) * Vector3d(0.0, 1.0, 0.0)));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_Y) = dw_y(i);

    dw_m = MatrixXd::Zero(3, 3);
    dw_m(0, 1) = -1;
    dw_m(1, 0) = 1;
    Vector3d dw_z = I_ptr_->inverse() *
      ((-dw_m * ((*I_ptr_) * w))
       - w_m * ((*I_ptr_) * Vector3d(0.0, 0.0, 1.0)));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_Z) = dw_z(i);

    (*A_ptr_) = (*A_ptr_) / control_freq_ + MatrixXd::Identity(x_size_, x_size_);
  }

  void SlqFiniteDiscreteControlQuadrotor::updateMatrixB(VectorXd *x_ptr, VectorXd *u_ptr){
    *B_ptr_ = MatrixXd::Zero(x_size_, u_size_);

    /* x, y, z */
    /* all 0 */

    /* v_x, v_y, v_z */
    /* d v_x = (2 * q_w * q_y + 2 * q_x * q_z) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_X, U_1) = (2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_Y] + 2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_Z]) / uav_mass_;
    /* d v_y = (2 * q_y * q_z - 2 * q_w * q_x) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_Y, U_1) = (2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Z] - 2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_X]) / uav_mass_;
    /* d v_z = (1 - 2 * q_x * q_x - 2 * q_y * q_y) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_Z, U_1) = (1 -2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_X] - 2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Y]) / uav_mass_;
    for (int i = V_X; i <= V_Z; ++i)
      for (int j = U_2; j <= U_4; ++j)
      (*B_ptr_)(i, j) = (*B_ptr_)(i, U_1);

    /* q_w, q_x, q_y, q_z */
    /* d q = 1/2 * q * [0, w_b]^T */
    /* all 0 */

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + M_para * [u1;u2;u3;u4]) */
    /* d w_u = I^-1 * M_para * d[u1;u2;u3;u4] */
    Vector3d dw_u1 = I_ptr_->inverse() * (*M_para_ptr_) * Vector4d(1.0, 0.0, 0.0, 0.0);
    for (int i = 0; i < 3; ++i)
      (*B_ptr_)(W_X + i, U_1) = dw_u1(i);

    Vector3d dw_u2 = I_ptr_->inverse() * (*M_para_ptr_) * Vector4d(0.0, 1.0, 0.0, 0.0);
    for (int i = 0; i < 3; ++i)
      (*B_ptr_)(W_X + i, U_2) = dw_u2(i);

    Vector3d dw_u3 = I_ptr_->inverse() * (*M_para_ptr_) * Vector4d(0.0, 0.0, 1.0, 0.0);
    for (int i = 0; i < 3; ++i)
      (*B_ptr_)(W_X + i, U_3) = dw_u3(i);

    Vector3d dw_u4 = I_ptr_->inverse() * (*M_para_ptr_) * Vector4d(0.0, 0.0, 0.0, 1.0);
    for (int i = 0; i < 3; ++i)
      (*B_ptr_)(W_X + i, U_4) = dw_u4(i);

    (*B_ptr_) = (*B_ptr_) / control_freq_;
  }

  void SlqFiniteDiscreteControlQuadrotor::updateNewState(VectorXd *new_x_ptr, VectorXd *x_ptr, VectorXd *u_ptr){
    VectorXd dev_x = VectorXd::Zero(x_size_);
    /* x, y, z */
    dev_x(P_X) = (*x_ptr)(V_X);
    dev_x(P_Y) = (*x_ptr)(V_Y);
    dev_x(P_Z) = (*x_ptr)(V_Z);

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += ((*u_ptr)[i] + (*un_ptr_)[i]);
    Quaterniond qw((*x_ptr)[Q_W], (*x_ptr)[Q_X], (*x_ptr)[Q_Y], (*x_ptr)[Q_Z]);
    MatrixXd rot = qw.normalized().toRotationMatrix();
    Vector3d acc = rot * Vector3d(0, 0, u / uav_mass_) - Vector3d(0, 0, 9.78);
    dev_x(V_X) = acc(0);
    dev_x(V_Y) = acc(1);
    dev_x(V_Z) = acc(2);

    /* q_w, q_x, q_y, q_z */
    /* d q = 1/2 * q * [0, R * w_b]^T (the multiply opeation is under quaternion multiply) */
    MatrixXd rot4 = MatrixXd::Zero(4, 4);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        rot4(i+1, j+1) = rot(i, j);
    MatrixXd q_m = MatrixXd::Zero(4, 4);
    q_m << qw.w(), -qw.x(), -qw.y(), -qw.z(),
      qw.x(), qw.w(), qw.z(), -qw.y(),
      qw.y(), -qw.z(), qw.w(), qw.x(),
      qw.z(), qw.y(), -qw.x(), qw.w();
    Vector4d w_b4(0, (*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    Vector4d dev_q = 0.5 * q_m * rot4 * w_b4;
    dev_x(Q_W) = dev_q(0);
    dev_x(Q_X) = dev_q(1);
    dev_x(Q_Y) = dev_q(2);
    dev_x(Q_Z) = dev_q(3);

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + M_para * [u1;u2;u3;u4]), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd w_m = MatrixXd::Zero(3, 3);
    w_m << 0, -(*x_ptr)[W_Z], (*x_ptr)[W_Y],
      (*x_ptr)[W_Z], 0, -(*x_ptr)[W_X],
      -(*x_ptr)[W_Y], (*x_ptr)[W_X], 0;
    Vector3d dw;
    dw = I_ptr_->inverse() * (-w_m * ((*I_ptr_) * w) + (*M_para_ptr_) * ((*u_ptr) + (*un_ptr_)));
    for (int i = 0; i < 3; ++i)
      dev_x(W_X + i) = dw(i);

    dev_x = dev_x / control_freq_;
    *new_x_ptr = dev_x + *x_ptr;
    // test
    //*new_x_ptr = stateAddition(x_ptr, &dev_x);

    normalizeQuaternion(new_x_ptr);
  }

  void SlqFiniteDiscreteControlQuadrotor::updateEulerMatrixAB(VectorXd *x_ptr, VectorXd *u_ptr){
    updateEulerMatrixA(x_ptr, u_ptr);
    updateEulerMatrixB(x_ptr, u_ptr);
  }

  void SlqFiniteDiscreteControlQuadrotor::updateEulerMatrixA(VectorXd *x_ptr, VectorXd *u_ptr){
    *A_ptr_ = MatrixXd::Zero(x_size_, x_size_);

    /* x, y, z */
    (*A_ptr_)(P_X, V_X) = 1;
    (*A_ptr_)(P_Y, V_Y) = 1;
    (*A_ptr_)(P_Z, V_Z) = 1;

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += ((*u_ptr)[i] + (*un_ptr_)[i]);
      // test: real u
      // u += (*u_ptr)[i];

    /* u' = u / m */
    u = u / uav_mass_;
    /* d v_x = (sin y * sin r + cos y * sin p * cos r) * u' */
    (*A_ptr_)(V_X, E_R) = (sin((*x_ptr)[E_Y]) * cos((*x_ptr)[E_R]) -
                           cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * sin((*x_ptr)[E_R])) * u;
    (*A_ptr_)(V_X, E_P) = cos((*x_ptr)[E_Y]) * cos((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]) * u;
    (*A_ptr_)(V_X, E_Y) = (cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) -
                           sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R])) * u;
    /* d v_y = (-cos y * sin r + sin y * sin p * cos r) * u' */
    (*A_ptr_)(V_Y, E_R) = (-cos((*x_ptr)[E_Y]) * cos((*x_ptr)[E_R]) -
                           sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]))* u;
    (*A_ptr_)(V_Y, E_P) = sin((*x_ptr)[E_Y]) * cos((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]) * u;
    (*A_ptr_)(V_Y, E_Y) = (sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) +
                           cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))* u;
    /* d v_z = (cos p * cos r) * u' */
    (*A_ptr_)(V_Z, E_P) = -sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]) * u;
    (*A_ptr_)(V_Z, E_R) = -cos((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]) * u;

    /* e_r, e_p, e_y */
    /* d e = R_e * w_b */
    Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd R_e = MatrixXd::Zero(3, 3);
    R_e << 1, tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]), tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]), -sin((*x_ptr)[E_R]),
      0, sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    Vector3d d_e_w_x = R_e * Vector3d(1.0, 0, 0);
    Vector3d d_e_w_y = R_e * Vector3d(0, 1.0, 0);
    Vector3d d_e_w_z = R_e * Vector3d(0, 0, 1.0);
    MatrixXd R_e_r = MatrixXd::Zero(3, 3);
    R_e_r << 0, tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]), -tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]),
      0, -sin((*x_ptr)[E_R]), -cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), -sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    Vector3d d_e_e_r = R_e_r * w;
    MatrixXd R_e_p = MatrixXd::Zero(3, 3);
    double d_cosp = sin((*x_ptr)[E_P]) / pow(cos((*x_ptr)[E_P]), 2.0);
    R_e_p << 0, sin((*x_ptr)[E_R]) / pow(cos((*x_ptr)[E_P]), 2.0), cos((*x_ptr)[E_R]) / pow(cos((*x_ptr)[E_P]), 2.0),
      0, 0, 0,
      0, sin((*x_ptr)[E_R]) * d_cosp, cos((*x_ptr)[E_R]) * d_cosp;
    Vector3d d_e_e_p = R_e_p * w;
    for (int i = E_R; i <= E_Y; ++i){
      (*A_ptr_)(i, W_X) = d_e_w_x(i - E_R);
      (*A_ptr_)(i, W_Y) = d_e_w_y(i - E_R);
      (*A_ptr_)(i, W_Z) = d_e_w_z(i - E_R);
      (*A_ptr_)(i, E_R) = d_e_e_r(i - E_R);
      (*A_ptr_)(i, E_P) = d_e_e_p(i - E_R);
    }

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + tau), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    /* d w_w = I^-1 * (- d(w^) * (Iw) - (w^) * (I * d(w))) */
    MatrixXd w_m = MatrixXd::Zero(3, 3);
    w_m << 0, -(*x_ptr)[W_Z], (*x_ptr)[W_Y],
      (*x_ptr)[W_Z], 0, -(*x_ptr)[W_X],
      -(*x_ptr)[W_Y], (*x_ptr)[W_X], 0;
    MatrixXd dw_m = MatrixXd::Zero(3, 3);
    dw_m(1, 2) = -1;
    dw_m(2, 1) = 1;
    Vector3d dw_x = I_ptr_->inverse() *
      ((-dw_m * ((*I_ptr_) * w))
       - w_m * ((*I_ptr_) * Vector3d(1.0, 0.0, 0.0)));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_X) = dw_x(i);

    dw_m = MatrixXd::Zero(3, 3);
    dw_m(0, 2) = 1;
    dw_m(2, 0) = -1;
    Vector3d dw_y = I_ptr_->inverse() *
      ((-dw_m * ((*I_ptr_) * w))
       - w_m * ((*I_ptr_) * Vector3d(0.0, 1.0, 0.0)));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_Y) = dw_y(i);

    dw_m = MatrixXd::Zero(3, 3);
    dw_m(0, 1) = -1;
    dw_m(1, 0) = 1;
    Vector3d dw_z = I_ptr_->inverse() *
      ((-dw_m * ((*I_ptr_) * w))
       - w_m * ((*I_ptr_) * Vector3d(0.0, 0.0, 1.0)));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_Z) = dw_z(i);

    (*A_ptr_) = (*A_ptr_) / control_freq_ + MatrixXd::Identity(x_size_, x_size_);
  }

  void SlqFiniteDiscreteControlQuadrotor::updateEulerMatrixB(VectorXd *x_ptr, VectorXd *u_ptr){
    *B_ptr_ = MatrixXd::Zero(x_size_, u_size_);

    /* x, y, z */
    /* all 0 */

    /* v_x, v_y, v_z */
    /* d v_x = (sin y * sin r + cos y * sin p * cos r) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_X, U_1) = (sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) +
                           cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      / uav_mass_;
    /* d v_y = (-cos y * sin r + sin y * sin p * cos r) * (u1 + u2 + u3 + u4) / m  */
    (*B_ptr_)(V_Y, U_1) = (-cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) +
                           sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      / uav_mass_;
    /* d v_z = (cos p * cos r) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_Z, U_1) = (cos((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      / uav_mass_;
    for (int i = V_X; i <= V_Z; ++i)
      for (int j = U_2; j <= U_4; ++j)
      (*B_ptr_)(i, j) = (*B_ptr_)(i, U_1);

    /* e_r, e_p, e_y */
    /* all 0 */

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + M_para * [u1;u2;u3;u4]) */
    /* d w_u = I^-1 * M_para * d[u1;u2;u3;u4] */
    Vector3d dw_u1 = I_ptr_->inverse() * (*M_para_ptr_) * Vector4d(1.0, 0.0, 0.0, 0.0);
    for (int i = 0; i < 3; ++i)
      (*B_ptr_)(W_X + i, U_1) = dw_u1(i);

    Vector3d dw_u2 = I_ptr_->inverse() * (*M_para_ptr_) * Vector4d(0.0, 1.0, 0.0, 0.0);
    for (int i = 0; i < 3; ++i)
      (*B_ptr_)(W_X + i, U_2) = dw_u2(i);

    Vector3d dw_u3 = I_ptr_->inverse() * (*M_para_ptr_) * Vector4d(0.0, 0.0, 1.0, 0.0);
    for (int i = 0; i < 3; ++i)
      (*B_ptr_)(W_X + i, U_3) = dw_u3(i);

    Vector3d dw_u4 = I_ptr_->inverse() * (*M_para_ptr_) * Vector4d(0.0, 0.0, 0.0, 1.0);
    for (int i = 0; i < 3; ++i)
      (*B_ptr_)(W_X + i, U_4) = dw_u4(i);

    (*B_ptr_) = (*B_ptr_) / control_freq_;
  }

  void SlqFiniteDiscreteControlQuadrotor::updateEulerNewState(VectorXd *new_x_ptr, VectorXd *x_ptr, VectorXd *u_ptr){
    VectorXd dev_x = VectorXd::Zero(x_size_);
    /* x, y, z */
    dev_x(P_X) = (*x_ptr)(V_X);
    dev_x(P_Y) = (*x_ptr)(V_Y);
    dev_x(P_Z) = (*x_ptr)(V_Z);

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
       u += ((*u_ptr)[i] + (*un_ptr_)[i]);
      // test: real u
      // u += (*u_ptr)[i];
    /* d v_x = (sin y * sin r + cos y * sin p * cos r) * (u1 + u2 + u3 + u4) / m */
    dev_x(V_X) = (sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) +
                  cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      * u / uav_mass_;
    /* d v_y = (-cos y * sin r + sin y * sin p * cos r) * (u1 + u2 + u3 + u4) / m  */
    dev_x(V_Y) = (-cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) +
                  sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      * u / uav_mass_;
    /* d v_z = (cos p * cos r) * (u1 + u2 + u3 + u4) / m */
    dev_x(V_Z) = (cos((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      * u / uav_mass_ - 9.78;

    /* e_r, e_p, e_y */
    /* d e = R_e * w_b */
    Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd R_e = MatrixXd::Zero(3, 3);
    R_e << 1, tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]), tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]), -sin((*x_ptr)[E_R]),
      0, sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    Vector3d d_e = R_e * w;
    for (int i = E_R; i <= E_Y; ++i)
      dev_x(i) = d_e(i - E_R);

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + M_para * [u1;u2;u3;u4]), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    MatrixXd w_m = MatrixXd::Zero(3, 3);
    w_m << 0, -(*x_ptr)[W_Z], (*x_ptr)[W_Y],
      (*x_ptr)[W_Z], 0, -(*x_ptr)[W_X],
      -(*x_ptr)[W_Y], (*x_ptr)[W_X], 0;
    Vector3d dw;
    dw = I_ptr_->inverse() * (-w_m * ((*I_ptr_) * w) + (*M_para_ptr_) * ((*u_ptr) + (*un_ptr_)));
    // test: real u
    // dw = I_ptr_->inverse() * (-w_m * ((*I_ptr_) * w) + (*M_para_ptr_) * (*u_ptr));
    for (int i = 0; i < 3; ++i)
      dev_x(W_X + i) = dw(i);

    dev_x = dev_x / control_freq_;
    *new_x_ptr = dev_x + *x_ptr;
    // test
    //*new_x_ptr = stateAddition(x_ptr, &dev_x);
  }

  void SlqFiniteDiscreteControlQuadrotor::normalizeQuaternion(VectorXd *new_x_ptr){
    double q_sum = 0.0;
    for (int i = Q_W; i <= Q_Z; ++i)
      q_sum += pow((*new_x_ptr)(i), 2.0);
    q_sum = sqrt(q_sum);
    if (q_sum == 0.0){
      (*new_x_ptr)(Q_W) = 1.0;
      for (int i = Q_X; i <= Q_Z; ++i)
        (*new_x_ptr)(i) = 0.0;
    }
    else{
      for (int i = Q_W; i <= Q_Z; ++i)
        (*new_x_ptr)(i) = (*new_x_ptr)(i) / q_sum;
    }
  }

  double SlqFiniteDiscreteControlQuadrotor::getSystemEnergy(){
    double energy_sum = 0.0;
    VectorXd u = u_vec_[0];
    VectorXd x = x_vec_[0];
    for (int i = 0; i < iteration_times_; ++i){
      // test: add weight for waypoint
      updateQWeight(i * end_time_ / iteration_times_);

      double state_energy = x.transpose() * (*Q_ptr_) * x;
      state_energy *= 0.5;
      Vector4d u_real = u + (*un_ptr_);
      // test: real u
      // Vector4d u_real = u;
      double control_energy = (u_real.transpose() * (*R_ptr_) * u_real);
      control_energy *= 0.5;
      energy_sum += state_energy + control_energy;

      if (quaternion_mode_)
        updateMatrixAB(&x, &u);
      else
        updateEulerMatrixAB(&x, &u);
      VectorXd new_x(x_size_);
      if (quaternion_mode_){
        updateNewState(&new_x, &x, &u);
        normalizeQuaternion(&new_x);
      }
      else
        updateEulerNewState(&new_x, &x, &u);
      x = new_x;
      if (i != iteration_times_ - 1)
        u = u_vec_[i+1];
    }
    // tf
    (*P_ptr_) = MatrixXd::Zero(x_size_, x_size_);
    for (int i = 0; i < 3; ++i)
      (*P_ptr_)(i, i) = 100.0;
    for (int i = 3; i < x_size_; ++i)
      (*P_ptr_)(i, i) = 1.0;
    double state_tf_energy = 0.5 * x.transpose() * (*P_ptr_) * x;
    energy_sum += state_tf_energy;

    if (debug_)
      std::cout << "alpha: " << alpha_ << ",  energy: " << energy_sum << "\n";
    return energy_sum;
  }

  bool SlqFiniteDiscreteControlQuadrotor::feedforwardConverged(){
    double fw_max = 0.0;
    for (int i = 0; i < iteration_times_; ++i){
      double control_sum = 0.0;
      for (int j = 0; j < u_size_; ++j){
        control_sum += pow((u_fw_vec_[i])(j), 2.0);
      }
      control_sum = sqrt(control_sum);
      if (control_sum > fw_max)
        fw_max = control_sum;
    }
    if (fw_max < 0.1)
      return true;
    else
      return false;
  }

  VectorXd SlqFiniteDiscreteControlQuadrotor::stateAddition(VectorXd *x1_ptr, VectorXd *x2_ptr){
    VectorXd result(x_size_);
    result = (*x1_ptr) + (*x2_ptr);
    Vector4d q = quationAddition(Vector4d((*x1_ptr)(Q_W), (*x1_ptr)(Q_X),
                                          (*x1_ptr)(Q_Y), (*x1_ptr)(Q_Z)),
                                 Vector4d((*x2_ptr)(Q_W), (*x2_ptr)(Q_X),
                                          (*x2_ptr)(Q_Y), (*x2_ptr)(Q_Z)));
    for (int i = Q_W; i <= Q_Z; ++i)
      result(i) = q(i - Q_W);

    normalizeQuaternion(&result);
    return result;
  }

  VectorXd SlqFiniteDiscreteControlQuadrotor::stateSubtraction(VectorXd *x1_ptr, VectorXd *x2_ptr){
    VectorXd result(x_size_);
    result = (*x1_ptr) - (*x2_ptr);
    Vector4d q = quationAddition(Vector4d((*x1_ptr)(Q_W), (*x1_ptr)(Q_X),
                                          (*x1_ptr)(Q_Y), (*x1_ptr)(Q_Z)),
                                 Vector4d((*x2_ptr)(Q_W), -(*x2_ptr)(Q_X),
                                          -(*x2_ptr)(Q_Y), -(*x2_ptr)(Q_Z)));
    for (int i = Q_W; i <= Q_Z; ++i)
      result(i) = q(i - Q_W);

    normalizeQuaternion(&result);
    return result;
  }

  Vector4d SlqFiniteDiscreteControlQuadrotor::quationAddition(Vector4d q1, Vector4d q2){
    Vector4d result;
    result(0) = q1(0) * q2(0) - q1(1) * q2(1) -
      q1(2) * q2(2) - q1(3) * q2(3);
    result(1) = q1(0) * q2(1) + q1(1) * q2(0) +
      q1(2) * q2(3) - q1(3) * q2(2);
    result(2) = q1(0) * q2(2) - q1(1) * q2(3) +
      q1(2) * q2(0) + q1(3) * q2(1);
    result(3) = q1(0) * q2(3) + q1(1) * q2(2) -
      q1(2) * q2(1) + q1(3) * q2(0);
    return result;
  }

  VectorXd SlqFiniteDiscreteControlQuadrotor::getAbsoluteState(VectorXd *relative_x_ptr){
    if (quaternion_mode_)
      return stateAddition(relative_x_ptr, xn_ptr_);
    else
      // todo
      return (*relative_x_ptr + *xn_ptr_);
  }

  VectorXd SlqFiniteDiscreteControlQuadrotor::getRelativeState(VectorXd *absolute_x_ptr){
    if (quaternion_mode_)
      return stateSubtraction(absolute_x_ptr, xn_ptr_);
    else
      // todo
      return (*absolute_x_ptr - *xn_ptr_);
  }

  void SlqFiniteDiscreteControlQuadrotor::updateQWeight(double time){
    double rho = 1.0;
    double weight = exp(-rho / 2 * pow(time - end_time_, 2.0));
    for (int j = 0; j < 6; ++j)
      // (*Q_ptr_)(j, j) = 10.0 * weight + 1;
      (*Q_ptr_)(j, j) = 10.0 * weight;
    for (int j = 6; j < x_size_; ++j)
      // (*Q_ptr_)(j, j) = weight + 1;
      (*Q_ptr_)(j, j) = weight;
  }

  void SlqFiniteDiscreteControlQuadrotor::updateSLQEquations(){
    (*H_ptr_) = (*R_ptr_) + B_ptr_->transpose() * (*P_ptr_) * (*B_ptr_);
    (*G_ptr_) = B_ptr_->transpose() * (*P_ptr_) * (*A_ptr_);
    (*g_ptr_) = (*r_ptr_) + B_ptr_->transpose() * (*p_ptr_);
    (*K_ptr_) = -(H_ptr_->inverse() * (*G_ptr_));
    (*l_ptr_) = -(H_ptr_->inverse() * (*g_ptr_));
    (*P_ptr_) = (*Q_ptr_) + A_ptr_->transpose() * (*P_ptr_) * (*A_ptr_)
      + K_ptr_->transpose() * (*H_ptr_) * (*K_ptr_)
      + K_ptr_->transpose() * (*G_ptr_)
      + G_ptr_->transpose() * (*K_ptr_);
    (*p_ptr_) = (*q_ptr_) + A_ptr_->transpose() * (*p_ptr_)
      + K_ptr_->transpose() * (*H_ptr_) * (*l_ptr_)
      + K_ptr_->transpose() * (*g_ptr_)
      + G_ptr_->transpose() * (*l_ptr_);
  }

  void SlqFiniteDiscreteControlQuadrotor::FDLQR(){
    *x_ptr_ = x_vec_[0];
    *u_ptr_ = u_vec_[0];
    if (quaternion_mode_)
      updateMatrixAB(x_ptr_, u_ptr_);
    else
      updateEulerMatrixAB(x_ptr_, u_ptr_);

    std::vector<MatrixXd> F_vec;
    MatrixXd P(x_size_, x_size_);
    P = *Q_ptr_;
    for (int i = 0; i < iteration_times_; ++i){
      MatrixXd F = MatrixXd::Zero(u_size_, x_size_);
      F = ((*R_ptr_) + B_ptr_->transpose() * P * (*B_ptr_)).inverse()
        * (B_ptr_->transpose() * P * (*A_ptr_));
      P = A_ptr_->transpose() * P * (*A_ptr_)
        - (A_ptr_->transpose() * P * (*B_ptr_)) * F
        + (*Q_ptr_);
      F_vec.push_back(F);
    }

    VectorXd x = getRelativeState(x0_ptr_);
    for (int i = iteration_times_ - 1; i >= 0; --i){
      VectorXd u = -F_vec[i] * x;
      VectorXd new_x(x_size_);
      updateEulerNewState(&new_x, &x, &u);
      x = new_x;
      x_vec_[iteration_times_ - i] = x;

      // Guarantee control is in bound
      checkControlInputFeasible(&u);
      u_vec_[iteration_times_ - i] = u;

      if ((i % 50 == 0 || i == iteration_times_ - 1) && debug_){
        printStateInfo(&x, i);
        printControlInfo(&u, i);
      }
    }
    ROS_INFO("[SLQ] LQR init finished");
  }

  void SlqFiniteDiscreteControlQuadrotor::checkControlInputFeasible(VectorXd *u){
    for (int j = 0; j < u_size_; ++j){
      if ((*u)(j) + (*un_ptr_)(j) < uav_rotor_thrust_min_)
        (*u)(j) = uav_rotor_thrust_min_ - (*un_ptr_)(j);
      else if ((*u)(j) + (*un_ptr_)(j) > uav_rotor_thrust_max_)
        (*u)(j) = uav_rotor_thrust_max_ - (*un_ptr_)(j);
    }
  }

  void SlqFiniteDiscreteControlQuadrotor::printStateInfo(VectorXd *x, int id){
    VectorXd new_absolute_x;
    new_absolute_x = getAbsoluteState(x);
    std::cout << "[debug] id[" << id << "]print current state:\n";
    for (int j = 0; j < x_size_; ++j)
      std::cout << new_absolute_x(j) << ", ";
    std::cout << "\n";
  }

  void SlqFiniteDiscreteControlQuadrotor::printControlInfo(VectorXd *u, int id){
    std::cout << "[debug] id[" << id << "]print current u:\n";
    for (int j = 0; j < u_size_; ++j)
      std::cout << (*u)(j) + (*un_ptr_)(j) << ", ";
    std::cout << "\n";
  }

  void SlqFiniteDiscreteControlQuadrotor::printMatrixAB(){
    std::cout << "examine A:";
    for (int i = 0; i < x_size_; ++i){
      std::cout << "\n";
      for (int j = 0; j < x_size_; ++j){
        std::cout << (*A_ptr_)(i, j) << ", ";
      }
      std::cout << ";";
    }
    std::cout << "\n\nexamine B:";
    for (int i = 0; i < x_size_; ++i){
      std::cout << "\n";
      for (int j = 0; j < u_size_; ++j){
        std::cout << (*B_ptr_)(i, j) << ", ";
      }
      std::cout << ";";
    }
    std::cout << "\n\n";
  }

}


