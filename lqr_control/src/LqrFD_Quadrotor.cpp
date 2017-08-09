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

#include <lqr_control/LqrFD_Quadrotor.h>
namespace lqr_discrete{
  void LqrFiniteDiscreteControlQuadrotor::initLQR(double freq, double period, VectorXd *x0, VectorXd *xn){
    control_freq_ = freq;
    if (floor(freq * period) < freq * period){
      end_time_ = (floor(freq * period) + 1.0) / freq;
      iteration_times_ = floor(freq * period) + 1;
    }
    else{
      end_time_ = period;
      iteration_times_ = floor(freq * period);
    }
    x_size_ = 13;
    u_size_ = 4;
    A_ptr_ = new MatrixXd(x_size_, x_size_);
    B_ptr_ = new MatrixXd(x_size_, u_size_);
    Q_ptr_ = new MatrixXd(x_size_, x_size_);
    R_ptr_ = new MatrixXd(u_size_, u_size_);
    x0_ptr_ = new VectorXd(x_size_);
    xn_ptr_ = new VectorXd(x_size_);
    x_ptr_ = new VectorXd(x_size_);
    u_ptr_ = new VectorXd(u_size_);

    for (int i = 0; i < x_size_; ++i)
      (*x0_ptr_)(i) = (*x0)(i);
    for (int i = 0; i < x_size_; ++i)
      (*xn_ptr_)(i) = (*xn)(i);

    /* init Q and R matrice */
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < x_size_; ++j)
        (*Q_ptr_)(i, j) = 0.0;
    for (int i = 0; i < x_size_; ++i)
      (*Q_ptr_)(i, i) = 10.0;

    for (int i = 0; i < u_size_; ++i)
      for (int j = 0; j < u_size_; ++j)
        (*R_ptr_)(i, j) = 0.0;
    for (int i = 0; i < u_size_; ++i)
      (*R_ptr_)(i, i) = 1.0;

    /* uav property from paper eth15-slq-window */
    I_ptr_ = new MatrixXd(3, 3);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        (*I_ptr_)(i, j) = 0.0;
    (*I_ptr_)(0, 0) = 0.03;
    (*I_ptr_)(1, 1) = 0.03;
    (*I_ptr_)(2, 2) = 0.05;

    uav_mass_ = 0.5;
    M_para_ptr_ = new MatrixXd(3, 4);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 4; ++j)
        (*M_para_ptr_)(i, j) = 0.0;
    double r_uav = 0.3, c_rf = 0.025;
    (*M_para_ptr_)(0, 1) = -r_uav;
    (*M_para_ptr_)(0, 3) = r_uav;
    (*M_para_ptr_)(1, 0) = r_uav;
    (*M_para_ptr_)(1, 2) = -r_uav;
    (*M_para_ptr_)(2, 0) = (*M_para_ptr_)(2, 2) = -c_rf;
    (*M_para_ptr_)(2, 1) = (*M_para_ptr_)(2, 3) = c_rf;

    /* Assume initial and final state is still, namely dx = [v, a] = 0 */
    u0_ptr_ = new VectorXd(u_size_);
    un_ptr_ = new VectorXd(u_size_);
    for (int i = 0; i < 4; ++i)
      (*u0_ptr_)(i) = (*un_ptr_)(i) = uav_mass_ * 9.78 / 4.0;

    debug_ = true;
  }

  void LqrFiniteDiscreteControlQuadrotor::updateAll(){
    MatrixXd P = MatrixXd::Zero(x_size_, x_size_);
    // todo: assume N is zero, namely do not have x^T *N*u in cost function
    MatrixXd N = MatrixXd::Zero(x_size_, u_size_);
    P = *Q_ptr_;
    *x_ptr_ = *x0_ptr_;

    // todo: assume start point the quadrotor is hovering
    for (int i = 0; i < 4; ++i)
      (*u_ptr_)(i) = uav_mass_ * 9.78 / 4.0;
    updateMatrixAB();
    for (int i = 0; i < iteration_times_; ++i){
      MatrixXd F = MatrixXd::Zero(u_size_, x_size_);
      F = ((*R_ptr_) + B_ptr_->transpose() * P * (*B_ptr_)).inverse()
        * (B_ptr_->transpose() * P * (*A_ptr_) + N.transpose());
      P = A_ptr_->transpose() * P * (*A_ptr_)
        - (A_ptr_->transpose() * P * (*B_ptr_) + N) * F
        + (*Q_ptr_);
      F_vec_.push_back(F);
    }

    // debug
    if (debug_){
      std::cout << "\n\nexamine A:";
      for (int i = 0; i < x_size_; ++i){
        std::cout << "\n";
        for (int j = 0; j < x_size_; ++j){
          std::cout << (*A_ptr_)(i, j) << ", ";
        }
      }
      std::cout << "\n\nexamine B:";
      for (int i = 0; i < x_size_; ++i){
        std::cout << "\n";
        for (int j = 0; j < u_size_; ++j){
          std::cout << (*B_ptr_)(i, j) << ", ";
        }
      }
    }

    x_vec_.push_back(*x0_ptr_);
    u_vec_.push_back(*u0_ptr_);
    for (int i = 0; i < iteration_times_; ++i){
      VectorXd u = -F_vec_[iteration_times_-1-i] * (*x_ptr_);
      // debug
      if (debug_){
        if (i % 100 == 0){
          std::cout << "\nexamine x:\n";
          for (int j = 0; j < x_size_; ++j)
            std::cout << (*x_ptr_)(j) << ", ";
          std::cout << "\nexamine u:\n";
          for (int j = 0; j < u_size_; ++j)
            std::cout << u(j) << ", ";
          std::cout << "\n\nexamine F:";
          for (int k = 0; k < u_size_; ++k){
            std::cout << "\n";
            for (int j = 0; j < x_size_; ++j){
              std::cout << F_vec_[iteration_times_-1-i](k, j) << ", ";
            }
          }
        }
      }

      u_vec_.push_back(u);
      // method 1: x(k+1) = A*x(k) + B*u(k)
      VectorXd x = (*A_ptr_) * (*x_ptr_)
        + (*B_ptr_) * u;
      x_vec_.push_back(x);
      *x_ptr_ = x;

      // method 2: x(k+1) = f(x(k), u(k))
      // VectorXd *x = new VectorXd(x_size_);
      // updateNewState(x);
      // x_vec_.push_back(x);
      // *x_ptr_ = *x;

    }
  }

  void LqrFiniteDiscreteControlQuadrotor::updateMatrixAB(){
    updateMatrixA();
    updateMatrixB();
  }

  void LqrFiniteDiscreteControlQuadrotor::updateMatrixA(){
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < x_size_; ++j){
        (*A_ptr_)(i, j) = 0.0;
      }
    /* x, y, z */
    (*A_ptr_)(P_X, V_X) = 1;
    (*A_ptr_)(P_Y, V_Y) = 1;
    (*A_ptr_)(P_Z, V_Z) = 1;

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += (*u_ptr_)[i];
    /* u' = u / m */
    u = u / uav_mass_;
    /* d v_x = (2 * q_w * q_y + 2 * q_x * q_z) * u' */
    (*A_ptr_)(V_X, Q_W) = 2 * (*x_ptr_)[Q_Y] * u;
    (*A_ptr_)(V_X, Q_X) = 2 * (*x_ptr_)[Q_Z] * u;
    (*A_ptr_)(V_X, Q_Y) = 2 * (*x_ptr_)[Q_W] * u;
    (*A_ptr_)(V_X, Q_Z) = 2 * (*x_ptr_)[Q_X] * u;
    /* d v_y = (2 * q_w * q_x + 2 * q_y * q_z) * u' */
    (*A_ptr_)(V_Y, Q_W) = 2 * (*x_ptr_)[Q_X] * u;
    (*A_ptr_)(V_Y, Q_X) = 2 * (*x_ptr_)[Q_W] * u;
    (*A_ptr_)(V_Y, Q_Y) = 2 * (*x_ptr_)[Q_Z] * u;
    (*A_ptr_)(V_Y, Q_Z) = 2 * (*x_ptr_)[Q_Y] * u;
    /* d v_z = (1 - 2 * q_x * q_x - 2 * q_y * q_y) * u' */
    (*A_ptr_)(V_Z, Q_X) = -2 * (*x_ptr_)[Q_X] * u;
    (*A_ptr_)(V_Z, Q_Y) = -2 * (*x_ptr_)[Q_Y] * u;

    /* q_w, q_x, q_y, q_z */
    /* d q = 1/2 * q * [0, w]^T (the multiply opeation is under quaternion multiply) */
    /* d q_w = 1/2 (q_w * 0 - q_x * w_x - q_y * w_y - q_z * w_z) */
    (*A_ptr_)(Q_W, Q_X) = - (*x_ptr_)[W_X] / 2.0;
    (*A_ptr_)(Q_W, Q_Y) = - (*x_ptr_)[W_Y] / 2.0;
    (*A_ptr_)(Q_W, Q_Z) = - (*x_ptr_)[W_Z] / 2.0;
    /* d q_x = 1/2 (q_w * w_x + q_x * 0 + q_y * w_z - q_z * w_y) */
    (*A_ptr_)(Q_X, Q_W) = (*x_ptr_)[W_X] / 2.0;
    (*A_ptr_)(Q_X, Q_Y) = (*x_ptr_)[W_Z] / 2.0;
    (*A_ptr_)(Q_X, Q_Z) = - (*x_ptr_)[W_Y] / 2.0;
    /* d q_y = 1/2 (q_w * w_y - q_x * w_z + q_y * 0 + q_z * w_x) */
    (*A_ptr_)(Q_Y, Q_W) = (*x_ptr_)[W_Y] / 2.0;
    (*A_ptr_)(Q_Y, Q_X) = - (*x_ptr_)[W_Z] / 2.0;
    (*A_ptr_)(Q_Y, Q_Z) = (*x_ptr_)[W_X] / 2.0;
    /* d q_z = 1/2 (q_w * w_z + q_x * w_y - q_y * w_x + q_z * 0) */
    (*A_ptr_)(Q_Z, Q_W) = (*x_ptr_)[W_Z] / 2.0;
    (*A_ptr_)(Q_Z, Q_X) = (*x_ptr_)[W_Y] / 2.0;
    (*A_ptr_)(Q_Z, Q_Y) = (*x_ptr_)[W_X] / 2.0;

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + tau), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    /* d w_w = I^-1 * (- d(w^) * (Iw) - (w^) * (I * d(w))) */
    Vector3d w((*x_ptr_)[W_X], (*x_ptr_)[W_Y], (*x_ptr_)[W_Z]);
    MatrixXd w_m = MatrixXd::Zero(3, 3);
    w_m(0, 1) = -(*x_ptr_)[W_Z]; w_m(0, 2) = (*x_ptr_)[W_Y];
    w_m(1, 0) = (*x_ptr_)[W_Z]; w_m(1, 2) = -(*x_ptr_)[W_X];
    w_m(2, 0) = -(*x_ptr_)[W_Y]; w_m(2, 1) = (*x_ptr_)[W_X];
    MatrixXd dw_m = MatrixXd::Zero(3, 3);
    dw_m(1, 2) = -1;
    dw_m(2, 1) = 1;
    Vector3d dw_x = I_ptr_->inverse() * (-dw_m * ((*I_ptr_) * w))
      - w_m * ((*I_ptr_) * Vector3d(1.0, 0.0, 0.0));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_X) = dw_x(i);

    dw_m = MatrixXd::Zero(3, 3);
    dw_m(0, 2) = 1;
    dw_m(2, 0) = -1;
    Vector3d dw_y = I_ptr_->inverse() * (-dw_m * ((*I_ptr_) * w))
      - w_m * ((*I_ptr_) * Vector3d(0.0, 1.0, 0.0));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_Y) = dw_y(i);

    dw_m = MatrixXd::Zero(3, 3);
    dw_m(0, 1) = -1;
    dw_m(1, 0) = 1;
    Vector3d dw_z = I_ptr_->inverse() * (-dw_m * ((*I_ptr_) * w))
      - w_m * ((*I_ptr_) * Vector3d(0.0, 0.0, 1.0));
    for (int i = 0; i < 3; ++i)
      (*A_ptr_)(W_X + i, W_Z) = dw_z(i);

    (*A_ptr_) = (*A_ptr_) / control_freq_ + MatrixXd::Identity(x_size_, x_size_);
  }

  void LqrFiniteDiscreteControlQuadrotor::updateMatrixB(){
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < u_size_; ++j){
        (*B_ptr_)(i, j) = 0.0;
      }
    /* x, y, z */
    /* all 0 */

    /* v_x, v_y, v_z */
    /* d v_x = (2 * q_w * q_y + 2 * q_x * q_z) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_X, U_1) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_Y] + 2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_X, U_2) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_Y] + 2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_X, U_3) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_Y] + 2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_X, U_4) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_Y] + 2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_Z]) / uav_mass_;
    /* d v_y = (2 * q_w * q_x + 2 * q_y * q_z) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_Y, U_1) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_X] + 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_Y, U_2) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_X] + 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_Y, U_3) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_X] + 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_Y, U_4) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_X] + 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Z]) / uav_mass_;
    /* d v_z = (1 - 2 * q_x * q_x - 2 * q_y * q_y) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_Z, U_1) = (1 -2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_X] - 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Y]) / uav_mass_;
    (*B_ptr_)(V_Z, U_2) = (1 -2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_X] - 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Y]) / uav_mass_;
    (*B_ptr_)(V_Z, U_3) = (1 -2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_X] - 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Y]) / uav_mass_;
    (*B_ptr_)(V_Z, U_4) = (1 -2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_X] - 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Y]) / uav_mass_;

    /* q_w, q_x, q_y, q_z */
    /* d q = 1/2 * q * [0, w]^T */
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
  void LqrFiniteDiscreteControlQuadrotor::updateNewState(VectorXd *new_x_ptr){
    VectorXd dev_x(x_size_);
    for (int i = 0; i < x_size_; ++i)
      dev_x(i) = 0.0;
    /* x, y, z */
    dev_x(P_X) = (*x_ptr_)(V_X);
    dev_x(P_Y) = (*x_ptr_)(V_Y);
    dev_x(P_Z) = (*x_ptr_)(V_Z);

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += (*u_ptr_)[i];
    /* u' = u / m */
    u = u / uav_mass_;
    /* d v_x = (2 * q_w * q_y + 2 * q_x * q_z) * u' */
    dev_x(V_X) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_Y]
                  + 2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_Z]) * u;
    /* d v_y = (2 * q_w * q_x + 2 * q_y * q_z) * u' */
    dev_x(V_Y) = (2 * (*x_ptr_)[Q_W] * (*x_ptr_)[Q_X]
                  + 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Z]) * u;
    /* d v_z = (1 - 2 * q_x * q_x - 2 * q_y * q_y) * u' */
    dev_x(V_Z) = (1 - 2 * (*x_ptr_)[Q_X] * (*x_ptr_)[Q_X]
                  - 2 * (*x_ptr_)[Q_Y] * (*x_ptr_)[Q_Y]) * u - 9.78;

    /* q_w, q_x, q_y, q_z */
    /* d q = 1/2 * q * [0, w]^T (the multiply opeation is under quaternion multiply) */
    /* d q_w = 1/2 (q_w * 0 - q_x * w_x - q_y * w_y - q_z * w_z) */
    dev_x(Q_W) = (-(*x_ptr_)[Q_X] * (*x_ptr_)[W_X]
                  -(*x_ptr_)[Q_Y] * (*x_ptr_)[W_Y]
                  -(*x_ptr_)[Q_Z] * (*x_ptr_)[W_Z])/ 2.0;
    /* d q_x = 1/2 (q_w * w_x + q_x * 0 + q_y * w_z - q_z * w_y) */
    dev_x(Q_X) = ((*x_ptr_)[Q_W] * (*x_ptr_)[W_X]
                  +(*x_ptr_)[Q_Y] * (*x_ptr_)[W_Z]
                  -(*x_ptr_)[Q_Z] * (*x_ptr_)[W_Y])/ 2.0;
    /* d q_y = 1/2 (q_w * w_y - q_x * w_z + q_y * 0 + q_z * w_x) */
    dev_x(Q_Y) = ((*x_ptr_)[Q_W] * (*x_ptr_)[W_Y]
                  -(*x_ptr_)[Q_X] * (*x_ptr_)[W_Z]
                  +(*x_ptr_)[Q_Z] * (*x_ptr_)[W_X])/ 2.0;
    /* d q_z = 1/2 (q_w * w_z + q_x * w_y - q_y * w_x + q_z * 0) */
    dev_x(Q_Z) = ((*x_ptr_)[Q_W] * (*x_ptr_)[W_Z]
                  +(*x_ptr_)[Q_X] * (*x_ptr_)[W_Y]
                  -(*x_ptr_)[Q_Y] * (*x_ptr_)[W_X])/ 2.0;

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + M_para * [u1;u2;u3;u4]), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    /* d w_w = I^-1 * (- d(w^) * (Iw) - (w^) * (I * d(w))) */
    Vector3d w((*x_ptr_)[W_X], (*x_ptr_)[W_Y], (*x_ptr_)[W_Z]);
    MatrixXd w_m = MatrixXd::Zero(3, 3);
    w_m(0, 1) = -(*x_ptr_)[W_Z]; w_m(0, 2) = (*x_ptr_)[W_Y];
    w_m(1, 0) = (*x_ptr_)[W_Z]; w_m(1, 2) = -(*x_ptr_)[W_X];
    w_m(2, 0) = -(*x_ptr_)[W_Y]; w_m(2, 1) = (*x_ptr_)[W_X];
    Vector3d dw;
    dw = I_ptr_->inverse() * (-w_m * (*I_ptr_) * w + (*M_para_ptr_) * (*u_ptr_));
    for (int i = 0; i < 3; ++i)
      dev_x(W_X + i) = dw(i);

    *new_x_ptr = dev_x / control_freq_ + *x_ptr_;
  }
}


