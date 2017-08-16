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
    lqr_controller_ptr_ = new LqrFiniteDiscreteControlQuadrotor(nh_, nhp_);
    lqr_controller_ptr_->initLQR(freq, period, x0, xn);

    control_freq_ = freq;
    if (floor(freq * period) < freq * period){
      end_time_ = (floor(freq * period) + 1.0) / freq;
      iteration_times_ = floor(freq * period) + 1;
    }
    else{
      end_time_ = period;
      iteration_times_ = floor(freq * period);
    }
    std::cout << "[SLQ] Trajectory period: " << period
              << ", Itetation times: " << iteration_times_ << "\n";
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
    for (int i = 0; i < 6; ++i)
      (*Q_ptr_)(i, i) = 10.0;
    for (int i = 6; i < x_size_; ++i)
      (*Q_ptr_)(i, i) = 1.0;
    for (int i = 0; i < u_size_; ++i)
      for (int j = 0; j < u_size_; ++j)
        (*R_ptr_)(i, j) = 0.0;
    for (int i = 0; i < u_size_; ++i)
      (*R_ptr_)(i, i) = 10000.0;

    /* uav property from paper eth15-slq-window */
    I_ptr_ = new MatrixXd(3, 3);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        (*I_ptr_)(i, j) = 0.0;
    (*I_ptr_)(0, 0) = 0.03;
    (*I_ptr_)(1, 1) = 0.03;
    (*I_ptr_)(2, 2) = 0.05;

    uav_mass_ = 1.0;
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
    for (int i = 0; i < 4; ++i){
      (*u0_ptr_)(i) = uav_mass_ * 9.78 / 4.0;
      (*un_ptr_)(i) = uav_mass_ * 9.78 / 4.0;
    }

    /* SLQ special initialization */
    // todo: assume start point the quadrotor is hovering
    VectorXd x_init(x_size_), u_init(u_size_);
    x_init = getRelativeState(x0_ptr_);
    u_init = VectorXd::Zero(u_size_);
    // for (int i = 0; i < u_size_; ++i)
    //   u_init(i) = uav_mass_ * 9.78 / 4.0;

    for (int i = 0; i < iteration_times_; ++i){
      x_vec_.push_back(x_init);
      u_vec_.push_back(u_init);
    }

    // todo: initialization these matrix
    P_ptr_ = new MatrixXd(x_size_, x_size_);
    *P_ptr_ = *Q_ptr_;

    p_ptr_ = new VectorXd(x_size_);
    // x^T * p. initial value depends on system paramater
    (*p_ptr_) = VectorXd::Zero(x_size_);

    H_ptr_ = new MatrixXd(u_size_, u_size_);
    // method 1:
    (*H_ptr_) = MatrixXd::Zero(u_size_, u_size_);
    // method 2:
    // lqr_controller_ptr_->updateAll();
    // std::cout << "\n[SLQ]: LQR initial value is generated.\n";
    // getRicattiH();
    // std::cout << "\n[SLQ]: Ricatti matrice is generated.\n";

    G_ptr_ = new MatrixXd(u_size_, x_size_);
    (*G_ptr_) = MatrixXd::Zero(u_size_, x_size_);

    K_ptr_ = new MatrixXd(u_size_, x_size_);
    (*K_ptr_) = MatrixXd::Zero(u_size_, x_size_);

    g_ptr_ = new VectorXd(u_size_);
    (*g_ptr_) = VectorXd::Zero(u_size_);

    l_ptr_ = new VectorXd(u_size_);
    (*l_ptr_) = VectorXd::Zero(u_size_);

    r_ptr_ = new VectorXd(u_size_);
    (*r_ptr_) = VectorXd::Zero(u_size_);
    // test: u is (u_uav - un)
    // control energy: (u + un)' * R * (u + un) = u'Ru + 2u'R un
    (*r_ptr_) = (*R_ptr_) * (*un_ptr_);

    alpha_ = 1.0;

    debug_ = true;
  }

  void SlqFiniteDiscreteControlQuadrotor::getRicattiH(){
    VectorXd x_tf(x_size_);
    VectorXd u_tf(u_size_);
    std::cout << lqr_controller_ptr_->x_vec_.size() << "]\n";
    x_tf = lqr_controller_ptr_->x_vec_[lqr_controller_ptr_->x_vec_.size() - 1];
    u_tf = lqr_controller_ptr_->u_vec_[lqr_controller_ptr_->u_vec_.size() - 1];
    LqrInfiniteDiscreteControlQuadrotor *lqr_ID_controller_ptr;
    lqr_ID_controller_ptr = new LqrInfiniteDiscreteControlQuadrotor(nh_, nhp_);
    lqr_ID_controller_ptr->initLQR(control_freq_, iteration_times_, x0_ptr_, xn_ptr_);
    if (debug_){
      std::cout << "\nx and u in tf:\n";
      for (int i = 0; i < x_size_; ++i)
        std::cout << x_tf(i) << ", ";
      std::cout << "\n";
      for (int i = 0; i < u_size_; ++i)
        std::cout << u_tf(i) << ", ";
      std::cout << "\n";
    }
    *(lqr_ID_controller_ptr->x_ptr_) = x_tf;
    *(lqr_ID_controller_ptr->u_ptr_) = u_tf;
    lqr_ID_controller_ptr->updateMatrixAB();
    lqr_ID_controller_ptr->ricatti(H_ptr_);
    if (debug_){
      std::cout << "\n\nH matrice:";
      for (int i = 0; i < x_size_; ++i){
        std::cout << "\n";
        for (int j = 0; j < x_size_; ++j){
          std::cout << (*H_ptr_)(i, j) << ", ";
        }
      }
    }
  }

  void SlqFiniteDiscreteControlQuadrotor::updateAll(){
  }

  void SlqFiniteDiscreteControlQuadrotor::iterativeOptimization(){
    *P_ptr_ = *Q_ptr_;
    *p_ptr_ = VectorXd::Zero(x_size_);

    // test: ARE from matlab
    // for (int i = 0; i < x_size_; ++i)
    //   (*P_ptr_)(i, i) = 100;
    // (*P_ptr_)(6, 6) = 37142;
    // (*P_ptr_)(7, 7) = 38643;
    // (*P_ptr_)(8, 8) = 38643;
    // (*P_ptr_)(9, 9) = 1;
    // (*P_ptr_)(10, 10) = 9662;
    // (*P_ptr_)(11, 11) = 9662;
    // (*P_ptr_)(12, 12) = 1;

    std::vector<VectorXd> u_fw_vec;
    std::vector<VectorXd> u_fb_vec;

    for (int i = iteration_times_ - 1; i >= 0; --i){
      // test: add weight for waypoint
      double ru = 1.0;
      double weight = exp(-ru / 2 * pow(i * 5.0 / double(iteration_times_) - 5.0, 2.0));
      for (int j = 0; j < 6; ++j)
        (*Q_ptr_)(j, j) = 10.0 * weight + 1;
      for (int j = 6; j < x_size_; ++j)
        (*Q_ptr_)(j, j) = weight + 1;

      *x_ptr_ = x_vec_[i];
      *u_ptr_ = u_vec_[i];
      updateMatrixAB(x_ptr_, u_ptr_);
      (*H_ptr_) = (*R_ptr_) + B_ptr_->transpose() * (*P_ptr_) * (*B_ptr_);
      (*G_ptr_) = B_ptr_->transpose() * (*P_ptr_) * (*A_ptr_);
      (*g_ptr_) = (*r_ptr_) + B_ptr_->transpose() * (*p_ptr_);
      (*K_ptr_) = -(H_ptr_->inverse() * (*G_ptr_));
      (*l_ptr_) = -(H_ptr_->inverse() * (*g_ptr_));
      (*P_ptr_) = (*Q_ptr_) + A_ptr_->transpose() * (*P_ptr_) * (*A_ptr_)
        + K_ptr_->transpose() * (*H_ptr_) * (*K_ptr_)
        + K_ptr_->transpose() * (*G_ptr_)
        + G_ptr_->transpose() * (*K_ptr_);
      (*p_ptr_) = VectorXd::Zero(x_size_) + A_ptr_->transpose() * (*p_ptr_)
        + K_ptr_->transpose() * (*H_ptr_) * (*l_ptr_)
        + K_ptr_->transpose() * (*g_ptr_)
        + G_ptr_->transpose() * (*l_ptr_);

      u_fb_vec.push_back((*K_ptr_) * (*x_ptr_));
      u_fw_vec.push_back((*l_ptr_));
    }

    // todo
    alpha_ = 1.0;
    double alpha_candidate = 1.0, energy_min = -1.0;
    if (feedforwardConverged(u_fw_vec))
      std::cout << "[SLQ] feedforward converge.";
    else{
      while (1){
        bool u_flag = true;
        for (int i = 0; i < iteration_times_; ++i){
          VectorXd new_u = u_vec_[i] + alpha_ * u_fw_vec[iteration_times_ - 1 - i]
            + u_fb_vec[iteration_times_ - 1 - i];
          // test: u is (u_uav - un)
          new_u = new_u + (*un_ptr_);

          for (int j = 0; j < u_size_; ++j){
            if (new_u(j) < 0 || new_u(j) > uav_mass_ * 9.78 / 4.0 * 3.0){
              u_flag = false;
              break;
            }
          }
          if (!u_flag){
            std::cout << "alpha: " << alpha_ << ", out of dynamic limitation\n";
            break;
          }
        }
        if (alpha_ < 1.0/32.0)
          break;
        else if (u_flag){
          double energy = getSystemEnergy(u_fw_vec, u_fb_vec);
          if (energy_min < 0 || energy < energy_min){
            energy_min = energy;
            alpha_candidate = alpha_;
          }
        }
        alpha_ /= 2.0;
      }
    }

    std::cout << "\nAlpha selected: " << alpha_candidate << "\n\n";

    for (int i = 0; i < iteration_times_; ++i){
      u_vec_[i] = u_vec_[i] + alpha_candidate * u_fw_vec[iteration_times_ - 1 - i]
        + u_fb_vec[iteration_times_ - 1 - i];
    }

    // u is updated in backward way
    *u_ptr_ = u_vec_[0];
    *x_ptr_ = x_vec_[0];
    for (int i = 0; i < iteration_times_ - 1; ++i){
      updateMatrixAB(x_ptr_, u_ptr_);
      VectorXd new_x(x_size_);
      // method 1:
      updateNewState(&new_x, x_ptr_, u_ptr_);
      // method 2:
      // new_x = (*A_ptr_) * x_vec_[i] + (*B_ptr_) * *u_ptr_;
      normalizeQuaternion(&new_x);

      if ((i+1) % 100 <= 2){
        if (debug_){
          VectorXd new_absolute_x = getAbsoluteState(&new_x);
          std::cout << "\n\n[debug] id[" << i << "]print current state:\n";
          for (int j = 0; j < x_size_; ++j)
            std::cout << new_absolute_x(j) << ", ";
          std::cout << "\n[debug] id[" << i << "]print current u:\n";
          for (int j = 0; j < u_size_; ++j)
            //std::cout << (*u_ptr_)(j) << ", ";
            // test: u is (u_uav - un)
            std::cout << (*u_ptr_)(j) + (*un_ptr_)(j) << ", ";
          std::cout << "\n";
        }
      }

      *u_ptr_ = u_vec_[i + 1];
      x_vec_[i + 1] = new_x;
      *x_ptr_ = x_vec_[i + 1];


      // u_vec_[i] = u_vec_[i] + alpha_ * (*l_ptr_) + (*K_ptr_) * (*x_ptr_ - VectorXd::Zero(x_size_));
      // // method 1:
      // updateNewState(x_ptr_);
      // // method 2:
      // // x_vec_[i] = (*A_ptr_) * x_vec_[i] + (*B_ptr_) * *u_ptr_;
      // normalizeQuaternion(x_ptr_);
      // if (i != iteration_times_ - 1){
      //   *u_ptr_ = u_vec_[i + 1];
      //   x_vec_[i + 1] = *x_ptr_;
      // }

    }

    // test output A and B
    VectorXd x = getRelativeState(xn_ptr_);
    VectorXd u = VectorXd::Zero(u_size_);
    updateMatrixAB(&x, &u);
    std::cout << "\n\nexamine A:";
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
    
  }

  void SlqFiniteDiscreteControlQuadrotor::updateMatrixAB(VectorXd *x_ptr, VectorXd *u_ptr){
    updateMatrixA(x_ptr, u_ptr);
    updateMatrixB(x_ptr, u_ptr);
  }

  void SlqFiniteDiscreteControlQuadrotor::updateMatrixA(VectorXd *x_ptr, VectorXd *u_ptr){
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
      u += (*u_ptr)[i];
    // test: u is (u_uav - un)
    for (int i = 0; i < u_size_; ++i)
      u += (*un_ptr_)[i];

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
    /* d q = 1/2 * q * [0, w]^T (the multiply opeation is under quaternion multiply) */
    /* d q_w = 1/2 (q_w * 0 - q_x * w_x - q_y * w_y - q_z * w_z) */
    (*A_ptr_)(Q_W, Q_X) = - (*x_ptr)[W_X] / 2.0;
    (*A_ptr_)(Q_W, Q_Y) = - (*x_ptr)[W_Y] / 2.0;
    (*A_ptr_)(Q_W, Q_Z) = - (*x_ptr)[W_Z] / 2.0;
    (*A_ptr_)(Q_W, W_X) = - (*x_ptr)[Q_X] / 2.0;
    (*A_ptr_)(Q_W, W_Y) = - (*x_ptr)[Q_Y] / 2.0;
    (*A_ptr_)(Q_W, W_Z) = - (*x_ptr)[Q_Z] / 2.0;
    /* d q_x = 1/2 (q_w * w_x + q_x * 0 + q_y * w_z - q_z * w_y) */
    (*A_ptr_)(Q_X, Q_W) = (*x_ptr)[W_X] / 2.0;
    (*A_ptr_)(Q_X, Q_Y) = (*x_ptr)[W_Z] / 2.0;
    (*A_ptr_)(Q_X, Q_Z) = - (*x_ptr)[W_Y] / 2.0;
    (*A_ptr_)(Q_X, W_X) = (*x_ptr)[Q_W] / 2.0;
    (*A_ptr_)(Q_X, W_Y) = -(*x_ptr)[Q_Z] / 2.0;
    (*A_ptr_)(Q_X, W_Z) = (*x_ptr)[Q_Y] / 2.0;
    /* d q_y = 1/2 (q_w * w_y - q_x * w_z + q_y * 0 + q_z * w_x) */
    (*A_ptr_)(Q_Y, Q_W) = (*x_ptr)[W_Y] / 2.0;
    (*A_ptr_)(Q_Y, Q_X) = - (*x_ptr)[W_Z] / 2.0;
    (*A_ptr_)(Q_Y, Q_Z) = (*x_ptr)[W_X] / 2.0;
    (*A_ptr_)(Q_Y, W_X) = (*x_ptr)[Q_Z] / 2.0;
    (*A_ptr_)(Q_Y, W_Y) = (*x_ptr)[Q_W] / 2.0;
    (*A_ptr_)(Q_Y, W_Z) = - (*x_ptr)[Q_X] / 2.0;
    /* d q_z = 1/2 (q_w * w_z + q_x * w_y - q_y * w_x + q_z * 0) */
    (*A_ptr_)(Q_Z, Q_W) = (*x_ptr)[W_Z] / 2.0;
    (*A_ptr_)(Q_Z, Q_X) = (*x_ptr)[W_Y] / 2.0;
    (*A_ptr_)(Q_Z, Q_Y) = (*x_ptr)[W_X] / 2.0;
    (*A_ptr_)(Q_Z, W_X) = - (*x_ptr)[Q_Y] / 2.0;
    (*A_ptr_)(Q_Z, W_Y) = (*x_ptr)[Q_X] / 2.0;
    (*A_ptr_)(Q_Z, W_Z) = (*x_ptr)[Q_W] / 2.0;

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + tau), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    /* d w_w = I^-1 * (- d(w^) * (Iw) - (w^) * (I * d(w))) */
    Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd w_m = MatrixXd::Zero(3, 3);
    w_m(0, 1) = -(*x_ptr)[W_Z]; w_m(0, 2) = (*x_ptr)[W_Y];
    w_m(1, 0) = (*x_ptr)[W_Z]; w_m(1, 2) = -(*x_ptr)[W_X];
    w_m(2, 0) = -(*x_ptr)[W_Y]; w_m(2, 1) = (*x_ptr)[W_X];
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
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < u_size_; ++j){
        (*B_ptr_)(i, j) = 0.0;
      }
    /* x, y, z */
    /* all 0 */

    /* v_x, v_y, v_z */
    /* d v_x = (2 * q_w * q_y + 2 * q_x * q_z) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_X, U_1) = (2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_Y] + 2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_X, U_2) = (2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_Y] + 2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_X, U_3) = (2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_Y] + 2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_Z]) / uav_mass_;
    (*B_ptr_)(V_X, U_4) = (2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_Y] + 2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_Z]) / uav_mass_;
    /* d v_y = (2 * q_y * q_z - 2 * q_w * q_x) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_Y, U_1) = (2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Z] - 2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_X]) / uav_mass_;
    (*B_ptr_)(V_Y, U_2) = (2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Z] - 2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_X]) / uav_mass_;
    (*B_ptr_)(V_Y, U_3) = (2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Z] - 2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_X]) / uav_mass_;
    (*B_ptr_)(V_Y, U_4) = (2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Z] - 2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_X]) / uav_mass_;
    /* d v_z = (1 - 2 * q_x * q_x - 2 * q_y * q_y) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_Z, U_1) = (1 -2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_X] - 2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Y]) / uav_mass_;
    (*B_ptr_)(V_Z, U_2) = (1 -2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_X] - 2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Y]) / uav_mass_;
    (*B_ptr_)(V_Z, U_3) = (1 -2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_X] - 2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Y]) / uav_mass_;
    (*B_ptr_)(V_Z, U_4) = (1 -2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_X] - 2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Y]) / uav_mass_;

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

  void SlqFiniteDiscreteControlQuadrotor::updateNewState(VectorXd *new_x_ptr, VectorXd *x_ptr, VectorXd *u_ptr){
    VectorXd dev_x = VectorXd::Zero(x_size_);
    /* x, y, z */
    dev_x(P_X) = (*x_ptr)(V_X);
    dev_x(P_Y) = (*x_ptr)(V_Y);
    dev_x(P_Z) = (*x_ptr)(V_Z);

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += (*u_ptr)[i];
    // test: u is (u_uav - un)
    for (int i = 0; i < u_size_; ++i)
      u += (*un_ptr_)[i];

    /* u' = u / m */
    u = u / uav_mass_;
    /* d v_x = (2 * q_w * q_y + 2 * q_x * q_z) * u' */
    dev_x(V_X) = (2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_Y]
                  + 2 * (*x_ptr)[Q_X] * (*x_ptr)[Q_Z]) * u;
    /* d v_y = (2 * q_y * q_z - 2 * q_w * q_x) * u' */
    dev_x(V_Y) = (-2 * (*x_ptr)[Q_W] * (*x_ptr)[Q_X]
                  + 2 * (*x_ptr)[Q_Y] * (*x_ptr)[Q_Z]) * u;
    /* d v_z = (q_w ^2 - q_x ^2 - q_y ^2 + q_z ^2) * u' */
    dev_x(V_Z) = (pow((*x_ptr)[Q_W], 2.0) - pow((*x_ptr)[Q_X], 2.0)
                      - pow((*x_ptr)[Q_Y], 2.0) + pow((*x_ptr)[Q_Z], 2.0))
                  * u - 9.78;

    /* q_w, q_x, q_y, q_z */
    /* d q = 1/2 * q * [0, w]^T (the multiply opeation is under quaternion multiply) */
    /* d q_w = 1/2 (q_w * 0 - q_x * w_x - q_y * w_y - q_z * w_z) */
    dev_x(Q_W) = (-(*x_ptr)[Q_X] * (*x_ptr)[W_X]
                  -(*x_ptr)[Q_Y] * (*x_ptr)[W_Y]
                  -(*x_ptr)[Q_Z] * (*x_ptr)[W_Z])/ 2.0;
    /* d q_x = 1/2 (q_w * w_x + q_x * 0 + q_y * w_z - q_z * w_y) */
    dev_x(Q_X) = ((*x_ptr)[Q_W] * (*x_ptr)[W_X]
                  +(*x_ptr)[Q_Y] * (*x_ptr)[W_Z]
                  -(*x_ptr)[Q_Z] * (*x_ptr)[W_Y])/ 2.0;
    /* d q_y = 1/2 (q_w * w_y - q_x * w_z + q_y * 0 + q_z * w_x) */
    dev_x(Q_Y) = ((*x_ptr)[Q_W] * (*x_ptr)[W_Y]
                  -(*x_ptr)[Q_X] * (*x_ptr)[W_Z]
                  +(*x_ptr)[Q_Z] * (*x_ptr)[W_X])/ 2.0;
    /* d q_z = 1/2 (q_w * w_z + q_x * w_y - q_y * w_x + q_z * 0) */
    dev_x(Q_Z) = ((*x_ptr)[Q_W] * (*x_ptr)[W_Z]
                  +(*x_ptr)[Q_X] * (*x_ptr)[W_Y]
                  -(*x_ptr)[Q_Y] * (*x_ptr)[W_X])/ 2.0;

    /* w_x, w_y, w_z */
    /* d w = I^-1 * (- (w^) * (Iw) + M_para * [u1;u2;u3;u4]), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd w_m = MatrixXd::Zero(3, 3);
    w_m(0, 1) = -(*x_ptr)[W_Z]; w_m(0, 2) = (*x_ptr)[W_Y];
    w_m(1, 0) = (*x_ptr)[W_Z]; w_m(1, 2) = -(*x_ptr)[W_X];
    w_m(2, 0) = -(*x_ptr)[W_Y]; w_m(2, 1) = (*x_ptr)[W_X];
    Vector3d dw;
    //dw = I_ptr_->inverse() * (-w_m * (*I_ptr_) * w + (*M_para_ptr_) * (*u_ptr));
    // test: u is (u_uav - un)
    dw = I_ptr_->inverse() * (-w_m * (*I_ptr_) * w + (*M_para_ptr_) * ((*u_ptr) + (*un_ptr_)));
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

  double SlqFiniteDiscreteControlQuadrotor::getSystemEnergy(std::vector<VectorXd> &u_fw_vec, std::vector<VectorXd> &u_fb_vec){
    double energy_sum = 0.0;
    std::vector<VectorXd> u_vec;
    for (int i = 0; i < iteration_times_; ++i){
      u_vec.push_back(u_vec_[i] + alpha_ * u_fw_vec[iteration_times_ - 1 - i]
                      + u_fb_vec[iteration_times_ - 1 - i]);
    }
    VectorXd u = u_vec[0];
    VectorXd x = x_vec_[0];
    for (int i = 0; i < iteration_times_; ++i){
      // test: add weight for waypoint
      double ru = 1.0;
      double weight = exp(-ru / 2 * pow(i * 5.0 / iteration_times_ - 5.0, 2.0));
      MatrixXd Q = MatrixXd::Zero(x_size_, x_size_);
      for (int j = 0; j < 3; ++j)
        Q(j, j) = 100.0 * weight;
      for (int j = 3; j < x_size_; ++j)
        Q(j, j) = weight;
      double state_energy = x.transpose() * Q * x;
      state_energy *= 0.5;
      u = u + (*un_ptr_);
      double control_energy = (u.transpose() * (*R_ptr_) * u);
      control_energy *= 0.5;
      energy_sum += state_energy + control_energy;

      updateMatrixAB(&x, &u);
      VectorXd new_x(x_size_);
      updateNewState(&new_x, &x, &u);
      normalizeQuaternion(&new_x);
      x = new_x;
      if (i != iteration_times_ - 1)
        u = u_vec[i+1];
  }
    // tf
    (*Q_ptr_) = MatrixXd::Zero(x_size_, x_size_);
    for (int i = 0; i < 3; ++i)
      (*Q_ptr_)(i, i) = 100.0;
    for (int i = 3; i < x_size_; ++i)
      (*Q_ptr_)(i, i) = 1.0;
    double state_tf_energy = 0.5 * x.transpose() * (*Q_ptr_) * x;
    energy_sum += state_tf_energy;

    std::cout << "alpha: " << alpha_ << ",  energy: " << energy_sum << "\n";
    return energy_sum;
  }

  bool SlqFiniteDiscreteControlQuadrotor::feedforwardConverged(std::vector<VectorXd> &u_fw_vec){
    double fw_max = 0.0;
    for (int i = 0; i < iteration_times_; ++i){
      double control_sum = 0.0;
      for (int j = 0; j < u_size_; ++j){
        control_sum += pow((u_fw_vec[i])(j), 2.0);
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
    return stateAddition(relative_x_ptr, xn_ptr_);
  }

  VectorXd SlqFiniteDiscreteControlQuadrotor::getRelativeState(VectorXd *absolute_x_ptr){
    return stateSubtraction(absolute_x_ptr, xn_ptr_);
  }
}


