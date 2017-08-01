// -*- mode: c++ -*-
/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2015, JSK Lab
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

#include <lqr_control/LqrFD.h>

namespace lqr_finite_discrete
{
  LqrFiniteDiscreteControl::LqrFiniteDiscreteControl(ros::NodeHandle nh, ros::NodeHandle nhp): nh_(nh), nhp_(nhp){
  }

  void LqrFiniteDiscreteControl::initLQR(double freq, double period, MatrixXd *A, MatrixXd *B, MatrixXd *Q, MatrixXd *R, VectorXd *s0){
    control_freq_ = freq;
    if (floor(freq * period) < freq * period){
      end_time_ = (floor(freq * period) + 1.0) / freq;
      iteration_times_ = floor(freq * period) + 1;
    }
    else{
      end_time_ = period;
      iteration_times_ = floor(freq * period);
    }
    s_size_ = Q_ptr_->rows();
    u_size_ = R_ptr_->rows();
    A_ptr_ = new MatrixXd(s_size_, s_size_);
    B_ptr_ = new MatrixXd(u_size_, u_size_);
    Q_ptr_ = new MatrixXd(s_size_, s_size_);
    R_ptr_ = new MatrixXd(u_size_, u_size_);
    for (int i = 0; i < s_size_; ++i)
      for (int j = 0; j < s_size_; ++j){
        (*A_ptr_)(i, j) = (*A)(i, j);
        (*Q_ptr_)(i, j) = (*Q)(i, j);
      }
    for (int i = 0; i < u_size_; ++i)
      for (int j = 0; j < u_size_; ++j){
        (*B_ptr_)(i, j) = (*B)(i, j);
        (*R_ptr_)(i, j) = (*R)(i, j);
      }
    for (int i = 0; i < s_size_; ++i)
      (*s0_ptr_)[i] = (*s0)[i];
  }

  void LqrFiniteDiscreteControl::backwardIteration(){
    VectorXd P = MatrixXd::Zero(s_size_, s_size_);
    // todo: assume N is zero, namely do not have x^T *N*u in cost function
    VectorXd N = MatrixXd::Zero(s_size_, u_size_);
    P = *Q_ptr_;
    for (int i = 0; i < iteration_times_; ++i){
      MatrixXd *F_ptr = new MatrixXd(u_size_, s_size_);
      MatrixXd Fk_1 = ((*R_ptr_) + B_ptr_->transpose() * P * (*B_ptr_)).inverse();
      MatrixXd Fk_2 = B_ptr_->transpose() * P * (*A_ptr_) + N.transpose();
      (*F_ptr) = Fk_1 * Fk_2;
      P = A_ptr_->transpose() * P * (*A_ptr_)
        - (A_ptr_->transpose() * P * (*B_ptr_) + N.transpose()) * (*F_ptr)
        + (*Q_ptr_);
      F_ptr_vec_.push_back(F_ptr);
    }
  }

  void LqrFiniteDiscreteControl::forwardUpdateControlValue(){
    x_ptr_vec_.push_back(s0_ptr_);
    for (int i = F_ptr_vec_.size() - 1; i >= 0; --i){
      VectorXd *x_ptr = new VectorXd(u_size_);
      VectorXd *u_ptr = new VectorXd(u_size_);
      (*u_ptr) = -(*(F_ptr_vec_[i])) * (*(x_ptr_vec_[F_ptr_vec_.size() - 1 - i]));
      u_ptr_vec_.push_back(u_ptr);
      (*x_ptr) = (*A_ptr_) * (*(x_ptr_vec_[F_ptr_vec_.size() - 1 - i]))
        +(*B_ptr_) * (*u_ptr);
      x_ptr_vec_.push_back(x_ptr);
    }
  }
}
