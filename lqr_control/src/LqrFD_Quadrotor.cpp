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
namespace lqr_finite_discrete{
  void LqrFiniteDiscreteControlQuadrotor::initLQR(double freq, double period, VectorXd *x0){
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
    B_ptr_ = new MatrixXd(u_size_, u_size_);
    Q_ptr_ = new MatrixXd(x_size_, x_size_);
    R_ptr_ = new MatrixXd(u_size_, u_size_);
    x0_ptr_ = new VectorXd(x_size_);
    u_ptr_ = new VectorXd(u_size_);

    for (int i = 0; i < x_size_; ++i)
      x0_ptr_[i] = x0[i];
  }

  void LqrFiniteDiscreteControlQuadrotor::updateMatrixA(){
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < x_size_; ++j){
        (*A_ptr_)(i, j) = 0.0;
        (*Q_ptr_)(i, j) = 0.0;
      }
    for (int i = 0; i < u_size_; ++i)
      for (int j = 0; j < u_size_; ++j){
        (*B_ptr_)(i, j) = 0.0;
        (*R_ptr_)(i, j) = 0.0;
      }
    /* x, y, z */
    (*A_ptr_)(P_X, V_X) = 1;
    (*A_ptr_)(P_Y, V_Y) = 1;
    (*A_ptr_)(P_Z, V_Z) = 1;

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += (*u_ptr_)[i];
    /* d v_x = (2 * q_w * q_y + 2 * q_x * q_z) * (thrust / m) */
    (*A_ptr_)(V_X, Q_W) = 2 * (*x0_ptr_)[Q_Y] * u;
    (*A_ptr_)(V_X, Q_X) = 2 * (*x0_ptr_)[Q_Z] * u;
    (*A_ptr_)(V_X, Q_Y) = 2 * (*x0_ptr_)[Q_W] * u;
    (*A_ptr_)(V_X, Q_Z) = 2 * (*x0_ptr_)[Q_X] * u;
    /* d v_y = (2 * q_w * q_x + 2 * q_y * q_z) * (thrust / m) */
    (*A_ptr_)(V_Y, Q_W) = 2 * (*x0_ptr_)[Q_X] * u;
    (*A_ptr_)(V_Y, Q_X) = 2 * (*x0_ptr_)[Q_W] * u;
    (*A_ptr_)(V_Y, Q_Y) = 2 * (*x0_ptr_)[Q_Z] * u;
    (*A_ptr_)(V_Y, Q_Z) = 2 * (*x0_ptr_)[Q_Y] * u;
    /* d v_z = (1 - 2 * q_x * q_x - 2 * q_y * q_y) * (thrust / m) */
    (*A_ptr_)(V_Z, Q_X) = -2 * (*x0_ptr_)[Q_X] * u;
    (*A_ptr_)(V_Z, Q_Y) = -2 * (*x0_ptr_)[Q_Y] * u;

    /* q_w, q_x, q_y, q_z */
    /* d q = 1/2 * q * [0, w]^T */
    /* d q_w = q_w * 0 - q_x * w_x - q_y * w_y - q_z * w_z */
    (*A_ptr_)(Q_W, Q_X) = - (*x0_ptr_)[W_X];
    (*A_ptr_)(Q_W, Q_Y) = - (*x0_ptr_)[W_Y];
    (*A_ptr_)(Q_W, Q_Z) = - (*x0_ptr_)[W_Z];
    /* d q_x = q_w * w_x + q_x * 0 + q_y * w_z - q_z * w_y */
    (*A_ptr_)(Q_X, Q_W) = (*x0_ptr_)[W_X];
    (*A_ptr_)(Q_X, Q_Y) = (*x0_ptr_)[W_Z];
    (*A_ptr_)(Q_X, Q_Z) = - (*x0_ptr_)[W_Y];
    /* d q_y = q_w * w_y - q_x * w_z + q_y * 0 + q_z * w_x */
    (*A_ptr_)(Q_Y, Q_W) = (*x0_ptr_)[W_Y];
    (*A_ptr_)(Q_Y, Q_X) = - (*x0_ptr_)[W_Z];
    (*A_ptr_)(Q_Y, Q_Z) = (*x0_ptr_)[W_X];
    /* d q_z = q_w * w_z + q_x * w_y - q_y * w_x + q_z * 0 */
    (*A_ptr_)(Q_Z, Q_W) = (*x0_ptr_)[W_Z];
    (*A_ptr_)(Q_Z, Q_X) = (*x0_ptr_)[W_Y];
    (*A_ptr_)(Q_Z, Q_Y) = (*x0_ptr_)[W_X];

    /* w_x, w_y, w_z */
    /* d w = J^-1 * (- (w^) * (Jw) + tau), w^ = [0, -w_z, w_y; w_z, 0, -w_x; -w_y, w_x, 0] */
    /* d w = J^-1 * (- d(w^) * (Jw) - (w^) * (J * d(w))) */
  }

}


