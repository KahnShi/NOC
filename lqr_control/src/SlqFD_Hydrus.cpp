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

#include <lqr_control/SlqFD_Hydrus.h>
namespace lqr_discrete{
  void SlqFiniteDiscreteControlHydrus::initSLQ(double freq, std::vector<double> *time_ptr, std::vector<VectorXd> *waypoints_ptr){
    /* Ros service */
    dare_client_ = nh_.serviceClient<lqr_control::Dare>("dare_solver");

    control_freq_ = freq;
    time_ptr_ = time_ptr;
    double period = (*time_ptr)[time_ptr->size() - 1] - (*time_ptr)[0];
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

    waypoints_ptr_ = waypoints_ptr;

    /* hydrus */
    link_length_ = 0.44;
    n_links_ = 4;
    for (int i = 0; i < n_links_; ++i)
      link_weight_vec_.push_back(0.55);
    hydrus_weight_ = 0.0;
    for (int i = 0; i < n_links_; ++i)
      hydrus_weight_ += link_weight_vec_[i];
    joint_ptr_ = new VectorXd(n_links_ - 1);
    I_ptr_ = new MatrixXd(3, 3);
    *I_ptr_ = MatrixXd::Zero(3, 3);
    (*I_ptr_)(0, 0) = 0.0001;
    (*I_ptr_)(1, 1) = 0.0001;
    (*I_ptr_)(2, 2) = 0.0002;

    x_size_ = 12;
    u_size_ = 4;
    A_ptr_ = new MatrixXd(x_size_, x_size_);
    B_ptr_ = new MatrixXd(x_size_, u_size_);
    Q0_ptr_ = new MatrixXd(x_size_, x_size_);
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

    *x0_ptr_ = (*waypoints_ptr)[0];
    *xn_ptr_ = (*waypoints_ptr)[waypoints_ptr->size() - 1];

    /* init Q and R matrice */
    *Q0_ptr_ = MatrixXd::Zero(x_size_, x_size_);
    for (int i = 0; i <= P_Z; ++i)
      (*Q0_ptr_)(i, i) = 10.0;
    for (int i = V_X; i <= V_Z; ++i)
      (*Q0_ptr_)(i, i) = 10.0;
    for (int i = W_X; i < x_size_; ++i)
      (*Q0_ptr_)(i, i) = 1.0;
    // test: weight on z
    (*Q0_ptr_)(P_Z, P_Z) = (*Q0_ptr_)(V_Z, V_Z) = 100.0;
    *Q_ptr_ = (*Q0_ptr_);

    *R_ptr_ = 400 * MatrixXd::Identity(u_size_, u_size_);

    *M_para_ptr_ = MatrixXd::Zero(3, 4);
    double r_uav = 0.3, c_rf = 0.016;
    (*M_para_ptr_)(0, 1) = -r_uav;
    (*M_para_ptr_)(0, 3) = r_uav;
    (*M_para_ptr_)(1, 0) = r_uav;
    (*M_para_ptr_)(1, 2) = -r_uav;
    (*M_para_ptr_)(2, 0) = (*M_para_ptr_)(2, 2) = -c_rf;
    (*M_para_ptr_)(2, 1) = (*M_para_ptr_)(2, 3) = c_rf;

    uav_rotor_thrust_min_ = 0.0;
    uav_rotor_thrust_max_ = (hydrus_weight_ * 9.78 / 4.0) * 3.0;

    /* Assume initial and final state is still, namely dx = [v, a] = 0 */
    u0_ptr_ = new VectorXd(u_size_);
    un_ptr_ = new VectorXd(u_size_);
    for (int i = 0; i < 4; ++i){
      (*u0_ptr_)(i) = hydrus_weight_ * 9.78 / 4.0;
      (*un_ptr_)(i) = hydrus_weight_ * 9.78 / 4.0;
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
      VectorXd cur_joint = getCurrentJoint(i/control_freq_);
      joint_vec_.push_back(cur_joint);
      getHydrusLinksCenter(&cur_joint);
      getHydrusInertialTensor(&cur_joint, i);
      u_fw_vec_.push_back(u_init);
      u_fb_vec_.push_back(u_init);
      K_vec_.push_back(MatrixXd::Zero(u_size_, x_size_));
    }

    line_search_steps_ = 4;

    debug_ = true;
    FDLQR();
    getRiccatiH();
    ROS_INFO("[SLQ] Initialization finished.");
  }

  void SlqFiniteDiscreteControlHydrus::getRiccatiH(){
    // method 1: use real state at time tf (get from initial LQR result)
    *x_ptr_ = x_vec_[iteration_times_];
    *u_ptr_ = u_vec_[iteration_times_];
    // method 2: use ideal state at time tf
    //*x_ptr_ = VectorXd::Zero(x_size_);
    //*u_ptr_ = VectorXd::Zero(u_size_);
    if (debug_){
      printStateInfo(x_ptr_, iteration_times_);
      printControlInfo(u_ptr_, iteration_times_);
    }
    *joint_ptr_ = joint_vec_[iteration_times_];
    updateMatrixAB(iteration_times_);

    /* debug: print matrix A and B */
    // if (debug_){
    //   printMatrixAB();
    // }

    lqr_control::float64Array dat_A, dat_B, dat_Q, dat_R, dat_P;

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
    dat_Q.array.layout.dim.push_back(std_msgs::MultiArrayDimension());
    dat_Q.array.layout.dim.push_back(std_msgs::MultiArrayDimension());
    dat_Q.array.layout.dim[0].label = "height";
    dat_Q.array.layout.dim[1].label = "width";
    dat_Q.array.layout.dim[0].size = x_size_;
    dat_Q.array.layout.dim[1].size = x_size_;
    dat_Q.array.layout.dim[0].stride = x_size_ * x_size_;
    dat_Q.array.layout.dim[1].stride = x_size_;
    dat_Q.array.layout.data_offset = 0;
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < x_size_; ++j)
        dat_Q.array.data.push_back((*Q_ptr_)(i, j));
    dat_R.array.layout.dim.push_back(std_msgs::MultiArrayDimension());
    dat_R.array.layout.dim.push_back(std_msgs::MultiArrayDimension());
    dat_R.array.layout.dim[0].label = "height";
    dat_R.array.layout.dim[1].label = "width";
    dat_R.array.layout.dim[0].size = u_size_;
    dat_R.array.layout.dim[1].size = u_size_;
    dat_R.array.layout.dim[0].stride = u_size_ * u_size_;
    dat_R.array.layout.dim[1].stride = u_size_;
    dat_R.array.layout.data_offset = 0;
    for (int i = 0; i < u_size_; ++i)
      for (int j = 0; j < u_size_; ++j)
        dat_R.array.data.push_back((*R_ptr_)(i, j));

    lqr_control::Dare dare_srv;
    dare_srv.request.A = dat_A;
    dare_srv.request.B = dat_B;
    dare_srv.request.Q = dat_Q;
    dare_srv.request.R = dat_R;
    ROS_INFO("[SLQ] Matrix AB is sent to Riccati solver.");
    if (dare_client_.call(dare_srv))
      dat_P = dare_srv.response.P;
    else
      ROS_ERROR("[SLQ] No response from dare sever.");
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < x_size_; ++j)
        (*Riccati_P_ptr_)(i, j) = dat_P.array.data[i * x_size_ + j];
    ROS_INFO("[SLQ] Matrix P is received from Riccati solver.");

    /* debug: output matrix result from riccati eqation */
    // if (debug_){
    //   std::cout << "[Debug] print matrix P initial value:\n";
    //   for (int i = 0; i < x_size_; ++i){
    //     for (int j = 0; j < x_size_; ++j)
    //       std::cout << (*Riccati_P_ptr_)(i, j) << ", ";
    //     std::cout << "\n";
    //   }
    // }

    ROS_INFO("[SLQ] Get P matrice initial value from Ricatti function.");
  }

  void SlqFiniteDiscreteControlHydrus::iterativeOptimization(){
    // *p_ptr_ = VectorXd::Zero(x_size_);
    // *r_ptr_ = (*R_ptr_) * (*un_ptr_);
    *P_ptr_ = *Riccati_P_ptr_;
    *p_ptr_ = (*P_ptr_) * (x_vec_[iteration_times_]);
    // test: real u
    // *r_ptr_ = VectorXd::Zero(u_size_);

    for (int i = iteration_times_ - 1; i >= 0; --i){
      // add weight for waypoints
      std::vector<MatrixXd> W_vec;
      for (int j = 1; j < waypoints_ptr_->size(); ++j){
        MatrixXd W = MatrixXd::Zero(x_size_, x_size_);
        updateWaypointWeightMatrix(i * end_time_ / iteration_times_, (*time_ptr_)[j] - (*time_ptr_)[0], &W, j == (waypoints_ptr_->size() - 1));
        W_vec.push_back(W);
      }

      // update current Q and R matrix
      *Q_ptr_ = (*Q0_ptr_);
      for (int j = 1; j < waypoints_ptr_->size(); ++j)
        *Q_ptr_ = *Q_ptr_ + W_vec[j-1];

      *x_ptr_ = x_vec_[i];
      *u_ptr_ = u_vec_[i];
      *joint_ptr_ = joint_vec_[i];
      updateMatrixAB(i);

      *q_ptr_ = (*Q0_ptr_) * x_vec_[i];
      for (int j = 1; j < waypoints_ptr_->size(); ++j)
        *q_ptr_ = (*q_ptr_) +
          W_vec[j-1] * (*xn_ptr_ - (*waypoints_ptr_)[j] + x_vec_[i]);

      *r_ptr_ = (*R_ptr_) * (u_vec_[i]);
      // *r_ptr_ = (*R_ptr_) * ((u_vec_[i]) + (*un_ptr_));
      updateSLQEquations();

      Vector4d u_fb = (*K_ptr_) * (*x_ptr_);
      u_fb_vec_[i] = u_fb;
      u_fw_vec_[i] = (*l_ptr_);
      K_vec_[i] = (*K_ptr_);
    }

    /* Update control by finding the best alpha */
    alpha_ = 1.0;
    double alpha_candidate = 1.0, energy_min = -1.0;
    double search_rate = 2.0;
    if (feedforwardConverged() && debug_)
      std::cout << "[SLQ] feedforward converge.";
    else{
      for (int factor = 0; factor < line_search_steps_; ++factor){
        double energy_sum = 0.0;
        VectorXd cur_u(u_size_);
        VectorXd cur_x = x_vec_[0];
        for (int i = 0; i < iteration_times_; ++i){
          VectorXd cur_u(u_size_);
          cur_u = u_vec_[i] + alpha_ * u_fw_vec_[i]
            + K_vec_[i] * (cur_x - x_vec_[i]);
          checkControlInputFeasible(&cur_u);
          // calculate energy
          // add weight for waypoints
          std::vector<MatrixXd> W_vec;
          for (int j = 1; j < waypoints_ptr_->size(); ++j){
            MatrixXd W = MatrixXd::Zero(x_size_, x_size_);
            updateWaypointWeightMatrix(i * end_time_ / iteration_times_, (*time_ptr_)[j] - (*time_ptr_)[0], &W, j == (waypoints_ptr_->size() - 1));
            W_vec.push_back(W);
          }

          VectorXd real_u = cur_u + (*un_ptr_);
          // method 1: use "relative" u when calculting whole energy
          // energy_sum += (cur_x.transpose() * (*Q_ptr_) * cur_x
          //                + cur_u.transpose() * (*R_ptr_) * cur_u)(0);
          // method 2: use real u when calculting whole energy
          energy_sum += (real_u.transpose() * (*R_ptr_) * real_u)(0);

          energy_sum += (cur_x.transpose() * (*Q0_ptr_) * cur_x)(0);
          for (int j = 1; j < waypoints_ptr_->size(); ++j){
            VectorXd dx_pt = cur_x + (*xn_ptr_) - (*waypoints_ptr_)[j];
            energy_sum += (dx_pt.transpose() * W_vec[j-1] * dx_pt)(0);
          }

          VectorXd new_x(x_size_);
          VectorXd cur_joint = joint_vec_[i];
          updateNewState(&new_x, &cur_x, &cur_u, &cur_joint);
          cur_x = new_x;
        }
        energy_sum += (cur_x.transpose() * (*Riccati_P_ptr_) * cur_x)(0);

        // energy and alpha' relationships
        // std::cout << "[SLQ] Energy: " << energy_sum << ", alpha: " << alpha_ << "\n";

        if (energy_sum < energy_min || energy_min < 0){
          energy_min = energy_sum;
          alpha_candidate = alpha_;
        }
        alpha_ = alpha_ / search_rate;
      }
    }
    if (debug_)
      std::cout << "\nAlpha selected: " << alpha_candidate << "\n\n";

    // test K with every new state
    VectorXd cur_x(x_size_);
    cur_x = x_vec_[0];
    for (int i = 0; i < iteration_times_; ++i){
      VectorXd cur_u(u_size_);
      cur_u = u_vec_[i] + alpha_candidate * u_fw_vec_[i]
        + K_vec_[i] * (cur_x - x_vec_[i]);
      checkControlInputFeasible(&cur_u);
      VectorXd new_x(x_size_);
      VectorXd cur_joint = joint_vec_[i];
      updateNewState(&new_x, &cur_x, &cur_u, &cur_joint);
      if ((i % 100 == 0 || i == iteration_times_ - 1) && debug_){
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

  VectorXd SlqFiniteDiscreteControlHydrus::getCurrentJoint(double time, int order){
    VectorXd joint = VectorXd::Zero(n_links_ - 1);
    // keep quadrotor model, neglecting time, order
    for (int i = 0; i < n_links_ - 1; ++i)
      joint(i) = 3.14159 / 2.0;
    return joint;
  }

  void SlqFiniteDiscreteControlHydrus::updateMatrixAB(int time_id){
    updateMatrixA(time_id);
    updateMatrixB(time_id);
  }

  void SlqFiniteDiscreteControlHydrus::updateMatrixA(int time_id){
    *A_ptr_ = MatrixXd::Zero(x_size_, x_size_);

    VectorXd *x_ptr = new VectorXd(x_size_); *x_ptr = x_vec_[time_id];
    VectorXd *u_ptr = new VectorXd(u_size_); *u_ptr = u_vec_[time_id];
    VectorXd *joint_ptr = new VectorXd(n_links_ - 1); *joint_ptr = joint_vec_[time_id];

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
    u = u / hydrus_weight_;
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

  void SlqFiniteDiscreteControlHydrus::updateMatrixB(int time_id){
    *B_ptr_ = MatrixXd::Zero(x_size_, u_size_);

    VectorXd *x_ptr = new VectorXd(x_size_); *x_ptr = x_vec_[time_id];
    VectorXd *u_ptr = new VectorXd(u_size_); *u_ptr = u_vec_[time_id];
    VectorXd *joint_ptr = new VectorXd(n_links_ - 1); *joint_ptr = joint_vec_[time_id];

    /* x, y, z */
    /* all 0 */

    /* v_x, v_y, v_z */
    /* d v_x = (sin y * sin r + cos y * sin p * cos r) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_X, U_1) = (sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) +
                           cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      / hydrus_weight_;
    /* d v_y = (-cos y * sin r + sin y * sin p * cos r) * (u1 + u2 + u3 + u4) / m  */
    (*B_ptr_)(V_Y, U_1) = (-cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) +
                           sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      / hydrus_weight_;
    /* d v_z = (cos p * cos r) * (u1 + u2 + u3 + u4) / m */
    (*B_ptr_)(V_Z, U_1) = (cos((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      / hydrus_weight_;
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

  void SlqFiniteDiscreteControlHydrus::updateNewState(VectorXd *new_x_ptr, VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr){
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
      * u / hydrus_weight_;
    /* d v_y = (-cos y * sin r + sin y * sin p * cos r) * (u1 + u2 + u3 + u4) / m  */
    dev_x(V_Y) = (-cos((*x_ptr)[E_Y]) * sin((*x_ptr)[E_R]) +
                  sin((*x_ptr)[E_Y]) * sin((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      * u / hydrus_weight_;
    /* d v_z = (cos p * cos r) * (u1 + u2 + u3 + u4) / m */
    dev_x(V_Z) = (cos((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]))
      * u / hydrus_weight_ - 9.78;

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

  bool SlqFiniteDiscreteControlHydrus::feedforwardConverged(){
    double fw_max = 0.0;
    for (int i = 0; i < iteration_times_; ++i){
      double control_sum = 0.0;
      for (int j = 0; j < u_size_; ++j){
        control_sum += pow((u_fw_vec_[i])(j), 2.0);
      }
      if (control_sum > fw_max)
        fw_max = control_sum;
    }
    double fw_converge_threshold = (hydrus_weight_ * 9.78 / 4.0) * 0.1;
    if (fw_max < pow(fw_converge_threshold, 2))
      return true;
    else
      return false;
  }

  void SlqFiniteDiscreteControlHydrus::getHydrusInertialTensor(VectorXd *joint_ptr, int time_id){
    std::vector<MatrixXd> cur_I_vec;
    for (int i = 0; i < n_links_; ++i){
      MatrixXd cur_I = *I_ptr_;
      Vector3d center_pos = link_center_pos_local_vec_[time_id][i];
      cur_I(0, 0) += link_weight_vec_[i] * (pow(center_pos(1)/2.0, 2.0)
                                            + pow(center_pos(2)/2.0, 2.0));
      cur_I(1, 1) += link_weight_vec_[i] * (pow(center_pos(0)/2.0, 2.0)
                                            + pow(center_pos(2)/2.0, 2.0));
      cur_I(2, 2) += link_weight_vec_[i] * (pow(center_pos(0)/2.0, 2.0)
                                            + pow(center_pos(1)/2.0, 2.0));
      cur_I(0, 1) = -link_weight_vec_[i] * center_pos(0) * center_pos(1) / 4.0;
      cur_I(0, 2) = -link_weight_vec_[i] * center_pos(0) * center_pos(2) / 4.0;
      cur_I(1, 2) = -link_weight_vec_[i] * center_pos(1) * center_pos(2) / 4.0;
      cur_I(1, 0) = cur_I(0, 1);
      cur_I(2, 0) = cur_I(0, 2);
      cur_I(2, 1) = cur_I(1, 2);

      cur_I_vec.push_back(cur_I);
    }
    I_vec_.push_back(cur_I_vec);
  }

  void SlqFiniteDiscreteControlHydrus::getHydrusLinksCenter(VectorXd *joint_ptr){
    std::vector<Vector3d> links_center_vec;
    Vector3d link1_center(link_length_ / 2.0, 0, 0);
    links_center_vec.push_back(link1_center);
    Vector3d prev_link_end(link_length_, 0, 0);
    double joint_ang = 0.0;

    // only considering 2d hydrus
    for (int i = 1; i < n_links_; ++i){
      Vector3d link_center = Vector3d::Zero();
      joint_ang += (*joint_ptr)(i - 1);
      Matrix3d rot;
      rot << cos(joint_ang), -sin(joint_ang), 0,
        sin(joint_ang), cos(joint_ang), 0,
        0, 0, 1;
      link_center = prev_link_end + rot * Vector3d(link_length_ / 2.0, 0, 0);
      links_center_vec.push_back(link_center);
      prev_link_end = prev_link_end + rot * Vector3d(link_length_, 0, 0);
    }
    link_center_pos_local_vec_.push_back(links_center_vec);
  }

  VectorXd SlqFiniteDiscreteControlHydrus::stateAddition(VectorXd *x1_ptr, VectorXd *x2_ptr){
    VectorXd result(x_size_);
    result = (*x1_ptr) + (*x2_ptr);
    // todo: euler angle addition
    return result;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::stateSubtraction(VectorXd *x1_ptr, VectorXd *x2_ptr){
    VectorXd result(x_size_);
    result = (*x1_ptr) - (*x2_ptr);
    // todo: euler angle subtraction
    return result;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getAbsoluteState(VectorXd *relative_x_ptr){
    // todo
    return stateAddition(relative_x_ptr, xn_ptr_);
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getRelativeState(VectorXd *absolute_x_ptr){
    // todo
    return stateSubtraction(absolute_x_ptr, xn_ptr_);
  }

  void SlqFiniteDiscreteControlHydrus::updateWaypointWeightMatrix(double time, double end_time, MatrixXd *W_ptr, bool goal_flag){
    double rho = 1.0;
    double weight = exp(-rho / 2 * pow(time - end_time, 2.0));
    for (int j = 0; j < 6; ++j)
      (*W_ptr)(j, j) = 10.0 * weight;
    for (int j = 6; j < x_size_; ++j)
      (*W_ptr)(j, j) = 0.1 * weight;
    (*W_ptr)(E_Y, E_Y) = weight;
    (*W_ptr)(W_Y, W_Y) = weight;
    // test: weight on z
    (*W_ptr)(2, 2) = (*W_ptr)(5, 5) = 100.0 * weight;

    // test: more weight on mid state
    if (!goal_flag)
      *W_ptr = (*W_ptr) * 1.0;
  }

  void SlqFiniteDiscreteControlHydrus::updateSLQEquations(){
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

  void SlqFiniteDiscreteControlHydrus::FDLQR(){
    *x_ptr_ = x_vec_[0];
    *u_ptr_ = u_vec_[0];
    *joint_ptr_ = joint_vec_[0];
    updateMatrixAB(0);

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
      VectorXd cur_joint = joint_vec_[i];
      updateNewState(&new_x, &x, &u, &cur_joint);
      x = new_x;
      x_vec_[iteration_times_ - i] = x;

      // Guarantee control is in bound
      checkControlInputFeasible(&u);
      u_vec_[iteration_times_ - i] = u;

      if ((i % 100 == 0 || i == iteration_times_ - 1) && debug_){
        printStateInfo(&x, i);
        printControlInfo(&u, i);
      }
    }
    ROS_INFO("[SLQ] LQR init finished");
  }

  void SlqFiniteDiscreteControlHydrus::checkControlInputFeasible(VectorXd *u){
    for (int j = 0; j < u_size_; ++j){
      if ((*u)(j) + (*un_ptr_)(j) < uav_rotor_thrust_min_)
        (*u)(j) = uav_rotor_thrust_min_ - (*un_ptr_)(j);
      else if ((*u)(j) + (*un_ptr_)(j) > uav_rotor_thrust_max_)
        (*u)(j) = uav_rotor_thrust_max_ - (*un_ptr_)(j);
    }
  }

  void SlqFiniteDiscreteControlHydrus::printStateInfo(VectorXd *x, int id){
    VectorXd new_absolute_x;
    new_absolute_x = getAbsoluteState(x);
    std::cout << "[debug] id[" << id << "]print current state:\n";
    for (int j = 0; j < x_size_; ++j)
      std::cout << new_absolute_x(j) << ", ";
    std::cout << "\n";
  }

  void SlqFiniteDiscreteControlHydrus::printControlInfo(VectorXd *u, int id){
    std::cout << "[debug] id[" << id << "]print current u:\n";
    for (int j = 0; j < u_size_; ++j)
      std::cout << (*u)(j) + (*un_ptr_)(j) << ", ";
    std::cout << "\n";
  }

  void SlqFiniteDiscreteControlHydrus::printMatrixAB(){
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


