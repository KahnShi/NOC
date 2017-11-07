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
  void SlqFiniteDiscreteControlHydrus::initHydrus(){
    /* Ros service */
    dare_client_ = nh_.serviceClient<lqr_control::Dare>("/dare_solver");

    not_first_slq_flag_ = false;
    /* ros param */
    double R_para, Q_p_para, Q_v_para, Q_e_para, Q_w_para, Q_z_para;
    nhp_.param("transform_movement_flag", transform_movement_flag_, true);
    nhp_.param("R_para", R_para, 10.0);
    nhp_.param("Q_p_para", Q_p_para, 5.0);
    nhp_.param("Q_v_para", Q_v_para, 10.0);
    nhp_.param("Q_z_para", Q_z_para, 50.0);
    nhp_.param("Q_w_para", Q_w_para, 10.0);
    nhp_.param("Q_e_para", Q_e_para, 100.0);

    debug_ = false;

    /* hydrus */
    link_length_ = 0.6;
    n_links_ = 4;
    for (int i = 0; i < n_links_; ++i)
      link_weight_vec_.push_back(0.92);
    hydrus_weight_ = 0.0;
    for (int i = 0; i < n_links_; ++i)
      hydrus_weight_ += link_weight_vec_[i];
    joint_ptr_ = new VectorXd(n_links_ - 1);
    I_ptr_ = new Eigen::Matrix3d();
    *I_ptr_ = Eigen::Matrix3d::Zero(3, 3);
    // Inertial is too small to ignore.
    // (*I_ptr_)(0, 0) = 0.0001;
    // (*I_ptr_)(1, 1) = 0.0001;
    // (*I_ptr_)(2, 2) = 0.0002;
    double c_rf = -0.01676;
    M_z_ = VectorXd::Zero(n_links_);
    M_z_(0) = -c_rf;
    M_z_(1) = c_rf;
    M_z_(2) = -c_rf;
    M_z_(3) = c_rf;

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
    Riccati_P_ptr_ = new MatrixXd(x_size_, x_size_);
    IDlqr_F_ptr_ = new MatrixXd(u_size_, x_size_);
    P_ptr_ = new MatrixXd(x_size_, x_size_);
    p_ptr_ = new VectorXd(x_size_);
    H_ptr_ = new MatrixXd(u_size_, u_size_);
    G_ptr_ = new MatrixXd(u_size_, x_size_);
    K_ptr_ = new MatrixXd(u_size_, x_size_);
    g_ptr_ = new VectorXd(u_size_);
    l_ptr_ = new VectorXd(u_size_);
    r_ptr_ = new VectorXd(u_size_);
    q_ptr_ = new VectorXd(x_size_);
    waypoints_ptr_ = new std::vector<VectorXd>();
    time_ptr_ = new std::vector<double>();

    /* init Q and R matrice */
    *Q0_ptr_ = MatrixXd::Zero(x_size_, x_size_);
    for (int i = 0; i <= P_Z; ++i)
      (*Q0_ptr_)(i, i) = Q_p_para;
    for (int i = V_X; i <= V_Z; ++i)
      (*Q0_ptr_)(i, i) = Q_v_para;
    for (int i = W_X; i <= W_Z; ++i)
      (*Q0_ptr_)(i, i) = Q_w_para;
    for (int i = E_R; i <= E_Y; ++i)
      (*Q0_ptr_)(i, i) = Q_e_para;
    // test: weight on z
    (*Q0_ptr_)(P_Z, P_Z) = (*Q0_ptr_)(V_Z, V_Z) = Q_z_para;

    *R_ptr_ = R_para * MatrixXd::Identity(u_size_, u_size_);

    uav_rotor_thrust_min_ = 0.0;
    uav_rotor_thrust_max_ = 18.0;

    ROS_INFO("[SLQ] Hydrus init finished.");
  }

  void SlqFiniteDiscreteControlHydrus::initSLQ(double freq, std::vector<double> *time_ptr, std::vector<VectorXd> *waypoints_ptr){
    control_freq_ = freq;
    // Here we assume frequency is an odd integer
    control_high_freq_ = freq;
    control_low_freq_ = freq / 2.0;
    double period = (*time_ptr)[time_ptr->size() - 1] - (*time_ptr)[0];
    double high_freq_default_period = 1.0;
    if (period <= high_freq_default_period){
      double high_freq_period = period;
      if (floor(control_high_freq_ * high_freq_period) < control_high_freq_ * high_freq_period){
        end_time_ = (floor(control_high_freq_ * high_freq_period) + 1.0) / control_high_freq_;
        iteration_times_ = floor(control_high_freq_ * high_freq_period) + 1;
      }
      else{
        end_time_ = high_freq_period;
        iteration_times_ = floor(control_high_freq_ * high_freq_period);
      }
      high_freq_end_time_ = end_time_;
      low_freq_end_time_ = end_time_;
      high_freq_iteration_times_ = iteration_times_;
      low_freq_iteration_times_ = 0.0;
    }
    else{
      double low_freq_period = period - high_freq_default_period;
      high_freq_end_time_ = high_freq_default_period;
      high_freq_iteration_times_ = floor(control_high_freq_ * high_freq_end_time_);
      if (floor(control_low_freq_ * low_freq_period) < control_low_freq_ * low_freq_period){
        low_freq_iteration_times_ = floor(control_low_freq_ * low_freq_period) + 1;
        end_time_ = low_freq_iteration_times_ / control_low_freq_ + high_freq_end_time_;
        low_freq_end_time_ = end_time_;
        iteration_times_ = low_freq_iteration_times_ + high_freq_iteration_times_;
      }
      else{
        low_freq_iteration_times_ = floor(control_low_freq_ * low_freq_period);
        end_time_ = period;
        low_freq_end_time_ = end_time_;
        iteration_times_ = low_freq_iteration_times_ + high_freq_iteration_times_;
      }
    }

    if (debug_){
      std::cout << "[SLQ] Trajectory period: " << end_time_
                << ", Itetation times: " << iteration_times_ << "\n";
      std::cout << "[SLQ] Start position: " << (*waypoints_ptr)[0].transpose() << "\n";
      std::cout << "[SLQ] End position: " << (*waypoints_ptr)[waypoints_ptr->size()  - 1].transpose() << "\n";
      std::cout << "[SLQ] High frequency end time: " << high_freq_end_time_ << ", times: " << high_freq_iteration_times_ << "\n";
      std::cout << "[SLQ] Low frequency end time: " << low_freq_end_time_ << ", times: " << low_freq_iteration_times_ << "\n";
    }

    if (waypoints_ptr_->size())
      waypoints_ptr_->clear();
    for (int i = 0; i < waypoints_ptr->size(); ++i)
      waypoints_ptr_->push_back((*waypoints_ptr)[i]);
    if (time_ptr_->size())
      time_ptr_->clear();
    for (int i = 0; i < time_ptr->size(); ++i)
      time_ptr_->push_back((*time_ptr)[i]);

    *A_ptr_ = MatrixXd::Zero(x_size_, x_size_);
    *B_ptr_ = MatrixXd::Zero(x_size_, u_size_);
    *x_ptr_ = VectorXd::Zero(x_size_);
    *u_ptr_ = VectorXd::Zero(u_size_);
    *Riccati_P_ptr_ = MatrixXd::Zero(x_size_, x_size_);
    *P_ptr_ = MatrixXd::Zero(x_size_, x_size_);
    *p_ptr_ = VectorXd::Zero(x_size_);
    *H_ptr_ = MatrixXd::Zero(u_size_, u_size_);
    *G_ptr_ = MatrixXd::Zero(u_size_, x_size_);
    *K_ptr_ = MatrixXd::Zero(u_size_, x_size_);
    *g_ptr_ = VectorXd::Zero(u_size_);
    *l_ptr_ = VectorXd::Zero(u_size_);
    *r_ptr_ = VectorXd::Zero(u_size_);
    *q_ptr_ = VectorXd::Zero(x_size_);

    *x0_ptr_ = (*waypoints_ptr)[0];
    *xn_ptr_ = (*waypoints_ptr)[waypoints_ptr->size() - 1];
    *Q_ptr_ = (*Q0_ptr_);

    /* SLQ special initialization */
    // todo: assume start point the quadrotor is hovering
    VectorXd x_init(x_size_), u_init(u_size_);
    x_init = getRelativeState(x0_ptr_);
    u_init = VectorXd::Zero(u_size_);

    /* Clear assigned vector */
    if (not_first_slq_flag_){
      x_vec_.clear();
      u_vec_.clear();
      joint_vec_.clear();
      joint_dt_vec_.clear();
      joint_ddt_vec_.clear();
      I_vec_.clear();
      I_dt_vec_.clear();
      link_center_pos_local_vec_.clear();
      link_center_pos_local_dt_vec_.clear();
      cog_pos_local_vec_.clear();
      cog_pos_local_dt_vec_.clear();
      u_fw_vec_.clear();
      u_fb_vec_.clear();
      K_vec_.clear();
      un_vec_.clear();
    }
    else
      not_first_slq_flag_ = true;

    for (int i = 0; i <= iteration_times_; ++i){
      x_vec_.push_back(x_init);
      u_vec_.push_back(u_init);
      double cur_time;
      if (i <= high_freq_iteration_times_)
        cur_time = double(i) / control_high_freq_;
      else
        cur_time = high_freq_end_time_ +
          (i - high_freq_iteration_times_) / control_low_freq_;
      VectorXd cur_joint = getCurrentJoint(cur_time);
      VectorXd cur_joint_dt = getCurrentJoint(cur_time, 1);
      VectorXd cur_joint_ddt = getCurrentJoint(cur_time, 2);
      joint_vec_.push_back(cur_joint);
      joint_dt_vec_.push_back(cur_joint_dt);
      joint_ddt_vec_.push_back(cur_joint_ddt);
      getHydrusLinksCenter(&cur_joint);
      getHydrusLinksCenterDerivative(&cur_joint, &cur_joint_dt);
      updateHydrusCogPosition(i);
      updateHydrusCogPositionDerivative(i);
      getHydrusInertialTensor(&cur_joint, i);
      u_fw_vec_.push_back(u_init);
      u_fb_vec_.push_back(u_init);
      K_vec_.push_back(MatrixXd::Zero(u_size_, x_size_));
      VectorXd stable_u = getStableThrust(i);
      un_vec_.push_back(stable_u);
    }
    stable_u_last_ = un_vec_[iteration_times_];

    line_search_steps_ = 4;

    FDLQR();
    getRiccatiH();
    *IDlqr_F_ptr_ = (*R_ptr_ +
                     B_ptr_->transpose() * (*Riccati_P_ptr_) * (*B_ptr_)).inverse()
      * (B_ptr_->transpose() * (*Riccati_P_ptr_) * (*A_ptr_));

    if (debug_){
      double lqr_cost = calculateCostFunction();
      std::cout << "\n Cost: " << lqr_cost << "\n\n";
    }
    ROS_INFO("[SLQ] Initialization finished.");
  }

  void SlqFiniteDiscreteControlHydrus::getRiccatiH(){
    // method 1: use real state at time tf (get from initial LQR result)
    // *x_ptr_ = x_vec_[iteration_times_];
    // *u_ptr_ = u_vec_[iteration_times_];
    // method 2: use ideal state at time tf
    *x_ptr_ = VectorXd::Zero(x_size_);
    *u_ptr_ = VectorXd::Zero(u_size_);
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
    *P_ptr_ = *Riccati_P_ptr_;
    *p_ptr_ = VectorXd::Zero(x_size_);

    for (int i = iteration_times_ - 1; i >= 0; --i){
      // add weight for waypoints
      std::vector<MatrixXd> W_vec;
      for (int j = 1; j < waypoints_ptr_->size(); ++j){
        MatrixXd W = MatrixXd::Zero(x_size_, x_size_);
        double cur_time;
        if (i <= high_freq_iteration_times_)
          cur_time = double(i) / control_high_freq_;
        else
          cur_time = high_freq_end_time_ +
            (i - high_freq_iteration_times_) / control_low_freq_;
        updateWaypointWeightMatrix(cur_time, (*time_ptr_)[j] - (*time_ptr_)[0], &W, j == (waypoints_ptr_->size() - 1));
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

      *q_ptr_ = VectorXd::Zero(x_size_);
      for (int j = 1; j < waypoints_ptr_->size(); ++j)
        *q_ptr_ = (*q_ptr_) +
          2.0 * W_vec[j-1] * stateSubtraction(xn_ptr_, &((*waypoints_ptr_)[j]));

      *r_ptr_ = VectorXd::Zero(u_size_);
      updateSLQEquations();

      Vector4d u_fb = (*K_ptr_) * (*x_ptr_);
      u_fb_vec_[i] = u_fb;
      u_fw_vec_[i] = (*l_ptr_);
      K_vec_[i] = (*K_ptr_);
    }

    /* Update control by finding the best alpha */
    alpha_ = 1.0;
    alpha_candidate_ = 1.0;
    double energy_min = -1.0;
    double search_rate = 2.0;
    /* When there are no middle waypoints, feedforward term is 0. */
    bool alpha_iteration_flag = true;
    if (waypoints_ptr_->size() == 2 || feedforwardConverged()){
      alpha_iteration_flag = false;
      if (debug_)
        std::cout << "[SLQ] feedforward converge.";
    }
    else{
      for (int factor = 0; factor < line_search_steps_; ++factor){
        double energy_sum = 0.0;
        VectorXd cur_u(u_size_);
        VectorXd cur_x = x_vec_[0];
        for (int i = 0; i < iteration_times_; ++i){
          VectorXd cur_u(u_size_);
          cur_u = u_vec_[i] + alpha_ * u_fw_vec_[i]
            + K_vec_[i] * cur_x;
          checkControlInputFeasible(&cur_u, i);
          // calculate energy
          // add weight for waypoints
          std::vector<MatrixXd> W_vec;
          for (int j = 1; j < waypoints_ptr_->size(); ++j){
            MatrixXd W = MatrixXd::Zero(x_size_, x_size_);
            double cur_time;
            if (i <= high_freq_iteration_times_)
              cur_time = double(i) / control_high_freq_;
            else
              cur_time = high_freq_end_time_ +
                (i - high_freq_iteration_times_) / control_low_freq_;
            updateWaypointWeightMatrix(cur_time, (*time_ptr_)[j] - (*time_ptr_)[0], &W, j == (waypoints_ptr_->size() - 1));
            W_vec.push_back(W);
          }

          energy_sum += (cur_u.transpose() * (*R_ptr_) * cur_u)(0);
          energy_sum += (cur_x.transpose() * (*Q0_ptr_) * cur_x)(0);
          VectorXd real_x = getAbsoluteState(&cur_x);
          for (int j = 1; j < waypoints_ptr_->size(); ++j){
            VectorXd dx_pt = stateSubtraction(&real_x, &((*waypoints_ptr_)[j]));
            energy_sum += (dx_pt.transpose() * W_vec[j-1] * dx_pt)(0);
          }

          VectorXd new_x(x_size_);
          updateNewState(&new_x, &cur_x, &cur_u, i);
          cur_x = new_x;
        }
        energy_sum += (cur_x.transpose() * (*Riccati_P_ptr_) * cur_x)(0);

        // energy and alpha' relationships
        // std::cout << "[SLQ] Energy: " << energy_sum << ", alpha: " << alpha_ << "\n";

        if (energy_sum < energy_min || energy_min < 0){
          energy_min = energy_sum;
          alpha_candidate_ = alpha_;
        }
        alpha_ = alpha_ / search_rate;
      }
    }
    if (debug_ && alpha_iteration_flag){
      std::cout << "\nMininum energy: " << energy_min << "\n";
      std::cout << "Alpha selected: " << alpha_candidate_ << "\n\n";
    }

    // update new state
    VectorXd cur_x(x_size_);
    cur_x = x_vec_[0];
    for (int i = 0; i < iteration_times_; ++i){
      VectorXd cur_u(u_size_);
      cur_u = u_vec_[i] + alpha_candidate_ * u_fw_vec_[i]
        + K_vec_[i] * cur_x;
      checkControlInputFeasible(&cur_u, i);
      VectorXd new_x(x_size_);
      updateNewState(&new_x, &cur_x, &cur_u, i);
      if ((i % int(control_freq_) == 0 || i == iteration_times_ - 1) && debug_){
        printStateInfo(&cur_x, i);
        printControlInfo(&cur_u, i);
      }
      x_vec_[i] = cur_x;
      u_vec_[i] = cur_u;
      cur_x = new_x;
      if (i == iteration_times_ - 1)
        x_vec_[iteration_times_] = new_x;
    }
    if (debug_ && !alpha_iteration_flag){
      double traj_cost = calculateCostFunction();
      std::cout << "\n Cost: " << traj_cost << "\n\n";
    }
    infinite_feedback_update_flag_ = false;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getCurrentIdealPosition(double relative_time){
    int id;
    if (relative_time <= high_freq_end_time_)
      id = floor(relative_time * control_high_freq_);
    else
      id = int(high_freq_iteration_times_) +
        floor((relative_time - high_freq_end_time_) * control_low_freq_);

    if (id > iteration_times_ - 1){
      id = iteration_times_ - 1;
    }
    return getAbsoluteState(&(x_vec_[id]));
  }

  VectorXd SlqFiniteDiscreteControlHydrus::highFrequencyFeedbackControl(double relative_time, VectorXd *cur_real_x_ptr){
    // save for infinite state
    if (!infinite_feedback_update_flag_){
      infinite_feedback_update_flag_ = true;
      xn_last_ = *xn_ptr_;
    }

    // relative_time is (current time - start time)
    int id;
    if (relative_time <= high_freq_end_time_)
      id = floor(relative_time * control_high_freq_);
    else
      id = int(high_freq_iteration_times_) +
        floor((relative_time - high_freq_end_time_) * control_low_freq_);

    if (id > iteration_times_ - 1){
      id = iteration_times_ - 1;
      return infiniteFeedbackControl(cur_real_x_ptr);
    }
    VectorXd new_u = VectorXd::Zero(u_size_);
    VectorXd cur_x = getRelativeState(cur_real_x_ptr);
    new_u = u_vec_[id] + alpha_candidate_ * u_fw_vec_[id] + K_vec_[id] * (cur_x - x_vec_[id]);
    checkControlInputFeasible(&new_u, id);
    VectorXd stable_u = un_vec_[id];
    new_u = new_u + stable_u;
    return new_u;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::infiniteFeedbackControl(VectorXd *cur_real_x_ptr){
    VectorXd new_u = VectorXd::Zero(u_size_);
    VectorXd cur_x = stateSubtraction(cur_real_x_ptr, &xn_last_);
    new_u = -(*IDlqr_F_ptr_) * cur_x;
    checkControlInputFeasible(&new_u, iteration_times_);
    new_u = new_u + stable_u_last_;
    return new_u;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::highFrequencyLQRFeedbackControl(double relative_time, VectorXd *cur_real_x_ptr){
    VectorXd new_u = VectorXd::Zero(u_size_);
    VectorXd cur_x = getRelativeState(cur_real_x_ptr);
    // relative_time is (current time - start time)
    int id;
    if (relative_time <= high_freq_end_time_)
      id = floor(relative_time * control_high_freq_);
    else
      id = int(high_freq_iteration_times_) +
        floor((relative_time - high_freq_end_time_) * control_low_freq_);

    if (id > iteration_times_ - 1){
      id = iteration_times_;
      new_u = -(*IDlqr_F_ptr_) * cur_x;
      checkControlInputFeasible(&new_u, id);
      new_u = new_u + stable_u_last_;
    }
    else{
      new_u = -lqr_F_vec_[iteration_times_ - 1 - id] * cur_x;
      checkControlInputFeasible(&new_u, id);
      new_u = new_u + un_vec_[id];
    }
    return new_u;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getStableThrust(int time_id){
    VectorXd stable_u = VectorXd::Zero(u_size_);
    VectorXd g = VectorXd::Zero(4);
    g(3) = hydrus_weight_ * 9.78;
    MatrixXd H = MatrixXd::Zero(u_size_, u_size_);
    for (int i = 0; i < 4; ++i)
      H(3, i) = 1;
    // momentum
    MatrixXd H_minor = MatrixXd::Zero(u_size_ - 1, u_size_);
    for (int i = 0; i < 4; ++i){
      MatrixXd u_param = MatrixXd::Zero(u_size_ - 1, u_size_);
      u_param(u_size_ - 2, i) = 1.0;
      H_minor +=
        S_operation(link_center_pos_local_vec_[time_id][i] - cog_pos_local_vec_[time_id])
        * u_param;
    }
    // z momentum
    for (int i = 0; i < 4; ++i)
      H_minor(2, i) += M_z_(i);
    H.block<3, 4>(0, 0) = H_minor;

    // momentum from multi-link model
    Eigen::Vector3d w = Eigen::Vector3d::Zero();
    Eigen::Vector3d dw;
    Eigen::Vector3d momentum = Eigen::Vector3d::Zero();
    VectorXd dq = joint_dt_vec_[time_id];
    VectorXd ddq = joint_ddt_vec_[time_id];
    for (int i = 0; i < n_links_; ++i){
      MatrixXd JW_mat = getJacobianW(i);
      Eigen::Vector3d wi = w + VectorXdTo3d(JW_mat * dq);
      momentum +=
        I_vec_[time_id][i] * JW_mat * ddq
        + wi.cross(VectorXdTo3d(I_vec_[time_id][i] * wi))
        + I_dt_vec_[time_id][i] * wi;
    }
    //for (int i = 0; i < 3; ++i)
    //g(i) = momentum(i);

    /* lagrange mothod */
    // issue: min u_t * u; constraint: g = H * u  (stable point)
    //lamda: [4:0]
    // u = H_t * lamba
    // (H_  * H_t) * lamda = g
    // u = H_t * (H_ * H_t).inv * g
    Eigen::FullPivLU<Eigen::MatrixXd> solver((H * H.transpose()));
    Eigen::VectorXd lamda;
    lamda = solver.solve(g);
    stable_u = H.transpose() * lamda;
    return stable_u;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getCurrentJoint(double time, int order){
    VectorXd joint = VectorXd::Zero(n_links_ - 1);
    // keep quadrotor model, neglecting time, order
    if (order == 0){
      for (int i = 0; i < n_links_ - 1; ++i)
        joint(i) = PI / 2.0;
    }
    return joint;

    // example: end time is 6s: [0, 5] 1.57; [5, 5.5] 1.57-3.14*(t-5.0)^2; [5.5, 6] 3.14*(t-6.0)^2
    // double action_period = 2.0;
    // double action_start_time = end_time_ - action_period - 1.0;
    // if (transform_movement_flag_){
    //   if (order == 0){
    //     if (time > action_start_time + action_period)
    //       joint(2) = 0.0;
    //     else if (time > action_start_time + action_period / 2.0)
    //       joint(2) = 3.14 * pow((time - action_period - action_start_time) / action_period, 2.0);
    //     else if(time > action_start_time)
    //       joint(2) = 1.57 - 3.14 * pow((time - action_start_time) / action_period, 2.0);
    //   }
    //   else if (order == 1){
    //     if (time > action_start_time + action_period)
    //       joint(2) = 0.0;
    //     else if (time > action_start_time + action_period / 2.0)
    //       joint(2) = 3.14 * 2 * (time - action_period - action_start_time) / (action_period * action_period);
    //     else if(time > action_start_time)
    //       joint(2) = -3.14 * 2 * (time - action_start_time) / (action_period * action_period);
    //   }
    //   else if (order == 2){
    //     if (time > action_start_time + action_period)
    //       joint(2) = 0.0;
    //     else if (time > action_start_time + action_period / 2.0)
    //       joint(2) = 3.14 * 2 / (action_period * action_period);
    //     else if(time > action_start_time)
    //       joint(2) = -3.14 * 2 / (action_period * action_period);
    //   }
    // }

    double action_period = 2.0;
    double action_start_time = end_time_ - action_period - 1.0;
    double start_ang = PI / 2.0;
    double end_ang = 0.0;
    // example: sin function
    if (transform_movement_flag_){
      if (order == 0){
        if (time > action_start_time + action_period)
          joint(2) = end_ang;
        else if(time > action_start_time)
          joint(2) = (start_ang + end_ang) / 2.0
            + (start_ang - end_ang) / 2.0
            * cos(PI / action_period * (time - action_start_time));
      }
      else if (order == 1){
        if (time > action_start_time + action_period)
          joint(2) = 0.0;
        else if(time > action_start_time)
          joint(2) = -(start_ang - end_ang) / 2.0
            * sin(PI / action_period * (time - action_start_time))
            * PI / action_period;
      }
      else if (order == 2){
        if (time > action_start_time + action_period)
          joint(2) = 0.0;
        else if(time > action_start_time)
          joint(2) = -(start_ang - end_ang) / 2.0
            * cos(PI / action_period * (time - action_start_time))
            * pow(PI / action_period, 2.0);
      }
    }

    return joint;
  }

  Eigen::Matrix3d SlqFiniteDiscreteControlHydrus::getCurrentRotationMatrix(Eigen::Vector3d euler_angle, int order){
    Eigen::Matrix3d rot = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d rot_x, rot_y, rot_z;
    rot_x << 1, 0, 0,
      0, cos(euler_angle(0)), -sin(euler_angle(0)),
      0, sin(euler_angle(0)), cos(euler_angle(0));
    rot_y << cos(euler_angle(1)), 0, sin(euler_angle(1)),
      0, 1, 0,
      -sin(euler_angle(1)), 0, cos(euler_angle(1));
    rot_z << cos(euler_angle(2)), -sin(euler_angle(2)), 0,
      sin(euler_angle(2)), cos(euler_angle(2)), 0,
      0, 0, 1;
    if (order == 0){
      rot = rot_z * rot_y * rot_x;
    }
    else if(order == 1){
      Eigen::Matrix3d rot_x_d, rot_y_d, rot_z_d;
      rot_x_d << 1, 0, 0,
        0, -sin(euler_angle(0)), -cos(euler_angle(0)),
        0, cos(euler_angle(0)), -sin(euler_angle(0));
      rot_y_d << -sin(euler_angle(1)), 0, cos(euler_angle(1)),
        0, 1, 0,
        -cos(euler_angle(1)), 0, -sin(euler_angle(1));
      rot_z_d << -sin(euler_angle(2)), -cos(euler_angle(2)), 0,
        cos(euler_angle(2)), -sin(euler_angle(2)), 0,
        0, 0, 1;
      rot = rot_z_d * rot_y * rot_x
        + rot_z * rot_y_d * rot_x
        + rot_z * rot_y * rot_x_d;
    }
    else{
      ROS_WARN("[getCurrentRotationMatrix] order is too high.");
    }
    return rot;
  }

  void SlqFiniteDiscreteControlHydrus::updateMatrixAB(int time_id){
    updateMatrixA(time_id);
    updateMatrixB(time_id);
  }

  void SlqFiniteDiscreteControlHydrus::updateMatrixA(int time_id){
    *A_ptr_ = MatrixXd::Zero(x_size_, x_size_);

    VectorXd *x_ptr = new VectorXd(x_size_); *x_ptr = getAbsoluteState(&(x_vec_[time_id]));
    VectorXd *u_ptr = new VectorXd(u_size_); *u_ptr = u_vec_[time_id] + un_vec_[time_id];
    VectorXd *joint_ptr = new VectorXd(n_links_ - 1); *joint_ptr = joint_vec_[time_id];

    /* x, y, z */
    (*A_ptr_)(P_X, V_X) = 1;
    (*A_ptr_)(P_Y, V_Y) = 1;
    (*A_ptr_)(P_Z, V_Z) = 1;

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += (*u_ptr)[i];

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
    Eigen::Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd R_e = MatrixXd::Zero(3, 3);
    R_e << 1, tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]), tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]), -sin((*x_ptr)[E_R]),
      0, sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    Eigen::Vector3d d_e_w_x = R_e * Eigen::Vector3d(1.0, 0, 0);
    Eigen::Vector3d d_e_w_y = R_e * Eigen::Vector3d(0, 1.0, 0);
    Eigen::Vector3d d_e_w_z = R_e * Eigen::Vector3d(0, 0, 1.0);
    MatrixXd R_e_r = MatrixXd::Zero(3, 3);
    R_e_r << 0, tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]), -tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]),
      0, -sin((*x_ptr)[E_R]), -cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), -sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    Eigen::Vector3d d_e_e_r = R_e_r * w;
    MatrixXd R_e_p = MatrixXd::Zero(3, 3);
    double d_cosp = sin((*x_ptr)[E_P]) / pow(cos((*x_ptr)[E_P]), 2.0);
    R_e_p << 0, sin((*x_ptr)[E_R]) / pow(cos((*x_ptr)[E_P]), 2.0), cos((*x_ptr)[E_R]) / pow(cos((*x_ptr)[E_P]), 2.0),
      0, 0, 0,
      0, sin((*x_ptr)[E_R]) * d_cosp, cos((*x_ptr)[E_R]) * d_cosp;
    Eigen::Vector3d d_e_e_p = R_e_p * w;
    for (int i = E_R; i <= E_Y; ++i){
      (*A_ptr_)(i, W_X) = d_e_w_x(i - E_R);
      (*A_ptr_)(i, W_Y) = d_e_w_y(i - E_R);
      (*A_ptr_)(i, W_Z) = d_e_w_z(i - E_R);
      (*A_ptr_)(i, E_R) = d_e_e_r(i - E_R);
      (*A_ptr_)(i, E_P) = d_e_e_p(i - E_R);
    }

    /* w_x, w_y, w_z */
    /* d w = I.inv() * (sigma ri.cross(fi) + [0;0;fi * M_z(i)] - Ii*Jq_i*ddq - wi.cross(Ii * wi) - dIi * wi) */
    /* d w_w = I.inv() * (sigma - (d wi).cross(Ii * wi) - wi.cross(Ii * dwi) - dIi * dwi) */
    w = Eigen::Vector3d((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    Eigen::Matrix3d I_sum = Eigen::Matrix3d::Zero();
    for (int i = 0; i < n_links_; ++i)
      I_sum += I_vec_[time_id][i];
    Eigen::Matrix3d I_inv = I_sum.inverse();
    std::vector<Eigen::Vector3d> d_w_w_i_vec;
    for (int i = 0; i < 3; ++i){
      Eigen::Vector3d d_w_w_i = Eigen::Vector3d::Zero();
      Eigen::Vector3d dwi = Eigen::Vector3d::Zero(); dwi(i) = 1.0;
      for (int j = 0; j < n_links_; ++j){
        Eigen::Vector3d wj = w + VectorXdTo3d(getJacobianW(j) * joint_dt_vec_[time_id]);
        d_w_w_i = d_w_w_i + (-dwi.cross(VectorXdTo3d(I_vec_[time_id][j] * wj))
                             - wj.cross(VectorXdTo3d(I_vec_[time_id][j] * dwi))
                             - VectorXdTo3d(I_dt_vec_[time_id][i] * dwi));
      }
      d_w_w_i_vec.push_back(I_inv * d_w_w_i);
    }

    for (int i = W_X; i <= W_Z; ++i){
      (*A_ptr_)(i, W_X) = d_w_w_i_vec[0](i - W_X);
      (*A_ptr_)(i, W_Y) = d_w_w_i_vec[1](i - W_X);
      (*A_ptr_)(i, W_X) = d_w_w_i_vec[2](i - W_X);
    }

    if (time_id < high_freq_iteration_times_)
      (*A_ptr_) = (*A_ptr_) / control_high_freq_ + MatrixXd::Identity(x_size_, x_size_);
    else
      (*A_ptr_) = (*A_ptr_) / control_low_freq_ + MatrixXd::Identity(x_size_, x_size_);
  }

  void SlqFiniteDiscreteControlHydrus::updateMatrixB(int time_id){
    *B_ptr_ = MatrixXd::Zero(x_size_, u_size_);

    VectorXd *x_ptr = new VectorXd(x_size_); *x_ptr = getAbsoluteState(&(x_vec_[time_id]));
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
    /* d w = I.inv() * (sigma ri.cross(fi) + [0;0;fi * M_z(i)] - Ii*Jq_i*ddq - wi.cross(Ii * wi) - dIi * wi) */
    /* d w_u_i = I.inv() * (ri.cross(d fi) + [0;0;d fi * M_z(i)]) */
    Eigen::Matrix3d I_sum = Eigen::Matrix3d::Zero();
    for (int i = 0; i < n_links_; ++i)
      I_sum += I_vec_[time_id][i];
    Eigen::Matrix3d I_inv = I_sum.inverse();
    for (int i = 0; i < n_links_; ++i){
      Eigen::Vector3d dw_u_i = I_inv *
        ((link_center_pos_local_vec_[time_id][i] - cog_pos_local_vec_[time_id])
         .cross(Eigen::Vector3d(0, 0, 1.0))
         + Eigen::Vector3d(0, 0, M_z_(i)));
      for (int j = 0; j < 3; ++j)
        (*B_ptr_)(W_X + j, U_1 + i) = dw_u_i(j);
    }

    if (time_id < high_freq_iteration_times_)
      (*B_ptr_) = (*B_ptr_) / control_high_freq_;
    else
      (*B_ptr_) = (*B_ptr_) / control_low_freq_;
  }

  void SlqFiniteDiscreteControlHydrus::updateNewState(VectorXd *new_relative_x_ptr, VectorXd *relative_x_ptr, VectorXd *relative_u_ptr, int time_id){
    VectorXd dev_x = VectorXd::Zero(x_size_);
    VectorXd *x_ptr = new VectorXd(x_size_);
    VectorXd *u_ptr = new VectorXd(u_size_);
    *x_ptr = getAbsoluteState(relative_x_ptr);
    *u_ptr = *relative_u_ptr + un_vec_[time_id];

    VectorXd *joint_ptr = new VectorXd(n_links_ - 1);
    *joint_ptr = joint_vec_[time_id];

    /* x, y, z */
    dev_x(P_X) = (*x_ptr)(V_X);
    dev_x(P_Y) = (*x_ptr)(V_Y);
    dev_x(P_Z) = (*x_ptr)(V_Z);

    /* v_x, v_y, v_z */
    double u = 0.0;
    for (int i = 0; i < u_size_; ++i)
      u += (*u_ptr)[i];
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
    Eigen::Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd R_e = MatrixXd::Zero(3, 3);
    R_e << 1, tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]), tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]), -sin((*x_ptr)[E_R]),
      0, sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    Eigen::Vector3d d_e = R_e * w;
    for (int i = E_R; i <= E_Y; ++i)
      dev_x(i) = d_e(i - E_R);

    /* w_x, w_y, w_z */
    /* d w = I.inv() * (sigma ri.cross(fi) + [0;0;fi * M_z(i)] - Ii*Jq_i*ddq - wi.cross(Ii * wi) - dIi * wi) */
    Eigen::Vector3d dw;
    Eigen::Vector3d mid_result = Eigen::Vector3d::Zero();
    VectorXd dq = joint_dt_vec_[time_id];
    VectorXd ddq = joint_ddt_vec_[time_id];
    for (int i = 0; i < n_links_; ++i){
      MatrixXd JW_mat = getJacobianW(i);
      Eigen::Vector3d wi = w + VectorXdTo3d(JW_mat * dq);
      double fi = (*u_ptr)[i];
      mid_result +=
        (link_center_pos_local_vec_[time_id][i] - cog_pos_local_vec_[time_id]).
        cross(Eigen::Vector3d(0, 0, fi))
        + Eigen::Vector3d(0, 0, fi * M_z_(i))
        - I_vec_[time_id][i] * JW_mat * ddq
        - wi.cross(VectorXdTo3d(I_vec_[time_id][i] * wi))
        - I_dt_vec_[time_id][i] * wi;
    }
    Eigen::Matrix3d I_sum = Eigen::Matrix3d::Zero();
    for (int i = 0; i < n_links_; ++i)
      I_sum += I_vec_[time_id][i];
    dw = I_sum.inverse() * (mid_result);
    for (int i = 0; i < 3; ++i)
      dev_x(W_X + i) = dw(i);

    VectorXd new_x;
    if (time_id < high_freq_iteration_times_)
      new_x = dev_x / control_high_freq_ + *x_ptr;
    else
      new_x = dev_x / control_low_freq_ + *x_ptr;
    *new_relative_x_ptr = getRelativeState(&new_x);
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
    std::vector<Eigen::Matrix3d> cur_I_vec;
    std::vector<Eigen::Matrix3d> cur_I_dt_vec;
    for (int i = 0; i < n_links_; ++i){
      Eigen::Matrix3d cur_I = *I_ptr_;
      Eigen::Vector3d center_pos = link_center_pos_local_vec_[time_id][i] - cog_pos_local_vec_[time_id];
      cur_I(0, 0) += link_weight_vec_[i] * (pow(center_pos(1), 2.0)
                                            + pow(center_pos(2), 2.0));
      cur_I(1, 1) += link_weight_vec_[i] * (pow(center_pos(0), 2.0)
                                            + pow(center_pos(2), 2.0));
      cur_I(2, 2) += link_weight_vec_[i] * (pow(center_pos(0), 2.0)
                                            + pow(center_pos(1), 2.0));
      cur_I(0, 1) = -link_weight_vec_[i] * center_pos(0) * center_pos(1);
      cur_I(0, 2) = -link_weight_vec_[i] * center_pos(0) * center_pos(2);
      cur_I(1, 2) = -link_weight_vec_[i] * center_pos(1) * center_pos(2);
      cur_I(1, 0) = cur_I(0, 1);
      cur_I(2, 0) = cur_I(0, 2);
      cur_I(2, 1) = cur_I(1, 2);

      cur_I_vec.push_back(cur_I);

      Eigen::Matrix3d cur_I_dt = Eigen::Matrix3d::Zero();
      Eigen::Vector3d center_pos_dt = link_center_pos_local_dt_vec_[time_id][i] - cog_pos_local_dt_vec_[time_id];
      cur_I_dt(0, 0) = link_weight_vec_[i] * (2 * center_pos(1) * center_pos_dt(1)
                                              + 2 * center_pos(2) * center_pos_dt(2));
      cur_I_dt(1, 1) = link_weight_vec_[i] * (2 * center_pos(0) * center_pos_dt(0)
                                              + 2 * center_pos(2) * center_pos_dt(2));
      cur_I_dt(2, 2) = link_weight_vec_[i] * (2 * center_pos(0) * center_pos_dt(0)
                                              + 2 * center_pos(1) * center_pos_dt(1));
      cur_I_dt(0, 1) = -link_weight_vec_[i] * (center_pos(0) * center_pos_dt(1)
                                               + center_pos_dt(0) * center_pos(1));
      cur_I_dt(0, 2) = -link_weight_vec_[i] * (center_pos(0) * center_pos_dt(2)
                                               + center_pos_dt(0) * center_pos(2));
      cur_I_dt(1, 2) = -link_weight_vec_[i] * (center_pos(1) * center_pos_dt(2)
                                               + center_pos_dt(1) * center_pos(2));
      cur_I_dt(1, 0) = cur_I_dt(0, 1);
      cur_I_dt(2, 0) = cur_I_dt(0, 2);
      cur_I_dt(2, 1) = cur_I_dt(1, 2);
      cur_I_dt_vec.push_back(cur_I_dt);
    }
    I_vec_.push_back(cur_I_vec);
    I_dt_vec_.push_back(cur_I_dt_vec);
  }

  void SlqFiniteDiscreteControlHydrus::updateHydrusCogPosition(int time_id){
    Eigen::Vector3d cog_local_pos = Eigen::Vector3d::Zero();
    std::vector<Eigen::Vector3d> center_local_pos_vec = link_center_pos_local_vec_[time_id];
    for (int i = 0; i < n_links_; ++i){
      cog_local_pos = cog_local_pos + link_weight_vec_[i] * center_local_pos_vec[i];
    }
    cog_local_pos = cog_local_pos / hydrus_weight_;
    cog_pos_local_vec_.push_back(cog_local_pos);
  }

  void SlqFiniteDiscreteControlHydrus::updateHydrusCogPositionDerivative(int time_id){
    Eigen::Vector3d cog_local_pos_dt = Eigen::Vector3d::Zero();
    std::vector<Eigen::Vector3d> center_local_pos_dt_vec = link_center_pos_local_dt_vec_[time_id];
    for (int i = 0; i < n_links_; ++i){
      cog_local_pos_dt = cog_local_pos_dt + link_weight_vec_[i] * center_local_pos_dt_vec[i];
    }
    cog_local_pos_dt = cog_local_pos_dt / hydrus_weight_;
    cog_pos_local_dt_vec_.push_back(cog_local_pos_dt);
  }

  void SlqFiniteDiscreteControlHydrus::getHydrusLinksCenter(VectorXd *joint_ptr){
    std::vector<Eigen::Vector3d> links_center_vec;
    Eigen::Vector3d link1_center(link_length_ / 2.0, 0, 0);
    link1_center = link1_center;
    links_center_vec.push_back(link1_center);
    Eigen::Vector3d prev_link_end(link_length_, 0, 0);
    prev_link_end = prev_link_end;
    double joint_ang = 0.0;

    // only considering 2d hydrus
    for (int i = 1; i < n_links_; ++i){
      Eigen::Vector3d link_center = Eigen::Vector3d::Zero();
      joint_ang += (*joint_ptr)(i - 1);
      Eigen::Matrix3d rot;
      rot << cos(joint_ang), -sin(joint_ang), 0,
        sin(joint_ang), cos(joint_ang), 0,
        0, 0, 1;
      link_center = prev_link_end + rot * Eigen::Vector3d(link_length_ / 2.0, 0, 0);
      links_center_vec.push_back(link_center);
      prev_link_end = prev_link_end + rot * Eigen::Vector3d(link_length_, 0, 0);
    }
    link_center_pos_local_vec_.push_back(links_center_vec);
  }

  void SlqFiniteDiscreteControlHydrus::getHydrusLinksCenterDerivative(VectorXd *joint_ptr, VectorXd *joint_dt_ptr){
    std::vector<Eigen::Vector3d> links_center_dt_vec;
    Eigen::Vector3d link1_center_dt(0.0, 0, 0);
    links_center_dt_vec.push_back(link1_center_dt);
    Eigen::Vector3d prev_link_end_dt(0, 0, 0);
    double joint_ang = 0.0;
    double joint_ang_dt = 0.0;
    // only considering 2d hydrus
    for (int i = 1; i < n_links_; ++i){
      Eigen::Vector3d link_center_dt = Eigen::Vector3d::Zero();
      joint_ang += (*joint_ptr)(i - 1);
      joint_ang_dt += (*joint_dt_ptr)(i - 1);
      Eigen::Matrix3d rot_dt;
      rot_dt << -sin(joint_ang), -cos(joint_ang), 0,
        cos(joint_ang), -sin(joint_ang), 0,
        0, 0, 0;
      rot_dt = rot_dt * joint_ang_dt;
      link_center_dt = prev_link_end_dt + rot_dt * Eigen::Vector3d(link_length_ / 2.0, 0, 0);
      links_center_dt_vec.push_back(link_center_dt);
      prev_link_end_dt = prev_link_end_dt + rot_dt * Eigen::Vector3d(link_length_, 0, 0);
    }
    link_center_pos_local_dt_vec_.push_back(links_center_dt_vec);
  }

  MatrixXd SlqFiniteDiscreteControlHydrus::getJacobianW(int id){
    // eg. [0 0 0; 0 0 0; 0 0 0], [0 0 0; 0 0 0; 1 0 0],
    // [0 0 0; 0 0 0; 1 1 0], [0 0 0; 0 0 0; 1 1 1 ]
    MatrixXd JW_mat = MatrixXd::Zero(3, n_links_ - 1);
    for (int i = 0; i < id; ++i)
      JW_mat(2, i) = 1.0;
    return JW_mat;
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
      (*W_ptr)(j, j) = (*Q0_ptr_)(j, j) * weight;
    for (int j = 6; j < x_size_; ++j)
      (*W_ptr)(j, j) = (*Q0_ptr_)(j, j) * weight;

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

    MatrixXd P(x_size_, x_size_);
    P = *Q_ptr_;
    for (int i = 0; i < iteration_times_; ++i){
      MatrixXd F = MatrixXd::Zero(u_size_, x_size_);
      F = ((*R_ptr_) + B_ptr_->transpose() * P * (*B_ptr_)).inverse()
        * (B_ptr_->transpose() * P * (*A_ptr_));
      P = A_ptr_->transpose() * P * (*A_ptr_)
        - (A_ptr_->transpose() * P * (*B_ptr_)) * F
        + (*Q_ptr_);
      lqr_F_vec_.push_back(F);
    }

    VectorXd x = x_vec_[0];
    for (int i = iteration_times_ - 1; i >= 0; --i){
      VectorXd u = -lqr_F_vec_[i] * x;
      // Guarantee control is in bound
      checkControlInputFeasible(&u, i);

      VectorXd new_x(x_size_);
      updateNewState(&new_x, &x, &u, i);
      x = new_x;
      x_vec_[iteration_times_ - i] = x;
      u_vec_[iteration_times_ - i] = u;

      if ((i % 100 == 0 || i == iteration_times_ - 1) && debug_){
        printStateInfo(&x, i);
        printControlInfo(&u, i);
      }
    }
    ROS_INFO("[SLQ] LQR init finished");
  }

  double SlqFiniteDiscreteControlHydrus::calculateCostFunction(){
    double cost = 0.0;
    for (int i = 0; i < iteration_times_; ++i){
      // normal cost
      cost += (u_vec_[i].transpose() * (*R_ptr_) * u_vec_[i]
               + x_vec_[i].transpose() * (*Q0_ptr_) * x_vec_[i])(0);
      // waypoint cost
      VectorXd real_x = getAbsoluteState(&(x_vec_[i]));
      for (int j = 1; j < waypoints_ptr_->size(); ++j){
        MatrixXd W = MatrixXd::Zero(x_size_, x_size_);
        double cur_time;
        if (i <= high_freq_iteration_times_)
          cur_time = double(i) / control_high_freq_;
        else
          cur_time = high_freq_end_time_ +
            (i - high_freq_iteration_times_) / control_low_freq_;
        updateWaypointWeightMatrix(cur_time, (*time_ptr_)[j] - (*time_ptr_)[0], &W, j == (waypoints_ptr_->size() - 1));
        VectorXd dx_pt = stateSubtraction(&real_x, &((*waypoints_ptr_)[j]));
        cost += (dx_pt.transpose() * W * dx_pt)(0);
      }
    }
    // final time(tf) cost
    cost += (x_vec_[iteration_times_].transpose() * (*Riccati_P_ptr_) * x_vec_[iteration_times_])(0);
    return cost;
  }

  Eigen::Vector3d SlqFiniteDiscreteControlHydrus::VectorXdTo3d(VectorXd vec){
    Eigen::Vector3d vec3;
    for (int i = 0; i < 3; ++i)
      vec3(i) = vec(i);
    return vec3;
  }

  void SlqFiniteDiscreteControlHydrus::checkControlInputFeasible(VectorXd *u, int time_id){
    for (int j = 0; j < u_size_; ++j){
      if ((*u)(j) + un_vec_[time_id](j) < uav_rotor_thrust_min_)
        (*u)(j) = uav_rotor_thrust_min_ - un_vec_[time_id](j);
      else if ((*u)(j) + un_vec_[time_id](j) > uav_rotor_thrust_max_)
        (*u)(j) = uav_rotor_thrust_max_ - un_vec_[time_id](j);
    }
  }

  MatrixXd SlqFiniteDiscreteControlHydrus::S_operation(VectorXd vec){
    MatrixXd mat = MatrixXd::Zero(3, 3);
    mat << 0.0, -vec(2), vec(1),
      vec(2), 0.0, -vec(0),
      -vec(1), vec(0), 0.0;
    return mat;
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
      std::cout << (*u)(j) + un_vec_[id](j) << ", ";
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

  VectorXd SlqFiniteDiscreteControlHydrus::estimateFutureState(double relative_time){
    int id = floor(relative_time * control_freq_);
    if (id >= iteration_times_)
      return getAbsoluteState(&(x_vec_[iteration_times_]));
    else{
      VectorXd state_minor = getAbsoluteState(&(x_vec_[id]));
      VectorXd state_max = getAbsoluteState(&(x_vec_[id+1]));
      // simpliy average
      // return (state_minor + state_max) / 2.0;

      // get weight adding
      return (state_minor * ((id+1) / control_freq_ - relative_time)
              + state_max * (relative_time - id / control_freq_))
        * control_freq_;
    }
  }
}


