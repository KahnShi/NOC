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

    /* hydrus property */
    n_links_ = 4;
    I_ptr_ = new MatrixXd(3, 3);
    *I_ptr_ = MatrixXd::Zero(3, 3);
    (*I_ptr_)(0, 0) = 0.0001;
    (*I_ptr_)(1, 1) = 0.0001;
    (*I_ptr_)(2, 2) = 0.0002;

    /* link information */
    for (int i = 0; i < n_links_; ++i)
      link_mass_vec_.push_back(0.5);
    hydrus_mass_ = 0.0;
    for (int i = 0; i < n_links_; ++i)
      hydrus_mass_ += link_mass_vec_[i];
    link_length_ = 0.44;

    /* Initialize hydrus dynamic calculator */
    hydrus_dynamic_ptr_ = new HydrusDynamics(n_links_, link_length_, &link_mass_vec_, I_ptr_);

    /* rotor limits */
    uav_rotor_thrust_min_ = 0.0;
    uav_rotor_thrust_max_ = (hydrus_mass_ * 9.78 / 4.0) * 3.0;

    /* joint motor */
    joint_ptr_ = new VectorXd(3 * (n_links_-1)); // contains joint, d_joint, dd_joint

    /* Control information */
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
    for (int i = E_R; i <= E_Y; ++i)
      (*Q0_ptr_)(i, i) = 1.0;
    for (int i = DE_R; i <= DE_Y; ++i)
      (*Q0_ptr_)(i, i) = 1.0;
    // test: weight on z
    (*Q0_ptr_)(P_Z, P_Z) = (*Q0_ptr_)(V_Z, V_Z) = 100.0;
    *Q_ptr_ = (*Q0_ptr_);

    *R_ptr_ = 200 * MatrixXd::Identity(u_size_, u_size_);

    /* Assume initial and final state is still, namely dx = [v, a] = 0 */
    u0_ptr_ = new VectorXd(u_size_);
    un_ptr_ = new VectorXd(u_size_);
    for (int i = 0; i < 4; ++i){
      (*u0_ptr_)(i) = hydrus_mass_ * 9.78 / 4.0;
      (*un_ptr_)(i) = hydrus_mass_ * 9.78 / 4.0;
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

    line_search_steps_ = 4;

    debug_ = true;
    FDLQR();
    getRiccatiH();
    ROS_INFO("[SLQ] Initialization finished.");
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
      // assign value to joint
      *joint_ptr_ = VectorXd::Zero(3 * (n_links_-1));
      updateMatrixAB(x_ptr_, u_ptr_, joint_ptr_);

      *q_ptr_ = (*Q0_ptr_) * x_vec_[i];
      for (int j = 1; j < waypoints_ptr_->size(); ++j)
        *q_ptr_ = (*q_ptr_) +
          W_vec[j-1] * (*xn_ptr_ - (*waypoints_ptr_)[j] + x_vec_[i]);

      *r_ptr_ = (*R_ptr_) * (u_vec_[i]);
      // *r_ptr_ = (*R_ptr_) * ((u_vec_[i]) + (*un_ptr_));
      updateSLQEquations();

      VectorXd u_fb = (*K_ptr_) * (*x_ptr_);
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
          updateNewState(&new_x, &cur_x, &cur_u);
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
      updateNewState(&new_x, &cur_x, &cur_u);
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

  void SlqFiniteDiscreteControlHydrus::updateMatrixAB(VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr){
    hydrus_dynamic_ptr_->linaerizeState(x_ptr, u_ptr, joint_ptr, A_ptr_, B_ptr_);
    (*A_ptr_) = (*A_ptr_) / control_freq_ + MatrixXd::Identity(x_size_, x_size_);
  }

  void SlqFiniteDiscreteControlHydrus::updateNewState(VectorXd *new_x_ptr, VectorXd *x_ptr, VectorXd *u_ptr){
    // todo: examine
    VectorXd dev_x = hydrus_dynamic_ptr_->getStateDerivative();
    dev_x = dev_x / control_freq_;
    *new_x_ptr = dev_x + *x_ptr;
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
    // assign value to joint
    *joint_ptr_ = VectorXd::Zero(3 * (n_links_-1));
    updateMatrixAB(x_ptr_, u_ptr_, joint_ptr_);

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
    double fw_converge_threshold = (hydrus_mass_ * 9.78 / 4.0) * 0.1;
    if (fw_max < pow(fw_converge_threshold, 2))
      return true;
    else
      return false;
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
    (*W_ptr)(DE_Y, DE_Y) = weight;
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
    // assign value to joint
    *joint_ptr_ = VectorXd::Zero(3 * (n_links_-1));
    updateMatrixAB(x_ptr_, u_ptr_, joint_ptr_);

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
      updateNewState(&new_x, &x, &u);
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

  void SlqFiniteDiscreteControlHydrus::printMatrix(MatrixXd *mat_ptr, std::string mat_name){
    std::cout << "print " << mat_name << "\n";
    for (int i = 0; i < mat_ptr->rows(); ++i){
      for (int j = 0; j < mat_ptr->rows(); ++j){
        std::cout << (*mat_ptr)(i, j) << ", ";
      }
      std::cout << "\n";
    }
  }

  void SlqFiniteDiscreteControlHydrus::printMatrix(Matrix3d *mat_ptr, std::string mat_name){
    std::cout << "print " << mat_name << "\n";
    for (int i = 0; i < 3; ++i){
      for (int j = 0; j < 3; ++j){
        std::cout << (*mat_ptr)(i, j) << ", ";
      }
      std::cout << "\n";
    }
  }

  void SlqFiniteDiscreteControlHydrus::printVector(VectorXd *vec_ptr, std::string vec_name){
    std::cout << "print " << vec_name << "\n";
    for (int i = 0; i < vec_ptr->size(); ++i){
        std::cout << (*vec_ptr)(i) << ", ";
    }
    std::cout << "\n";
  }

  void SlqFiniteDiscreteControlHydrus::printVector(Vector3d *vec_ptr, std::string vec_name){
    std::cout << "print " << vec_name << "\n";
    for (int i = 0; i < 3; ++i){
        std::cout << (*vec_ptr)(i) << ", ";
    }
    std::cout << "\n";
  }

}


