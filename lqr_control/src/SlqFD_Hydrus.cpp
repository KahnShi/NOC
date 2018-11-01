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
  void SlqFiniteDiscreteControlHydrus::initHydrus(int baselink_id){
    verbose_ = false;
    baselink_id_ = baselink_id;
    std::cout << "[SlqFD_Hydrus] Baselink id: link " << baselink_id << "\n";

    /* Ros service */
    dare_client_ = nh_.serviceClient<lqr_control::Dare>("/dare_solver");

    not_first_slq_flag_ = false;
    /* ros param */
    nhp_.param("transform_movement_flag", transform_movement_flag_, true);
    nhp_.param("R_mid_para", R_mid_para_, 1.0);
    nhp_.param("Q_p_mid_para", Q_p_mid_para_, 50.0);
    nhp_.param("Q_v_mid_para", Q_v_mid_para_, 0.1);
    nhp_.param("Q_z_mid_para", Q_z_mid_para_, 50.0);
    nhp_.param("Q_w_mid_para", Q_w_mid_para_, 0.1);
    nhp_.param("Q_e_mid_para", Q_e_mid_para_, 5.0);
    nhp_.param("Q_yaw_mid_para", Q_yaw_mid_para_, 10.0);

    nhp_.param("manual_final_ocp_flag", manual_final_ocp_flag_, true);
    nhp_.param("Q_p_final_para", Q_p_final_para_, 500.0);
    nhp_.param("Q_v_final_para", Q_v_final_para_, 10.0);
    nhp_.param("Q_z_final_para", Q_z_final_para_, 500.0);
    nhp_.param("Q_w_final_para", Q_w_final_para_, 10.0);
    nhp_.param("Q_e_final_para", Q_e_final_para_, 300.0);
    nhp_.param("Q_yaw_final_para", Q_yaw_final_para_, 400.0);

    R_para_ = R_mid_para_;
    Q_p_para_ = Q_p_mid_para_;
    Q_v_para_ = Q_v_mid_para_;
    Q_z_para_ = Q_z_mid_para_;
    Q_w_para_ = Q_w_mid_para_;
    Q_e_para_ = Q_e_mid_para_;
    Q_yaw_para_ = Q_yaw_mid_para_;

    nhp_.param("verbose", debug_, false);
    nhp_.param("line_search_steps", line_search_steps_, 4);
    nhp_.param("line_search_mode", line_search_mode_, 2);// 0, standard mode; 1, final state priority mode; 2, final pose priority mode
    nhp_.param("n_links", n_links_, 4);
    nhp_.param("rotor_tilt_angle", rotor_tilt_ang_, -20.0);

    for (int i = 0; i < n_links_; ++i){
      link_weight_vec_.push_back(0.0);
      links_center_weight_link_frame_vec_.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
      links_center_on_link_frame_vec_.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
      s_tilts_.push_back(sin(pow(-1, i+1) * rotor_tilt_ang_ / 180.0 * PI));
      c_tilts_.push_back(cos(pow(-1, i+1) * rotor_tilt_ang_ / 180.0 * PI));
    }

    /* hydrus */
    double link_rod_mass, link_center_mass, joint_mass, protector_mass, protector_holder_mass, bat_mass, head_leg_mass, end_leg_mass, rotor_mass, fc_mass, gps_mass, pc_mass;
    nhp_.param("link_rod_mass", link_rod_mass, 0.099);
    nhp_.param("link_center_mass", link_center_mass, .211);
    nhp_.param("joint_mass", joint_mass, .19);
    nhp_.param("protector_mass", protector_mass, .05);
    nhp_.param("protector_holder_mass", protector_holder_mass, .0386);
    nhp_.param("bat_mass", bat_mass, 0.193);
    nhp_.param("head_leg_mass", head_leg_mass, 0.045);
    nhp_.param("end_leg_mass", end_leg_mass, 0.073);
    nhp_.param("rotor_mass", rotor_mass, .05);
    nhp_.param("fc_mass", fc_mass, 0.0457);
    nhp_.param("gps_mass", gps_mass, 0.044);
    nhp_.param("pc_mass", pc_mass, 0.164);

    // test
    link_center_mass += 0.03; // hard coding to make the mass equal

    double link_length, link_rod_length, protector_radius, link_joint_offset, rotor_radius, rotor_height, bat_offset, head_leg_offset, end_leg_offset;
    nhp_.param("link_length", link_length, .6);
    nhp_.param("link_rod_length", link_rod_length, .528);
    nhp_.param("protector_radius", protector_radius, .1925);
    nhp_.param("link_joint_offset", link_joint_offset, 0.3);
    nhp_.param("rotor_radius", rotor_radius, .02);
    nhp_.param("rotor_height", rotor_height, .05);
    nhp_.param("bat_offset", bat_offset, 0.0);
    nhp_.param("head_leg_offset", head_leg_offset, -0.25);
    nhp_.param("end_leg_offset", end_leg_offset, 0.235);

    /* create model from urdf file data */
    double protector_inertia = protector_mass * protector_radius * protector_radius;
    double link_rod_inertia = link_rod_mass * link_rod_length * link_rod_length /12;
    double protector_holder_inertia = protector_holder_mass * protector_radius * protector_radius;
    extraModel bat_model;
    bat_model.mass = bat_mass;
    bat_model.offset = Eigen::Vector3d(link_length / 2.0 + bat_offset, 0.0, -0.054);
    extraModel head_leg_model;
    head_leg_model.mass = head_leg_mass;
    head_leg_model.offset = Eigen::Vector3d(link_length * 0.5 + head_leg_offset, 0.0, -0.03);
    extraModel end_leg_model;
    end_leg_model.mass = end_leg_mass;
    end_leg_model.offset = Eigen::Vector3d(link_length * 0.5 + end_leg_offset, 0.0, -0.03);
    extraModel joint_model;
    joint_model.mass = joint_mass;
    joint_model.offset = Eigen::Vector3d(link_length * 0.5 + link_joint_offset, 0.0, 0.0);
    for (int i = 0; i < n_links_; ++i){
      hydrus_model_.link_vec.push_back(linkModel());
            hydrus_model_.link_vec[i].mass = link_rod_mass + link_center_mass + protector_mass + protector_holder_mass - rotor_mass + rotor_mass;
            hydrus_model_.link_vec[i].length = 0.6;
            hydrus_model_.link_vec[i].inertia_v = Eigen::Vector3d(protector_inertia / 2,
                                                                  link_rod_inertia + protector_holder_inertia + protector_inertia /2,
                                                                  link_rod_inertia + protector_holder_inertia + protector_inertia);
            hydrus_model_.link_vec[i].extra_vec.push_back(bat_model);
            if (i == 0)
              hydrus_model_.link_vec[i].extra_vec.push_back(head_leg_model);
            if (i == n_links_ - 1)
              hydrus_model_.link_vec[i].extra_vec.push_back(end_leg_model);
            if (i < n_links_ - 1)
              hydrus_model_.link_vec[i].extra_vec.push_back(joint_model);
            if (i == baselink_id_){
              extraModel fc_model;
              fc_model.mass = fc_mass;
              fc_model.offset = Eigen::Vector3d(link_length / 2.0 + 0.2281, -4.4*0.001, 0.03533)
                + Eigen::Vector3d(8.7*0.001, 4.4*0.001, -0.01);
              hydrus_model_.link_vec[i].extra_vec.push_back(fc_model);
              extraModel gps_model;
              gps_model.mass = gps_mass;
              gps_model.offset = Eigen::Vector3d(link_length / 2.0 + 0.206, 0.0, 0.18)
                + Eigen::Vector3d(0.0, 0.0, -0.08);
              hydrus_model_.link_vec[i].extra_vec.push_back(gps_model);
              extraModel pc_model;
              pc_model.mass = pc_mass;
              pc_model.offset = Eigen::Vector3d(link_length / 2.0 + 0.1565 + 0.063, 0.0, -0.0394)
                + Eigen::Vector3d(0.03*0.001, 0.72*0.001, 2.9 * 0.001);
              hydrus_model_.link_vec[i].extra_vec.push_back(pc_model);
            }
    }
    for (int i = 0; i < n_links_; ++i){
      link_weight_vec_[i] += hydrus_model_.link_vec[i].mass;
      links_center_weight_link_frame_vec_[i] = hydrus_model_.link_vec[i].mass * Eigen::Vector3d(hydrus_model_.link_vec[i].length / 2.0, 0.0, 0.0);
      for (int j = 0; j < hydrus_model_.link_vec[i].extra_vec.size(); ++j){
        link_weight_vec_[i] += hydrus_model_.link_vec[i].extra_vec[j].mass;
        links_center_weight_link_frame_vec_[i] += hydrus_model_.link_vec[i].extra_vec[j].mass * hydrus_model_.link_vec[i].extra_vec[j].offset;
      }
      links_center_on_link_frame_vec_[i] = links_center_weight_link_frame_vec_[i] / link_weight_vec_[i];
    }

    link_length_ = link_length;
    hydrus_weight_ = 0.0;
    for (int i = 0; i < n_links_; ++i)
      hydrus_weight_ += link_weight_vec_[i];
    // test
    std::cout << "\n\nhydrus_weight: " << hydrus_weight_ << "\n\n\n\n";
    propeller_pos_cog_offset_ = Eigen::Vector3d(0.0, 0.0, 0.10); // todo: be more accurate
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
    P0_ptr_ = new MatrixXd(x_size_, x_size_);
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
      (*Q0_ptr_)(i, i) = Q_p_para_;
    for (int i = V_X; i <= V_Z; ++i)
      (*Q0_ptr_)(i, i) = Q_v_para_;
    for (int i = W_X; i <= W_Z; ++i)
      (*Q0_ptr_)(i, i) = Q_w_para_;
    for (int i = E_R; i <= E_Y; ++i)
      (*Q0_ptr_)(i, i) = Q_e_para_;
    // test: weight on z
    (*Q0_ptr_)(P_Z, P_Z) = (*Q0_ptr_)(V_Z, V_Z) = Q_z_para_;

    (*R_ptr_).noalias() = R_para_ * MatrixXd::Identity(u_size_, u_size_);

    uav_rotor_thrust_min_ = 2.0;
    uav_rotor_thrust_max_ = 16.4;

    if (verbose_)
      ROS_INFO("[SLQ] Hydrus init finished.");
  }

  void SlqFiniteDiscreteControlHydrus::initSLQ(double freq, std::vector<double> *time_ptr, std::vector<VectorXd> *waypoints_ptr, TennisTaskDescriptor task_descriptor){
    if (verbose_)
      ROS_INFO("[SLQ] InitSLQ starts.");
    control_freq_ = freq;
    tennis_task_descriptor_ = task_descriptor;
    // Here we assume frequency is an odd integer
    slq_discrete_freq_ = freq;
    double period = (*time_ptr)[time_ptr->size() - 1] - (*time_ptr)[0];
    if (floor(slq_discrete_freq_ * period) < slq_discrete_freq_ * period){
      end_time_ = (floor(slq_discrete_freq_ * period) + 1.0) / slq_discrete_freq_;
      iteration_times_ = floor(slq_discrete_freq_ * period) + 1;
    }
    else{
      end_time_ = period;
      iteration_times_ = floor(slq_discrete_freq_ * period);
    }

    if (debug_){
      std::cout << "[SLQ] Trajectory period: " << end_time_
                << ", Itetation times: " << iteration_times_ << "\n";
      std::cout << "[SLQ] Start position: " << (*waypoints_ptr)[0].transpose() << "\n";
      std::cout << "[SLQ] End position: " << (*waypoints_ptr)[waypoints_ptr->size()  - 1].transpose() << "\n";
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

    x_vec_.resize(iteration_times_ + 1);
    u_vec_.resize(iteration_times_ + 1);
    joint_vec_.resize(iteration_times_ + 1);
    joint_dt_vec_.resize(iteration_times_ + 1);
    joint_ddt_vec_.resize(iteration_times_ + 1);
    I_vec_.resize(iteration_times_ + 1);
    I_dt_vec_.resize(iteration_times_ + 1);
    link_center_pos_local_vec_.resize(iteration_times_ + 1);
    link_center_pos_local_dt_vec_.resize(iteration_times_ + 1);
    link_center_pos_cog_vec_.resize(iteration_times_ + 1);
    cog_pos_local_vec_.resize(iteration_times_ + 1);
    cog_pos_local_dt_vec_.resize(iteration_times_ + 1);
    u_fw_vec_.resize(iteration_times_ + 1);
    K_vec_.resize(iteration_times_ + 1);
    un_vec_.resize(iteration_times_ + 1);
    lqr_F_vec_.resize(iteration_times_ + 1);

    if (verbose_)
      ROS_INFO("[SLQ] Assign vector starts.");
    #pragma omp parallel num_threads(4)
    {
      #pragma omp for
      for (int i = 0; i <= iteration_times_; ++i){
        double cur_time;
        cur_time = double(i) / slq_discrete_freq_;
        VectorXd cur_joint = getCurrentJoint(cur_time);
        VectorXd cur_joint_dt = getCurrentJoint(cur_time, 1);
        VectorXd cur_joint_ddt = getCurrentJoint(cur_time, 2);
        joint_vec_[i] = (cur_joint);
        joint_dt_vec_[i] = (cur_joint_dt);
        joint_ddt_vec_[i] = (cur_joint_ddt);
        getHydrusLinksCenter(&cur_joint, i);
        getHydrusLinksCenterDerivative(&cur_joint, &cur_joint_dt, i);
        getHydrusInertialTensor(&cur_joint, i);
        u_fw_vec_[i] = (u_init);
        K_vec_[i] = (MatrixXd::Zero(u_size_, x_size_));
        VectorXd stable_u = getStableThrust(i);
        un_vec_[i] = (stable_u);
        x_vec_[i] = (x_init);
      }
    }
    for (int i = 0; i <= iteration_times_; ++i)
      // u_vec_[i] = (un_vec_[0] - un_vec_[i]);
      u_vec_[i] = u_init; // all 0
    stable_u_last_ = un_vec_[iteration_times_];
    if (verbose_)
      ROS_INFO("[SLQ] Assign vector finished.");

    FDLQR();
    getRiccatiH();
    (*IDlqr_F_ptr_).noalias() = (*R_ptr_ +
                                 B_ptr_->transpose() * (*Riccati_P_ptr_) * (*B_ptr_)).inverse()
      * (B_ptr_->transpose() * (*Riccati_P_ptr_) * (*A_ptr_));

    if (manual_final_ocp_flag_){
      *P0_ptr_ = MatrixXd::Zero(x_size_, x_size_);
      for (int i = 0; i <= P_Z; ++i)
        (*P0_ptr_)(i, i) = Q_p_final_para_;
      for (int i = V_X; i <= V_Z; ++i)
        (*P0_ptr_)(i, i) = Q_v_final_para_;
      for (int i = W_X; i <= W_Z; ++i)
        (*P0_ptr_)(i, i) = Q_w_final_para_;
      for (int i = E_R; i <= E_Y; ++i)
        (*P0_ptr_)(i, i) = Q_e_final_para_;
      // test: weight on z
      (*P0_ptr_)(E_Y, E_Y) = Q_yaw_final_para_;
      (*P0_ptr_)(P_Z, P_Z) = Q_z_final_para_;
    }
    else{
      *P0_ptr_ = *Riccati_P_ptr_;
      if (debug_)
        std::cout << "\nRiccati P matrix: \n" << *P0_ptr_ << "\n\n";
    }

    if (debug_){
      double lqr_cost = calculateCostFunction();
      std::cout << "\n Cost: " << lqr_cost << "\n\n";
    }
    if (verbose_)
      ROS_INFO("[SLQ] InitSLQ finished with %d iterations.", iteration_times_);
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
    if (verbose_)
      ROS_INFO("[SLQ] Matrix AB is sent to Riccati solver.");
    if (dare_client_.call(dare_srv))
      dat_P = dare_srv.response.P;
    else
      ROS_ERROR("[SLQ] No response from dare sever.");
    for (int i = 0; i < x_size_; ++i)
      for (int j = 0; j < x_size_; ++j)
        (*Riccati_P_ptr_)(i, j) = dat_P.array.data[i * x_size_ + j];
    if (verbose_)
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

    if (verbose_)
      ROS_INFO("[SLQ] Get P matrice initial value from Ricatti function.");
  }

  void SlqFiniteDiscreteControlHydrus::iterativeOptimization(){
    *P_ptr_ = *P0_ptr_;
    //todo: judge the negative sign
    (*p_ptr_).noalias() = 2.0 * (*P_ptr_) * x_vec_[iteration_times_];

    for (int i = iteration_times_ - 1; i >= 0; --i){
      // add weight for waypoints
      double cur_time;
      cur_time = double(i) / slq_discrete_freq_;
      std::vector<MatrixXd> W_vec;
      for (int j = 1; j < waypoints_ptr_->size() - 1; ++j){
        MatrixXd W = MatrixXd::Zero(x_size_, x_size_);
        updateWaypointWeightMatrix(cur_time, (*time_ptr_)[j] - (*time_ptr_)[0], &W, j == (waypoints_ptr_->size() - 1));
        W_vec.push_back(W);
      }

      // update current Q and R matrix
      *Q_ptr_ = (*Q0_ptr_);
      for (int j = 1; j < waypoints_ptr_->size() - 1; ++j)
        (*Q_ptr_).noalias() += W_vec[j-1];

      *x_ptr_ = x_vec_[i];
      *u_ptr_ = u_vec_[i];
      *joint_ptr_ = joint_vec_[i];
      updateMatrixAB(i);

      //todo: judge the negative sign
      (*q_ptr_).noalias() = 2.0 * (*Q0_ptr_) * x_vec_[i];
      for (int j = 1; j < waypoints_ptr_->size() - 1; ++j)
        (*q_ptr_).noalias() +=
          2.0 * W_vec[j-1] * stateSubtraction(xn_ptr_, &((*waypoints_ptr_)[j]));

      //todo: judge the positive sign
      (*r_ptr_).noalias() = 2.0 * (*R_ptr_) * u_vec_[i];
      updateSLQEquations();

      u_fw_vec_[i] = (*l_ptr_);
      //todo: judge the positive sign
      K_vec_[i] = (*K_ptr_);
    }

    /* Update control by finding the best alpha */
    if (verbose_)
      ROS_INFO("[SLQ] Line search starts.");

    alpha_ = 1.0;
    alpha_candidate_ = 1.0;
    double energy_min = -1.0;
    double search_rate = 2.0;
    // openmp para
    std::vector<double> alpha_vec(line_search_steps_);
    std::vector<double> energy_vec(line_search_steps_);
    for (int i = 0; i < line_search_steps_; ++i){
      alpha_vec[i] = 1.0 / pow(search_rate, i);
      energy_vec[i] = 0.0;
    }
    /* When there are no middle waypoints, feedforward term is 0. */
    bool alpha_iteration_flag = true;
    if (feedforwardConverged()){
      alpha_iteration_flag = false;
      if (debug_)
        std::cout << "[SLQ] feedforward converge.\n";
    }
    else{
      #pragma omp parallel num_threads(line_search_steps_)
      {
        int id = omp_get_thread_num();
        double energy_sum = 0.0;
        VectorXd cur_u(u_size_);
        VectorXd cur_x = x_vec_[0];
        for (int i = 0; i < iteration_times_; ++i){
          VectorXd cur_u(u_size_);
          cur_u = u_vec_[i];
          cur_u.noalias() += alpha_vec[id] * u_fw_vec_[i];
          cur_u.noalias() += K_vec_[i] * cur_x;
          checkControlInputFeasible(&cur_u, i);
          // calculate energy
          // add weight for waypoints
          if (line_search_mode_ == 0){
            double cur_time;
            cur_time = double(i) / slq_discrete_freq_;
            std::vector<MatrixXd> W_vec;
            for (int j = 1; j < waypoints_ptr_->size() - 1; ++j){
              MatrixXd W = MatrixXd::Zero(x_size_, x_size_);
              updateWaypointWeightMatrix(cur_time, (*time_ptr_)[j] - (*time_ptr_)[0], &W, j == (waypoints_ptr_->size() - 1));
              W_vec.push_back(W);
            }

            energy_sum += (cur_u.transpose() * (*R_ptr_) * cur_u)(0);
            energy_sum += (cur_x.transpose() * (*Q0_ptr_) * cur_x)(0);
            VectorXd real_x = getAbsoluteState(&cur_x);
            for (int j = 1; j < waypoints_ptr_->size() - 1; ++j){
              VectorXd dx_pt = stateSubtraction(&real_x, &((*waypoints_ptr_)[j]));
              energy_sum += (dx_pt.transpose() * W_vec[j-1] * dx_pt)(0);
            }
          }

          VectorXd new_x(x_size_);
          updateNewState(&new_x, &cur_x, &cur_u, i);
          cur_x = new_x;
        }
        // when line_search_mode is 2, only considering final pose cost
        if (line_search_mode_ == 2){
          for (int i = P_X; i <= P_Z; ++i)
            energy_sum += pow(cur_x(i), 2) * (*P0_ptr_)(i, i);
          for (int i = E_R; i <= E_Y; ++i)
            energy_sum += pow(cur_x(i), 2) * (*P0_ptr_)(i, i);
        }
        else
          energy_sum += (cur_x.transpose() * (*P0_ptr_) * cur_x)(0);
        energy_vec[id] = energy_sum;
      }
    }
    // energy and alpha' relationships
    if (debug_){
      for (int i = 0; i < line_search_steps_; ++i)
        std::cout << "[SLQ] Energy: " << energy_vec[i] << ", alpha: " << alpha_vec[i] << "\n";
    }
    for (int i = 0; i < line_search_steps_; ++i){
      if (energy_vec[i] < energy_min || energy_min < 0){
        energy_min = energy_vec[i];
        alpha_candidate_ = alpha_vec[i];
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
      cur_u = u_vec_[i];
      cur_u.noalias() += alpha_candidate_ * u_fw_vec_[i];
      cur_u.noalias() += K_vec_[i] * (cur_x);
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
      if (i == iteration_times_ - 1){
        x_vec_[iteration_times_] = new_x;
        if (debug_)
          printStateInfo(&new_x, iteration_times_);
      }
    }
    if (debug_ && !alpha_iteration_flag){
      double traj_cost = calculateCostFunction();
      std::cout << "\n Cost: " << traj_cost << "\n\n";
    }
    infinite_feedback_update_flag_ = false;
    if (verbose_)
      ROS_INFO("[SLQ] Line search finished.");
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getCurrentIdealPosition(double relative_time){
    int id;
    id = floor(relative_time * slq_discrete_freq_);
    if (id > iteration_times_ - 1){
      id = iteration_times_ - 1;
    }
    return getAbsoluteState(&(x_vec_[id]));
  }

  VectorXd SlqFiniteDiscreteControlHydrus::slqFeedbackControl(double relative_time, VectorXd *cur_real_x_ptr){
    // save for infinite state
    if (!infinite_feedback_update_flag_){
      infinite_feedback_update_flag_ = true;
      xn_last_ = *xn_ptr_;
    }

    // relative_time is (current time - start time)
    int id;
    id = floor(relative_time * slq_discrete_freq_);
    if (id > iteration_times_ - 1){
      id = iteration_times_ - 1;
      return infiniteFeedbackControl(cur_real_x_ptr);
    }
    VectorXd new_u = VectorXd::Zero(u_size_);
    VectorXd cur_x = getRelativeState(cur_real_x_ptr);
    new_u = u_vec_[id];
    new_u.noalias() += alpha_candidate_ * u_fw_vec_[id];
    new_u.noalias() += K_vec_[id] * stateSubtraction(&cur_x, &(x_vec_[id]));
    checkControlInputFeasible(&new_u, id);
    VectorXd stable_u = un_vec_[id];
    new_u.noalias() += stable_u;
    return new_u;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::infiniteFeedbackControl(VectorXd *cur_real_x_ptr){
    VectorXd new_u = VectorXd::Zero(u_size_);
    VectorXd cur_x = stateSubtraction(cur_real_x_ptr, &xn_last_);
    new_u.noalias() = -(*IDlqr_F_ptr_) * cur_x;
    checkControlInputFeasible(&new_u, iteration_times_);
    new_u.noalias() += stable_u_last_;
    return new_u;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getStableThrust(int time_id){
    VectorXd stable_u = VectorXd::Zero(u_size_);
    // todo: currently simply average of the gravity to save computation, configuration needs to be considered
    // for (int i = 0; i < u_size_; ++i)
    //   stable_u(i) = hydrus_weight_ * 9.78 / u_size_;
    // return stable_u;


    // Eigen::MatrixXd P_dash = Eigen::MatrixXd::Zero(3, n_links_);
    // Eigen::VectorXd p_x(n_links_), p_y(n_links_), p_c(n_links_), p_m(n_links_);

    // for(int i = 0; i < n_links_; i++){
    //   p_y(i) =  link_center_pos_cog_vec_[time_id][i](1);
    //   p_x(i) = -link_center_pos_cog_vec_[time_id][i](0);
    //   // p_c(i) =  rotor_direction.at(i + 1) * m_f_rate_ ;
    //   p_m(i) =  1.0 / hydrus_weight_;
    // }
    // P_dash.row(0) = p_y;
    // P_dash.row(1) = p_x;
    // P_dash.row(2) = p_m;
    // Eigen::VectorXd g3(3);
    // g3 << 0, 0, 9.8;
    // Eigen::FullPivLU<Eigen::MatrixXd> solver((P_dash * P_dash.transpose()));
    // Eigen::VectorXd lamda;
    // lamda = solver.solve(g3);
    // stable_u = P_dash.transpose() * lamda;
    // return stable_u;


    Eigen::Matrix4d P;
    Eigen::VectorXd yaw_vec(n_links_);
    yaw_vec(baselink_id_) = 0.0;
    for (int i = baselink_id_ - 1; i >= 0; --i)
      yaw_vec(i) = yaw_vec(i + 1) - joint_vec_[time_id](i);
    for (int i = baselink_id_ + 1; i < n_links_; ++i)
      yaw_vec(i) = yaw_vec(i - 1) + joint_vec_[time_id](i - 1);
    for (int i = 0; i < n_links_; ++i){
      Eigen::Vector3d mid_result;
      mid_result.noalias() =
        (link_center_pos_cog_vec_[time_id][i] + propeller_pos_cog_offset_).
        cross(Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i], sin(yaw_vec(i)) * s_tilts_[i], c_tilts_[i]))
        + M_z_(i) * Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i], sin(yaw_vec(i)) * s_tilts_[i], c_tilts_[i]);
      for (int j = 0; j < 3; ++j)
        P(j, i) = mid_result(j);
    }
    for (int i = 0; i < n_links_; ++i)
      P(n_links_ - 1, i) = c_tilts_[i];
    // test
    std::cout  << P << "\n\n";
    Eigen::VectorXd g3(3);
    g3 << 0, 0, 9.8 * hydrus_weight_;

    /* method 1: directly calculate inverse matrix, which might not feasible */
    // stable_u = P.inverse() * g3;

    /* method 2: considering torque: x,y; force: z */
    // Eigen::MatrixXd P_dash = Eigen::MatrixXd::Zero(3, n_links_);
    // P_dash.row(0) = P.row(0);
    // P_dash.row(1) = P.row(1);
    // P_dash.row(2) = P.row(3);
    // Eigen::FullPivLU<Eigen::MatrixXd> solver((P_dash * P_dash.transpose()));
    // Eigen::VectorXd lamda;
    // lamda = solver.solve(g3);
    // stable_u = P_dash.transpose() * lamda;
    // // test
    // std::cout << time_id << ": " << stable_u.transpose() << "\n";
    // return stable_u;

    /* method 3: considering torque: x,y,z; force: z */
    Eigen::FullPivLU<Eigen::MatrixXd> solver((P * P.transpose()));
    Eigen::VectorXd lamda;
    Eigen::VectorXd g4(4);
    g4 << 0, 0, 0, 9.8 * hydrus_weight_;
    lamda = solver.solve(g4);
    stable_u = P.transpose() * lamda;
    return stable_u;

    /* method 4: considering torque: x,y,z; force: x,y,z */
    // Eigen::MatrixXd P_ext(6, n_links_);
    // for (int i = 0; i < 4; ++i)
    //   P_ext.row(i) = P.row(i);
    // for (int i = 0; i < n_links_; ++i){
    //   P_ext(4, i) = cos(yaw_vec(i)) * s_tilts_[i];
    //   P_ext(5, i) = sin(yaw_vec(i)) * s_tilts_[i];
    // }
    // Eigen::FullPivLU<Eigen::MatrixXd> solver((P_ext * P_ext.transpose()));
    // Eigen::VectorXd lamda;
    // Eigen::VectorXd g6(6);
    // g6 << 0, 0, 0, 9.8 * hydrus_weight_, 0, 0;
    // lamda = solver.solve(g6);
    // stable_u = P_ext.transpose() * lamda;
    // return stable_u;
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getCurrentJoint(double relative_time, int order){
    double time = relative_time + (*time_ptr_)[0];
    return getCurrentJointAbsoluteTime(time, order);
  }

  VectorXd SlqFiniteDiscreteControlHydrus::getCurrentJointAbsoluteTime(double time, int order){
    VectorXd joint = VectorXd::Zero(n_links_ - 1);
    // test
    // if (order == 0)
    //   joint << 0.785, 1.5708, 0.785;
    // return joint;
    
    
    int joint_id;
    if (tennis_task_descriptor_.hitting_hand == 0) // left hand
      joint_id = 0;
    else if (tennis_task_descriptor_.hitting_hand == 1) // right hand
      joint_id = 2;
    else if (tennis_task_descriptor_.hitting_hand == 2){ // no hand
      if (order == 0) // no transformation for no hand cases
        joint << 0.785, 1.5708, 0.785;
      return joint;
    }
    double dq = 5.0 * tennis_task_descriptor_.hitting_time; // joint velocity
    double racket_return_time = 1.2; //tennis_task_descriptor_.post_hitting_time
    // todo: temprarily assume transform action only in hitting time
    if (time > tennis_task_descriptor_.hitting_time + racket_return_time){
      // keep quadrotor model, neglecting time, order
      if (order == 0){
        for (int i = 0; i < n_links_ - 1; ++i)
          joint(i) = PI / 2.0;
        joint << 0.785, 1.5708, 0.785;
      }
      return joint;
    }
    else if (time > tennis_task_descriptor_.hitting_time){
      double tf = racket_return_time;
      double relative_time = time - tennis_task_descriptor_.hitting_time;
      if (order == 0){
        joint << 0.785, 1.5708, 0.785;
        joint(joint_id) = dq / pow(tf, 2) * pow(relative_time, 3) - 2 * dq / tf * pow(relative_time, 2) + dq * relative_time + 0.785;
      }
      else if (order == 1)
        joint(joint_id) = 3 * dq / pow(tf, 2) * pow(relative_time, 2) - 4 * dq / tf * relative_time + dq;
      else if (order == 2)
        joint(joint_id) = 6 * dq / pow(tf, 2) * relative_time - 4 * dq / tf;
      return joint;
    }
    else{
      double tf = tennis_task_descriptor_.hitting_time;
      if (order == 0){
        joint << 0.785, 1.5708, 0.785;
        joint(joint_id) = dq / pow(tf, 2) * pow(time, 3) - dq / tf * pow(time, 2) + 0.785;
      }
      else if (order == 1)
        joint(joint_id) = 3 * dq / pow(tf, 2) * pow(time, 2) - 2 * dq / tf * time;
      else if (order == 2)
        joint(joint_id) = 6 * dq / pow(tf, 2) * time - 2 * dq / tf;
      return joint;
    }

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
    VectorXd *u_ptr = new VectorXd(u_size_); (*u_ptr).noalias() = u_vec_[time_id] + un_vec_[time_id];
    VectorXd *joint_ptr = new VectorXd(n_links_ - 1); *joint_ptr = joint_vec_[time_id];

    /* x, y, z */
    (*A_ptr_)(P_X, V_X) = 1;
    (*A_ptr_)(P_Y, V_Y) = 1;
    (*A_ptr_)(P_Z, V_Z) = 1;

    /* v_x, v_y, v_z */
    Eigen::Vector3d f_sum_cog = Eigen::Vector3d::Zero();
    Eigen::VectorXd yaw_vec(n_links_);
    yaw_vec(baselink_id_) = 0.0;
    for (int i = baselink_id_ - 1; i >= 0; --i)
      yaw_vec(i) = yaw_vec(i + 1) - joint_vec_[time_id](i);
    for (int i = baselink_id_ + 1; i < n_links_; ++i)
      yaw_vec(i) = yaw_vec(i - 1) + joint_vec_[time_id](i - 1);
    for (int i = 0; i < n_links_; ++i){
      double fi = (*u_ptr)[i];
      f_sum_cog.noalias() += Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i] * fi, sin(yaw_vec(i)) * s_tilts_[i] * fi, c_tilts_[i] * fi);
    }
    Eigen::Matrix3d cog_rot_er_dt, cog_rot_ep_dt, cog_rot_ey_dt;
    cog_rot_ey_dt << -sin((*x_ptr)[E_R]), -cos((*x_ptr)[E_R]), 0,
      cos((*x_ptr)[E_R]), -sin((*x_ptr)[E_R]), 0,
      0, 0, 0;
    cog_rot_ep_dt << -sin((*x_ptr)[E_P]), 0, cos((*x_ptr)[E_P]),
      0, 0, 0,
      -cos((*x_ptr)[E_P]), 0, -sin((*x_ptr)[E_P]);
    cog_rot_er_dt << 0, 0, 0,
      0, -sin((*x_ptr)[E_Y]), -cos((*x_ptr)[E_Y]),
      0, cos((*x_ptr)[E_Y]), -sin((*x_ptr)[E_Y]);
    Eigen::Matrix3d cog_rot_dey = cog_rot_ey_dt
      * Eigen::AngleAxisd((*x_ptr)[E_P], Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd((*x_ptr)[E_R], Eigen::Vector3d::UnitX());
    Eigen::Matrix3d cog_rot_dep = Eigen::AngleAxisd((*x_ptr)[E_Y], Eigen::Vector3d::UnitZ())
      * cog_rot_ep_dt
      * Eigen::AngleAxisd((*x_ptr)[E_R], Eigen::Vector3d::UnitX());
    Eigen::Matrix3d cog_rot_der = Eigen::AngleAxisd((*x_ptr)[E_Y], Eigen::Vector3d::UnitZ())
      * Eigen::AngleAxisd((*x_ptr)[E_P], Eigen::Vector3d::UnitY())
      * cog_rot_er_dt;
    Eigen::Vector3d f_sum_w_der = cog_rot_der * f_sum_cog;
    Eigen::Vector3d f_sum_w_dep = cog_rot_dep * f_sum_cog;
    Eigen::Vector3d f_sum_w_dey = cog_rot_dey * f_sum_cog;
    for (int i = 0; i < 3; ++i) (*A_ptr_)(V_X + i, E_R) = f_sum_w_der(i);
    for (int i = 0; i < 3; ++i) (*A_ptr_)(V_X + i, E_P) = f_sum_w_dep(i);
    for (int i = 0; i < 3; ++i) (*A_ptr_)(V_X + i, E_Y) = f_sum_w_dey(i);

    /* e_r, e_p, e_y */
    /* d e = R_e * w_b */
    Eigen::Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd R_e = MatrixXd::Zero(3, 3);
    R_e << 1, tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]), tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]), -sin((*x_ptr)[E_R]),
      0, sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    MatrixXd R_e_r = MatrixXd::Zero(3, 3);
    R_e_r << 0, tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]), -tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]),
      0, -sin((*x_ptr)[E_R]), -cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), -sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    Eigen::Vector3d d_e_e_r; d_e_e_r.noalias() = R_e_r * w;
    MatrixXd R_e_p = MatrixXd::Zero(3, 3);
    double d_cosp = sin((*x_ptr)[E_P]) / pow(cos((*x_ptr)[E_P]), 2.0);
    R_e_p << 0, sin((*x_ptr)[E_R]) / pow(cos((*x_ptr)[E_P]), 2.0), cos((*x_ptr)[E_R]) / pow(cos((*x_ptr)[E_P]), 2.0),
      0, 0, 0,
      0, sin((*x_ptr)[E_R]) * d_cosp, cos((*x_ptr)[E_R]) * d_cosp;
    Eigen::Vector3d d_e_e_p; d_e_e_p.noalias() = R_e_p * w;
    for (int i = E_R; i <= E_Y; ++i){
      (*A_ptr_)(i, W_X) = R_e(i - E_R, 0);
      (*A_ptr_)(i, W_Y) = R_e(i - E_R, 1);
      (*A_ptr_)(i, W_Z) = R_e(i - E_R, 2);
      (*A_ptr_)(i, E_R) = d_e_e_r(i - E_R);
      (*A_ptr_)(i, E_P) = d_e_e_p(i - E_R);
    }

    /* w_x, w_y, w_z */
    /* d w = I.inv() * (sigma ri.cross(fi) + [0;0;fi * M_z(i)] - Ii*Jq_i*ddq - wi.cross(Ii * wi) - dIi * wi) */
    /* d w_w = I.inv() * (sigma - (d wi).cross(Ii * wi) - wi.cross(Ii * dwi) - dIi * dwi) */
    w = Eigen::Vector3d((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    Eigen::Matrix3d I_sum = Eigen::Matrix3d::Zero();
    for (int i = 0; i < n_links_; ++i)
      I_sum.noalias() += I_vec_[time_id][i];
    Eigen::Matrix3d I_inv = I_sum.inverse();
    std::vector<Eigen::Vector3d> d_w_w_i_vec;
    for (int i = 0; i < 3; ++i){
      Eigen::Vector3d d_w_w_i = Eigen::Vector3d::Zero();
      Eigen::Vector3d dwi = Eigen::Vector3d::Zero(); dwi(i) = 1.0;
      for (int j = 0; j < n_links_; ++j){
        Eigen::Vector3d wj = w + VectorXdTo3d(getJacobianW(j) * joint_dt_vec_[time_id]);
        d_w_w_i.noalias() -= dwi.cross(VectorXdTo3d(I_vec_[time_id][j] * wj));
        d_w_w_i.noalias() -= wj.cross(VectorXdTo3d(I_vec_[time_id][j] * dwi));
        d_w_w_i.noalias() -= VectorXdTo3d(I_dt_vec_[time_id][i] * dwi);
      }
      d_w_w_i_vec.push_back(I_inv * d_w_w_i);
    }

    for (int i = W_X; i <= W_Z; ++i){
      (*A_ptr_)(i, W_X) = d_w_w_i_vec[0](i - W_X);
      (*A_ptr_)(i, W_Y) = d_w_w_i_vec[1](i - W_X);
      (*A_ptr_)(i, W_X) = d_w_w_i_vec[2](i - W_X);
    }

    (*A_ptr_) = (*A_ptr_) / slq_discrete_freq_ + MatrixXd::Identity(x_size_, x_size_);
  }

  void SlqFiniteDiscreteControlHydrus::updateMatrixB(int time_id){
    *B_ptr_ = MatrixXd::Zero(x_size_, u_size_);

    VectorXd *x_ptr = new VectorXd(x_size_); *x_ptr = getAbsoluteState(&(x_vec_[time_id]));
    VectorXd *joint_ptr = new VectorXd(n_links_ - 1); *joint_ptr = joint_vec_[time_id];

    /* x, y, z */
    /* all 0 */

    // /* v_x, v_y, v_z */
    Eigen::Matrix3d cog_rot;
    cog_rot = Eigen::AngleAxisd((*x_ptr)[E_Y], Eigen::Vector3d::UnitZ())
      * Eigen::AngleAxisd((*x_ptr)[E_P], Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd((*x_ptr)[E_R], Eigen::Vector3d::UnitX());
    Eigen::VectorXd yaw_vec(n_links_);
    yaw_vec(baselink_id_) = 0.0;
    for (int i = baselink_id_ - 1; i >= 0; --i)
      yaw_vec(i) = yaw_vec(i + 1) - joint_vec_[time_id](i);
    for (int i = baselink_id_ + 1; i < n_links_; ++i)
      yaw_vec(i) = yaw_vec(i - 1) + joint_vec_[time_id](i - 1);
    for (int i = 0; i < n_links_; ++i){
      Eigen::Vector3d fi_dt_cog = Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i], sin(yaw_vec(i)) * s_tilts_[i], c_tilts_[i]);
      Eigen::Vector3d fi_dt_w = cog_rot * fi_dt_cog;
      for (int j = 0; j < 3; ++j)
        (*B_ptr_)(V_X + j, U_1 + i) = fi_dt_w(j);
    }

    /* e_r, e_p, e_y */
    /* all 0 */

    /* w_x, w_y, w_z */
    /* d w = I.inv() * (sigma ri.cross(fi) + [0;0;fi * M_z(i)] - Ii*Jq_i*ddq - wi.cross(Ii * wi) - dIi * wi) */
    /* d w_u_i = I.inv() * (ri.cross(d fi) + [0;0;d fi * M_z(i)]) */
    Eigen::Matrix3d I_sum = Eigen::Matrix3d::Zero();
    for (int i = 0; i < n_links_; ++i)
      I_sum.noalias() += I_vec_[time_id][i];
    Eigen::Matrix3d I_inv = I_sum.inverse();
    for (int i = 0; i < n_links_; ++i){
      Eigen::Vector3d dw_u_i;
      dw_u_i.noalias() = I_inv *
        ((link_center_pos_cog_vec_[time_id][i] + propeller_pos_cog_offset_)
         .cross(Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i], sin(yaw_vec(i)) * s_tilts_[i], c_tilts_[i]))
         + Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i], sin(yaw_vec(i)) * s_tilts_[i], c_tilts_[i]) * M_z_(i));
      for (int j = 0; j < 3; ++j)
        (*B_ptr_)(W_X + j, U_1 + i) = dw_u_i(j);
    }

    (*B_ptr_) = (*B_ptr_) / slq_discrete_freq_;
  }

  void SlqFiniteDiscreteControlHydrus::updateNewState(VectorXd *new_relative_x_ptr, VectorXd *relative_x_ptr, VectorXd *relative_u_ptr, int time_id){
    VectorXd dev_x = VectorXd::Zero(x_size_);
    VectorXd *x_ptr = new VectorXd(x_size_);
    VectorXd *u_ptr = new VectorXd(u_size_);
    *x_ptr = getAbsoluteState(relative_x_ptr);
    (*u_ptr).noalias() = *relative_u_ptr + un_vec_[time_id];

    VectorXd *joint_ptr = new VectorXd(n_links_ - 1);
    *joint_ptr = joint_vec_[time_id];

    /* x, y, z */
    dev_x(P_X) = (*x_ptr)(V_X);
    dev_x(P_Y) = (*x_ptr)(V_Y);
    dev_x(P_Z) = (*x_ptr)(V_Z);

    /* v_x, v_y, v_z */
    Eigen::Vector3d f_sum_cog = Eigen::Vector3d::Zero();
    Eigen::VectorXd yaw_vec(n_links_);
    yaw_vec(baselink_id_) = 0.0;
    for (int i = baselink_id_ - 1; i >= 0; --i)
      yaw_vec(i) = yaw_vec(i + 1) - joint_vec_[time_id](i);
    for (int i = baselink_id_ + 1; i < n_links_; ++i)
      yaw_vec(i) = yaw_vec(i - 1) + joint_vec_[time_id](i - 1);
    for (int i = 0; i < n_links_; ++i){
      double fi = (*u_ptr)[i];
      f_sum_cog.noalias() += Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i] * fi, sin(yaw_vec(i)) * s_tilts_[i] * fi, c_tilts_[i] * fi);
    }
    Eigen::Matrix3d cog_rot;
    cog_rot = Eigen::AngleAxisd((*x_ptr)[E_Y], Eigen::Vector3d::UnitZ())
      * Eigen::AngleAxisd((*x_ptr)[E_P], Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd((*x_ptr)[E_R], Eigen::Vector3d::UnitX());
    Eigen::Vector3d f_sum_w = cog_rot * f_sum_cog;
    dev_x(V_X) = f_sum_w(0) / hydrus_weight_;
    dev_x(V_Y) = f_sum_w(1) / hydrus_weight_;
    dev_x(V_Z) = f_sum_w(2) / hydrus_weight_ - 9.8;

    /* e_r, e_p, e_y */
    /* d e = R_e * w_b */
    Eigen::Vector3d w((*x_ptr)[W_X], (*x_ptr)[W_Y], (*x_ptr)[W_Z]);
    MatrixXd R_e = MatrixXd::Zero(3, 3);
    R_e << 1, tan((*x_ptr)[E_P]) * sin((*x_ptr)[E_R]), tan((*x_ptr)[E_P]) * cos((*x_ptr)[E_R]),
      0, cos((*x_ptr)[E_R]), -sin((*x_ptr)[E_R]),
      0, sin((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]), cos((*x_ptr)[E_R]) / cos((*x_ptr)[E_P]);
    Eigen::Vector3d d_e; d_e.noalias() = R_e * w;
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
      Eigen::Vector3d wi; wi.noalias() = w + VectorXdTo3d(JW_mat * dq);
      double fi = (*u_ptr)[i];
      mid_result.noalias() +=
        (link_center_pos_cog_vec_[time_id][i] + propeller_pos_cog_offset_).
        cross(Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i] * fi, sin(yaw_vec(i)) * s_tilts_[i] * fi, c_tilts_[i] * fi));
      mid_result.noalias() += Eigen::Vector3d(cos(yaw_vec(i)) * s_tilts_[i], sin(yaw_vec(i)) * s_tilts_[i], c_tilts_[i])
        * fi * M_z_(i); // too small to ignore
      mid_result.noalias() -= I_vec_[time_id][i] * JW_mat * ddq;
      mid_result.noalias() -= wi.cross(VectorXdTo3d(I_vec_[time_id][i] * wi));
      mid_result.noalias() -= I_dt_vec_[time_id][i] * wi;
    }
    Eigen::Matrix3d I_sum = Eigen::Matrix3d::Zero();
    for (int i = 0; i < n_links_; ++i)
      I_sum.noalias() += I_vec_[time_id][i];
    dw.noalias() = I_sum.inverse() * (mid_result);
    for (int i = 0; i < 3; ++i)
      dev_x(W_X + i) = dw(i);

    VectorXd new_x;
    new_x = dev_x / slq_discrete_freq_ + *x_ptr;
    *new_relative_x_ptr = getRelativeState(&new_x);
  }

  bool SlqFiniteDiscreteControlHydrus::feedforwardConverged(){
    double fw_max = 0.0;
    for (int i = 0; i < iteration_times_; ++i){
      for (int j = 0; j < u_size_; ++j){
        if (fw_max < fabs((u_fw_vec_[i])(j)))
          fw_max = fabs((u_fw_vec_[i])(j));
      }
    }
    double fw_converge_threshold = (hydrus_weight_ * 9.78 / 4.0) * 0.1;
    if (fw_max < fw_converge_threshold)
      return true;
    else
      return false;
  }

  void SlqFiniteDiscreteControlHydrus::getHydrusInertialTensor(VectorXd *joint_ptr, int time_id){
    std::vector<Eigen::Matrix3d> cur_I_vec;
    std::vector<Eigen::Matrix3d> cur_I_dt_vec;
    for (int i = 0; i < n_links_; ++i){
      Eigen::Matrix3d cur_I = *I_ptr_;
      Eigen::Vector3d center_pos = link_center_pos_cog_vec_[time_id][i];
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
    I_vec_[time_id] = (cur_I_vec);
    I_dt_vec_[time_id] = (cur_I_dt_vec);
  }

  void SlqFiniteDiscreteControlHydrus::getHydrusLinksCenter(VectorXd *joint_ptr, int time_id){
    std::vector<Eigen::Vector3d> links_center_vec, links_head_vec;
    Eigen::Vector3d cog_local_pos = Eigen::Vector3d::Zero();
    for (int i = 0; i < n_links_; ++i){
      links_center_vec.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
      links_head_vec.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
    }
    links_center_vec[baselink_id_] = Eigen::Vector3d(link_length_ / 2.0, 0, 0);
    cog_local_pos.noalias() += links_center_weight_link_frame_vec_[baselink_id_];
    Eigen::Vector3d prev_link_head(link_length_, 0, 0);
    links_head_vec[baselink_id_ + 1] = prev_link_head;
    double joint_ang = 0.0;

    // only considering 2d hydrus
    for (int i = baselink_id_ + 1; i < n_links_; ++i){
      Eigen::Vector3d link_center = Eigen::Vector3d::Zero();
      joint_ang += (*joint_ptr)(i - 1);
      Eigen::Matrix3d rot;
      rot << cos(joint_ang), -sin(joint_ang), 0,
        sin(joint_ang), cos(joint_ang), 0,
        0, 0, 1;
      link_center.noalias() = prev_link_head + rot * Eigen::Vector3d(link_length_ / 2.0, 0, 0);
      links_center_vec[i] = link_center;
      cog_local_pos.noalias() += prev_link_head * link_weight_vec_[i] + rot * links_center_weight_link_frame_vec_[i];
      prev_link_head.noalias() += rot * Eigen::Vector3d(link_length_, 0, 0);
      if (i < n_links_ - 1)
        links_head_vec[i + 1] = prev_link_head;
    }

    prev_link_head = Eigen::Vector3d(0.0, 0.0, 0.0);
    joint_ang = 0.0;
    for (int i = baselink_id_ - 1; i >= 0; --i){
      Eigen::Vector3d link_center = Eigen::Vector3d::Zero();
      joint_ang -= (*joint_ptr)(i);
      Eigen::Matrix3d rot;
      rot << cos(joint_ang), -sin(joint_ang), 0,
        sin(joint_ang), cos(joint_ang), 0,
        0, 0, 1;
      prev_link_head.noalias() += rot * Eigen::Vector3d(-link_length_, 0, 0);
      links_head_vec[i] = prev_link_head;
      link_center.noalias() = prev_link_head + rot * Eigen::Vector3d(link_length_ / 2.0, 0, 0);
      links_center_vec[i] = link_center;
      cog_local_pos.noalias() += prev_link_head * link_weight_vec_[i] + rot * links_center_weight_link_frame_vec_[i];
    }
    link_center_pos_local_vec_[time_id] = (links_center_vec);
    // cog:= (Sum{i} R_i * (Sum{j} pos_{i, j} * m_{i, j}) + offset_i * (Sum{j} m_{i, j}))
    //      / (Sum{i, j} m_{i, j}) ## i is link id, j is extra model id
    cog_local_pos = cog_local_pos / hydrus_weight_;
    cog_pos_local_vec_[time_id] = (cog_local_pos);

    std::vector<Eigen::Vector3d> link_center_cog_pos_vec;
    for (int i = 0; i < n_links_; ++i)
      link_center_cog_pos_vec.push_back(link_center_pos_local_vec_[time_id][i] - cog_local_pos);
    link_center_pos_cog_vec_[time_id] = (link_center_cog_pos_vec);
  }

  void SlqFiniteDiscreteControlHydrus::getHydrusLinksCenterDerivative(VectorXd *joint_ptr, VectorXd *joint_dt_ptr, int time_id){
    std::vector<Eigen::Vector3d> links_center_dt_vec, links_head_dt_vec;
    Eigen::Vector3d cog_local_vel = Eigen::Vector3d::Zero();
    for (int i = 0; i < n_links_; ++i){
      links_center_dt_vec.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
      links_head_dt_vec.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
    }

    Eigen::Vector3d prev_link_head_dt(0, 0, 0);
    double joint_ang = 0.0;
    double joint_ang_dt = 0.0;
    // only considering 2d hydrus
    for (int i = baselink_id_ + 1; i < n_links_; ++i){
      Eigen::Vector3d link_center_dt = Eigen::Vector3d::Zero();
      joint_ang += (*joint_ptr)(i - 1);
      joint_ang_dt += (*joint_dt_ptr)(i - 1);
      Eigen::Matrix3d rot_dt;
      rot_dt << -sin(joint_ang) * joint_ang_dt, -cos(joint_ang) * joint_ang_dt, 0,
        cos(joint_ang) * joint_ang_dt, -sin(joint_ang) * joint_ang_dt, 0,
        0, 0, 0;
      link_center_dt.noalias() = prev_link_head_dt + rot_dt * Eigen::Vector3d(link_length_ / 2.0, 0, 0);
      links_center_dt_vec[i] = link_center_dt;
      cog_local_vel.noalias() += prev_link_head_dt * link_weight_vec_[i] + rot_dt * links_center_weight_link_frame_vec_[i];
      prev_link_head_dt.noalias() += rot_dt * Eigen::Vector3d(link_length_, 0, 0);
      if (i < n_links_ - 1)
        links_head_dt_vec[i + 1] = prev_link_head_dt;
    }

    prev_link_head_dt = Eigen::Vector3d(0.0, 0.0, 0.0);
    joint_ang = 0.0;
    joint_ang_dt = 0.0;
    for (int i = baselink_id_ - 1; i >= 0; --i){
      Eigen::Vector3d link_center_dt = Eigen::Vector3d::Zero();
      joint_ang -= (*joint_ptr)(i);
      joint_ang_dt -= (*joint_dt_ptr)(i);
      Eigen::Matrix3d rot_dt;
      rot_dt << -sin(joint_ang) * joint_ang_dt, -cos(joint_ang) * joint_ang_dt, 0,
        cos(joint_ang) * joint_ang_dt, -sin(joint_ang) * joint_ang_dt, 0,
        0, 0, 0;
      prev_link_head_dt.noalias() += rot_dt * Eigen::Vector3d(-link_length_, 0, 0);
      links_head_dt_vec[i] = prev_link_head_dt;
      link_center_dt.noalias() = prev_link_head_dt + rot_dt * Eigen::Vector3d(link_length_ / 2.0, 0, 0);
      links_center_dt_vec[i] = link_center_dt;
      cog_local_vel.noalias() += prev_link_head_dt * link_weight_vec_[i] + rot_dt * links_center_weight_link_frame_vec_[i];
    }
    link_center_pos_local_dt_vec_[time_id] = (links_center_dt_vec);
    cog_local_vel = cog_local_vel / hydrus_weight_;
    cog_pos_local_dt_vec_[time_id] = (cog_local_vel);
  }

  MatrixXd SlqFiniteDiscreteControlHydrus::getJacobianW(int id){
    // eg. when baselink is 0: [0 0 0; 0 0 0; 0 0 0], [0 0 0; 0 0 0; 1 0 0],
    // [0 0 0; 0 0 0; 1 1 0], [0 0 0; 0 0 0; 1 1 1 ]

    // eg. when baselink is 1: [0 0 0; 0 0 0; -1 0 0], [0 0 0; 0 0 0; 0 0 0],
    // [0 0 0; 0 0 0; 0 1 0], [0 0 0; 0 0 0; 0 1 1 ]

    // eg. when baselink is 2: [0 0 0; 0 0 0; -1 -1 0], [0 0 0; 0 0 0; 0 -1 0],
    // [0 0 0; 0 0 0; 0 0 0], [0 0 0; 0 0 0; 0 0 1 ]
    MatrixXd JW_mat = MatrixXd::Zero(3, n_links_ - 1);
    if (baselink_id_ == id){}
    else if (id > baselink_id_){
      for (int i = baselink_id_; i < id; ++i)
        JW_mat(2, i) = 1.0;
    }
    else if (id < baselink_id_){
      for (int i = id; i < baselink_id_; ++i)
        JW_mat(2, i) = -1.0;
    }
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
    for (int i = x_size_ - 3; i < x_size_; ++i){
      while (result[i] > PI)
        result[i] -= 2 * PI;
      while (result[i] < -PI)
        result[i] += 2 * PI;
    }
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
    (*H_ptr_).noalias() = (*R_ptr_) + B_ptr_->transpose() * (*P_ptr_) * (*B_ptr_);
    (*G_ptr_).noalias() = B_ptr_->transpose() * (*P_ptr_) * (*A_ptr_);
    (*g_ptr_).noalias() = (*r_ptr_) + B_ptr_->transpose() * (*p_ptr_);
    (*K_ptr_).noalias() = -(H_ptr_->inverse() * (*G_ptr_));
    (*l_ptr_).noalias() = -(H_ptr_->inverse() * (*g_ptr_));

    (*P_ptr_) = A_ptr_->transpose() * (*P_ptr_) * (*A_ptr_);
    (*P_ptr_).noalias() += (*Q_ptr_);
    (*P_ptr_).noalias() += K_ptr_->transpose() * (*H_ptr_) * (*K_ptr_);
    (*P_ptr_).noalias() += K_ptr_->transpose() * (*G_ptr_);
    (*P_ptr_).noalias() += G_ptr_->transpose() * (*K_ptr_);

    (*p_ptr_) = A_ptr_->transpose() * (*p_ptr_);
    (*p_ptr_).noalias() += (*q_ptr_);
    (*p_ptr_).noalias() += K_ptr_->transpose() * (*H_ptr_) * (*l_ptr_);
    (*p_ptr_).noalias() += K_ptr_->transpose() * (*g_ptr_);
    (*p_ptr_).noalias() += G_ptr_->transpose() * (*l_ptr_);
  }

  void SlqFiniteDiscreteControlHydrus::FDLQR(){
    if (verbose_)
      ROS_INFO("[SLQ] LQR init starts.");
    *x_ptr_ = x_vec_[0];
    *u_ptr_ = u_vec_[0];
    *joint_ptr_ = joint_vec_[0];
    updateMatrixAB(0);

    MatrixXd P(x_size_, x_size_);
    P = *Q0_ptr_;
    for (int i = 0; i < iteration_times_; ++i){
      MatrixXd F = MatrixXd::Zero(u_size_, x_size_);
      F.noalias() = ((*R_ptr_) + B_ptr_->transpose() * P * (*B_ptr_)).inverse()
        * (B_ptr_->transpose() * P * (*A_ptr_));
      MatrixXd P_prev = P;
      P.noalias() = A_ptr_->transpose() * P_prev * (*A_ptr_);
      P.noalias() -= (A_ptr_->transpose() * P_prev * (*B_ptr_)) * F;
      P.noalias() += (*Q_ptr_);
      lqr_F_vec_[i] = (F);
    }

    VectorXd x = x_vec_[0];
    for (int i = iteration_times_ - 1; i >= 0; --i){
      VectorXd u; u.noalias() = -lqr_F_vec_[i] * x;
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
    if (verbose_)
      ROS_INFO("[SLQ] LQR init finished");
  }

  double SlqFiniteDiscreteControlHydrus::calculateCostFunction(){
    double cost = 0.0;
    for (int i = 0; i < iteration_times_; ++i){
      // normal cost
      cost += (u_vec_[i].transpose() * (*R_ptr_) * u_vec_[i]
               + x_vec_[i].transpose() * (*Q0_ptr_) * x_vec_[i])(0);
      // waypoint cost
      double cur_time;
      cur_time = double(i) / slq_discrete_freq_;
      VectorXd real_x = getAbsoluteState(&(x_vec_[i]));
      for (int j = 1; j < waypoints_ptr_->size() - 1; ++j){
        MatrixXd W = MatrixXd::Zero(x_size_, x_size_);
        updateWaypointWeightMatrix(cur_time, (*time_ptr_)[j] - (*time_ptr_)[0], &W, j == (waypoints_ptr_->size() - 1));
        VectorXd dx_pt = stateSubtraction(&real_x, &((*waypoints_ptr_)[j]));
        cost += (dx_pt.transpose() * W * dx_pt)(0);
      }
    }
    // final time(tf) cost
    cost += (x_vec_[iteration_times_].transpose() * (*P0_ptr_) * x_vec_[iteration_times_])(0);
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
      return (state_minor * ((id+1) / slq_discrete_freq_ - relative_time)
              + state_max * (relative_time - id / slq_discrete_freq_))
        * slq_discrete_freq_;
    }
  }
}


