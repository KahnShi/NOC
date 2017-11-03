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

#include <lqr_control/QuadrotorSimulator.h>
namespace quadrotor_simulator{
  QuadrotorSimulator::QuadrotorSimulator(ros::NodeHandle nh, ros::NodeHandle nhp): nh_(nh), nhp_(nhp){
    nhp_.param("anime_mode", anime_mode_, false);
    // controller_ptr_ = new SlqFiniteDiscreteControlQuadrotor(nh_, nhp_);
    controller_ptr_ = new SlqFiniteDiscreteControlHydrus(nh_, nhp_);
    controller_ptr_->initHydrus();

    pub_traj_path_ = nh_.advertise<nav_msgs::Path>("/lqr_path", 1);
    pub_traj_way_points_ = nh_.advertise<visualization_msgs::MarkerArray>("/end_points_markers", 1);
    pub_anime_joint_state_ = nh_.advertise<sensor_msgs::JointState>("/hydrusx/joint_states", 1);

    oc_iteration_times_ = 0;
    //sleep(1.0);
  }

  void QuadrotorSimulator::initQuadrotorSimulator(std::vector<VectorXd> *waypoints_ptr, std::vector<double> *time_ptr, double controller_freq){
    waypoints_ptr_ = waypoints_ptr;
    time_ptr_ = time_ptr;
    period_ = (*time_ptr)[time_ptr->size() - 1] - (*time_ptr)[0];
    controller_freq_ = controller_freq;
    // controller_ptr_->initLQR(controller_freq_, period_, start_state_ptr_, end_state_ptr_);
    controller_ptr_->initSLQ(controller_freq_, time_ptr_, waypoints_ptr_);
    visualizeTrajectory();
    ROS_INFO("[QuadrotorSimulator] init finished.");
  }

  void QuadrotorSimulator::planOptimalTrajectory(){
    // lqr
    // controller_ptr_->updateAll();
    // std::cout << "[QuadrotorSimulator] updateAll finished\n";

    // slq
    while (1){
      std::cout << "\n[QuadrotorSimulator] Iteration times: " << oc_iteration_times_ << "\n";
      std::cout << "[press 1 to continue, 2 to break, 5 to continue 5 groups, 11 to use anime visualization]\n";
      int id, repeat_times = 1;
      std::cin >> id;
      if (id == 2)
        break;
      else if (id == 5)
        repeat_times = 5;
      else if (id == 11 && anime_mode_){
        std::cout << "[Quadrotor Simulator] Anime mode starts\n";
        animizeTrajectory();
        repeat_times = 0;
      }
      else
        repeat_times = 1;
      for (int i = 0; i < repeat_times; ++i){
        controller_ptr_->iterativeOptimization();
        ++oc_iteration_times_;
        visualizeTrajectory();
      }
    }
  }

  void QuadrotorSimulator::visualizeTrajectory(){
    int n_waypoints = waypoints_ptr_->size();
    traj_.poses.clear();
    visualization_msgs::MarkerArray waypoints_markers;
    visualization_msgs::Marker point_marker;
    point_marker.ns = "waypoints";
    if (anime_mode_)
      point_marker.header.frame_id = std::string("/shi_world");
    else
      point_marker.header.frame_id = std::string("/world");
    point_marker.header.stamp = ros::Time().now();
    point_marker.action = visualization_msgs::Marker::ADD;
    point_marker.type = visualization_msgs::Marker::SPHERE;

    point_marker.id = 0;
    point_marker.pose.position.x = (*waypoints_ptr_)[0](0);
    point_marker.pose.position.y = (*waypoints_ptr_)[0](1);
    point_marker.pose.position.z = (*waypoints_ptr_)[0](2);
    point_marker.pose.orientation.x = 0.0;
    point_marker.pose.orientation.y = 0.0;
    point_marker.pose.orientation.z = 0.0;
    point_marker.pose.orientation.w = 1.0;
    point_marker.scale.x = 0.2;
    point_marker.scale.y = 0.2;
    point_marker.scale.z = 0.2;
    point_marker.color.a = 1;
    point_marker.color.r = 0.0f;
    point_marker.color.g = 1.0f;
    point_marker.color.b = 0.0f;
    waypoints_markers.markers.push_back(point_marker);

    point_marker.id = 1;
    point_marker.pose.position.x = (*waypoints_ptr_)[n_waypoints-1](0);
    point_marker.pose.position.y = (*waypoints_ptr_)[n_waypoints-1](1);
    point_marker.pose.position.z = (*waypoints_ptr_)[n_waypoints-1](2);
    point_marker.color.r = 1.0f;
    point_marker.color.g = 0.0f;
    point_marker.color.b = 0.0f;
    waypoints_markers.markers.push_back(point_marker);

    for (int i = 1; i <= n_waypoints - 2; ++i){
      point_marker.id += 1;
      point_marker.pose.position.x = (*waypoints_ptr_)[i](0);
      point_marker.pose.position.y = (*waypoints_ptr_)[i](1);
      point_marker.pose.position.z = (*waypoints_ptr_)[i](2);
      point_marker.color.r = 1.0f;
      point_marker.color.g = 1.0f;
      point_marker.color.b = 0.0f;
      waypoints_markers.markers.push_back(point_marker);
    }

    pub_traj_way_points_.publish(waypoints_markers);

    if (anime_mode_)
      traj_.header.frame_id = std::string("/shi_world");
    else
      traj_.header.frame_id = std::string("/world");
    traj_.header.stamp = ros::Time().now();
    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header = traj_.header;
    pose_stamped.pose.orientation.x = 0.0f;
    pose_stamped.pose.orientation.y = 0.0f;
    pose_stamped.pose.orientation.z = 0.0f;
    pose_stamped.pose.orientation.w = 1.0f;
    for (int i = 0; i < controller_ptr_->x_vec_.size(); ++i){
      VectorXd result = controller_ptr_->getAbsoluteState(&(controller_ptr_->x_vec_[i]));
      pose_stamped.pose.position.x = (result)(0);
      pose_stamped.pose.position.y = (result)(1);
      pose_stamped.pose.position.z = (result)(2);
      traj_.poses.push_back(pose_stamped);
    }
    pub_traj_path_.publish(traj_);
  }

  void QuadrotorSimulator::animizeTrajectory(){
    // init
    sensor_msgs::JointState joint_state;
    joint_state.name.push_back("joint1");
    joint_state.name.push_back("joint2");
    joint_state.name.push_back("joint3");
    for (int i = 0; i < 3; ++i){
      joint_state.position.push_back(0.0);
      joint_state.velocity.push_back(0.0);
    }
    tf::TransformBroadcaster br;

    for (int i = 0; i < int(controller_freq_ * period_); ++i){
      joint_state.header.seq = i;
      joint_state.header.stamp = ros::Time::now();
      for (int j = 0; j < 3; ++j){
        joint_state.position[j] = controller_ptr_->joint_vec_[i](j);
        joint_state.velocity[j] = controller_ptr_->joint_dt_vec_[i](j);
      }
      pub_anime_joint_state_.publish(joint_state);

      tf::Transform world2cog_transform_;
      geometry_msgs::Pose cog_pose;
      cog_pose.position.x = controller_ptr_->x_vec_[i](P_X) + (*controller_ptr_->xn_ptr_)(P_X);
      cog_pose.position.y = controller_ptr_->x_vec_[i](P_Y) + (*controller_ptr_->xn_ptr_)(P_Y);
      cog_pose.position.z = controller_ptr_->x_vec_[i](P_Z) + (*controller_ptr_->xn_ptr_)(P_Z);
      tf::Quaternion q = tf::createQuaternionFromRPY
        (controller_ptr_->x_vec_[i](E_R) + (*controller_ptr_->xn_ptr_)(E_R),
         controller_ptr_->x_vec_[i](E_P) + (*controller_ptr_->xn_ptr_)(E_P),
         controller_ptr_->x_vec_[i](E_Y) + (*controller_ptr_->xn_ptr_)(E_Y));
      cog_pose.orientation.w = q.w();
      cog_pose.orientation.x = q.x();
      cog_pose.orientation.y = q.y();
      cog_pose.orientation.z = q.z();
      tf::poseMsgToTF(cog_pose, world2cog_transform_);
      br.sendTransform(tf::StampedTransform(world2cog_transform_.inverse(), ros::Time::now(), "cog", "shi_world"));

      ros::Duration(1.0 / controller_freq_).sleep();
    }
  }
}


