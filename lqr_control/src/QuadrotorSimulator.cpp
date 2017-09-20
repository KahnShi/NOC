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
    controller1_ptr_ = new SlqFiniteDiscreteControlHydrus(nh_, nhp_);
    //controller1_ptr_->initHydrus();
    controller2_ptr_ = new SlqFiniteDiscreteControlHydrus(nh_, nhp_);
    //controller2_ptr_->initHydrus();
    controller3_ptr_ = new SlqFiniteDiscreteControlHydrus(nh_, nhp_);
    //controller3_ptr_->initHydrus();
    controller_ptr_vec_.push_back(controller1_ptr_);
    controller_ptr_vec_.push_back(controller2_ptr_);
    controller_ptr_vec_.push_back(controller3_ptr_);

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
    ROS_INFO("[QuadrotorSimulator] init finished.");
  }

  void QuadrotorSimulator::planOptimalTrajectory(){
    std::vector<std::vector<VectorXd> > waypoints_vec;
    for (int id = 0; id < waypoints_ptr_->size() - 1; ++id){
      std::vector<VectorXd> waypoints;
      waypoints.push_back((*waypoints_ptr_)[id]);
      waypoints.push_back((*waypoints_ptr_)[id+1]);
      waypoints_vec.push_back(waypoints);
    }
    for (int id = 0; id < waypoints_ptr_->size() - 1; ++id){
      controller_ptr_vec_[id]->initHydrus();
      controller_ptr_vec_[id]->plan_traj_id_ = id;
      std::vector<double> time_vec;
      time_vec.push_back(0);
      time_vec.push_back((*time_ptr_)[id+1] - (*time_ptr_)[id]);

      controller_ptr_vec_[id]->initSLQ(controller_freq_, &time_vec, &(waypoints_vec[id]));
      for (int i = 0; i < 4; ++i){
        controller_ptr_vec_[id]->iterativeOptimization();
        // todo
        visualizeTrajectory(id);
      }
      ROS_WARN("Trajectory optimization finished.\n");
      if (id < waypoints_ptr_->size() - 2){
        waypoints_vec[id+1][0] = controller_ptr_vec_[id]->getFinalAbsoluteState();
      }
    }
  }

  void QuadrotorSimulator::visualizeTrajectory(int id){
    int n_waypoints = waypoints_ptr_->size();
    traj_.poses.clear();
    visualization_msgs::MarkerArray waypoints_markers;
    visualization_msgs::Marker point_marker, wall_marker;
    point_marker.ns = "waypoints";
    if (anime_mode_)
      point_marker.header.frame_id = std::string("/shi_world");
    else
      point_marker.header.frame_id = std::string("/world");
    point_marker.header.stamp = ros::Time().now();
    point_marker.action = visualization_msgs::Marker::ADD;
    point_marker.type = visualization_msgs::Marker::SPHERE;

    wall_marker = point_marker;
    wall_marker.type = visualization_msgs::Marker::CUBE;

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

    wall_marker.id = point_marker.id + 1;
    wall_marker.pose.position.x = 6.1;
    wall_marker.pose.position.y = 4.0;
    wall_marker.pose.position.z = 0.0;
    wall_marker.pose.orientation.x = 0.0;
    wall_marker.pose.orientation.y = 0.0;
    wall_marker.pose.orientation.z = 0.0;
    wall_marker.pose.orientation.w = 1.0;
    wall_marker.scale.x = 10;
    wall_marker.scale.y = 0.2;
    wall_marker.scale.z = 5;
    wall_marker.color.a = 1;
    wall_marker.color.r = 0.6f;
    wall_marker.color.g = 0.45f;
    wall_marker.color.b = 0.15f;
    waypoints_markers.markers.push_back(wall_marker);

    wall_marker.id += 1;
    wall_marker.pose.position.x = -5.0;
    waypoints_markers.markers.push_back(wall_marker);

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
    for (int i = 0; i < controller_ptr_vec_[id]->x_vec_.size(); ++i){
      VectorXd result = controller_ptr_vec_[id]->getAbsoluteState(&(controller_ptr_vec_[id]->x_vec_[i]));
      pose_stamped.pose.position.x = (result)(0);
      pose_stamped.pose.position.y = (result)(1);
      pose_stamped.pose.position.z = (result)(2);
      traj_.poses.push_back(pose_stamped);
    }
    pub_traj_path_.publish(traj_);
  }

  void QuadrotorSimulator::animizeTrajectory(int id){
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
        joint_state.position[j] = controller_ptr_vec_[id]->joint_vec_[i](j);
        joint_state.velocity[j] = controller_ptr_vec_[id]->joint_dt_vec_[i](j);
      }
      pub_anime_joint_state_.publish(joint_state);

      tf::Transform world2cog_transform_;
      geometry_msgs::Pose cog_pose;
      cog_pose.position.x = controller_ptr_vec_[id]->x_vec_[i](P_X) + (*controller_ptr_vec_[id]->xn_ptr_)(P_X);
      cog_pose.position.y = controller_ptr_vec_[id]->x_vec_[i](P_Y) + (*controller_ptr_vec_[id]->xn_ptr_)(P_Y);
      cog_pose.position.z = controller_ptr_vec_[id]->x_vec_[i](P_Z) + (*controller_ptr_vec_[id]->xn_ptr_)(P_Z);
      tf::Quaternion q = tf::createQuaternionFromRPY
        (controller_ptr_vec_[id]->x_vec_[i](E_R) + (*controller_ptr_vec_[id]->xn_ptr_)(E_R),
         controller_ptr_vec_[id]->x_vec_[i](E_P) + (*controller_ptr_vec_[id]->xn_ptr_)(E_P),
         controller_ptr_vec_[id]->x_vec_[i](E_Y) + (*controller_ptr_vec_[id]->xn_ptr_)(E_Y));
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


