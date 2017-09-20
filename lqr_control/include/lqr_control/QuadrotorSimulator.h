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

#ifndef QUADROTOR_SIMULATOR_H
#define QUADROTOR_SIMULATOR_H

// #include <lqr_control/SlqFD_Quadrotor.h>
#include <lqr_control/SlqFD_Hydrus.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <std_msgs/Empty.h>
#include <sensor_msgs/JointState.h>

using namespace lqr_discrete;
namespace quadrotor_simulator{
  class QuadrotorSimulator{
  public:
    QuadrotorSimulator(ros::NodeHandle nh, ros::NodeHandle nhp);
    ~QuadrotorSimulator(){};
    ros::NodeHandle nh_;
    ros::NodeHandle nhp_;
    bool anime_mode_;
    // SlqFiniteDiscreteControlQuadrotor *controller_ptr_;
    SlqFiniteDiscreteControlHydrus *controller1_ptr_;
    SlqFiniteDiscreteControlHydrus *controller2_ptr_;
    SlqFiniteDiscreteControlHydrus *controller3_ptr_;
    std::vector<SlqFiniteDiscreteControlHydrus *> controller_ptr_vec_;
    std::vector<VectorXd> *waypoints_ptr_;
    std::vector<double> *time_ptr_;
    double period_;
    double controller_freq_;
    int oc_iteration_times_;

    nav_msgs::Path traj_;
    ros::Publisher pub_traj_path_;
    ros::Publisher pub_traj_way_points_;
    ros::Publisher pub_anime_joint_state_;

    void initQuadrotorSimulator(std::vector<VectorXd> *waypoints_ptr, std::vector<double> *time_ptr, double controller_freq);
    void planOptimalTrajectory();
    void visualizeTrajectory(int id);
    void animizeTrajectory(int id);
  };
}


#endif
