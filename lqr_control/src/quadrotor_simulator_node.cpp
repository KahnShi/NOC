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
using namespace quadrotor_simulator;

int main(int argc, char **argv)
{
  ros::init(argc, argv, "quadrotor_simulator");
  ros::NodeHandle nh;
  ros::NodeHandle nh_private("~");
  QuadrotorSimulator* quadrotor_sim_node = new QuadrotorSimulator(nh, nh_private);
  sleep(1.0);

  std::vector<VectorXd> way_pts_vec;
  std::vector<double> period_vec;
  VectorXd start_state = VectorXd::Zero(12);
  nh_private.param("start_x", start_state(0), 10.0);
  nh_private.param("start_y", start_state(1), 5.0);
  nh_private.param("start_z", start_state(2), 3.0);
  double start_time;
  nh_private.param("start_time", start_time, 0.0);
  start_state(E_Y) = 1.0;
  // start_state(V_X) = 1.0;
  // start_state(V_Y) = -1.0;
  // start_state(V_Z) = 1.0;
  // start_state(0) = 10.0;
  // start_state(1) = 5.0;
  // start_state(2) = 3.0;
  // start_state(E_R) = 0.0;
  // start_state(E_P) = 0.0;
  //start_state(E_Y) = -1.0;
  //start_state(Q_W) = 1.0;
  way_pts_vec.push_back(start_state);
  period_vec.push_back(start_time);

  // VectorXd mid_state = VectorXd::Zero(12);
  // mid_state(0) = 8.0;
  // mid_state(1) = 0.0;
  // mid_state(2) = 2.0;
  // mid_state(3) = -2.0;
  // mid_state(4) = -2.0;
  // mid_state(5) = -1.0;
  // way_pts_vec.push_back(mid_state);
  // period_vec.push_back(3.0);

  VectorXd end_state = VectorXd::Zero(12);
  nh_private.param("end_x", end_state(0), 12.0);
  nh_private.param("end_y", end_state(1), 3.0);
  nh_private.param("end_z", end_state(2), 6.0);
  double end_time;
  nh_private.param("end_time", end_time, 4.0);
  end_state(E_Y) = 2.0;
  // end_state(V_X) = 1.0;
  // end_state(V_Y) = -1.0;
  // end_state(V_Z) = 1.0;
  // end_state(0) = 12.0;
  // end_state(1) = 3.0;
  // end_state(2) = 6.0;
  //end_state(E_R) = 3.14 / 18.0;
  //end_state(E_P) = 3.14 / 18.0;
  // end_state(E_P) = 3.14 / 6.0;
  // end_state(E_R) = 3.14 / 6.0;
  //end_state(Q_W) = 1.0;
  way_pts_vec.push_back(end_state);
  period_vec.push_back(end_time);

  quadrotor_sim_node->initQuadrotorSimulator(&way_pts_vec, &period_vec, 20.0);
  quadrotor_sim_node->planOptimalTrajectory();
  quadrotor_sim_node->visualizeTrajectory();

  ros::spin();
  return 0;
}
