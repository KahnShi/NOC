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

  VectorXd start_state = VectorXd::Zero(13);
  start_state(0) = 5.0;
  start_state(1) = -10.0;
  start_state(2) = 3;
  start_state(Q_W) = 1.0;
  // start_state(Q_X) = 0.707;
  // start_state(Q_Y) = 0.0;
  // start_state(Q_Z) = 0.0;
  VectorXd end_state = VectorXd::Zero(13);
  end_state(0) = 0.0;
  end_state(1) = 0.0;
  end_state(2) = 0.0;
  end_state(Q_W) = 1.0;
  quadrotor_sim_node->initQuadrotorSimulator(&start_state, &end_state, 10.0, 100.0);
  quadrotor_sim_node->planOptimalTrajectory();
  quadrotor_sim_node->visualizeTrajectory();

  ros::spin();
  return 0;
}
