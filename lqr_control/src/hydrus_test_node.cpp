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
using namespace lqr_discrete;

int main(int argc, char **argv)
{
  ros::init(argc, argv, "quadrotor_simulator");
  ros::NodeHandle nh;
  ros::NodeHandle nh_private("~");
  SlqFiniteDiscreteControlHydrus* hydrus_node = new SlqFiniteDiscreteControlHydrus(nh, nh_private);

  std::vector<VectorXd> way_pts_vec;
  std::vector<double> period_vec;
  VectorXd start_state = VectorXd::Zero(12);
  start_state(0) = 10.0;
  start_state(1) = 5.0;
  start_state(2) = 3.0;
  // start_state(E_R) = 0.0;
  // start_state(E_P) = 0.0;
  // start_state(E_Y) = 0.0;
  //start_state(Q_W) = 1.0;
  way_pts_vec.push_back(start_state);
  period_vec.push_back(0.0);

  VectorXd mid_state = VectorXd::Zero(12);
  mid_state(0) = 8.0;
  mid_state(1) = 0.0;
  mid_state(2) = 2.0;
  mid_state(3) = -2.0;
  mid_state(4) = -2.0;
  mid_state(5) = -1.0;
  way_pts_vec.push_back(mid_state);
  period_vec.push_back(3.0);

  VectorXd end_state = VectorXd::Zero(12);
  end_state(0) = 5.0;
  end_state(1) = -5.0;
  end_state(2) = 0.0;
  // end_state(E_R) = 3.14 / 6.0;
  // end_state(E_P) = 3.14 / 6.0;
  // end_state(E_P) = 3.14 / 6.0;
  // end_state(E_R) = 3.14 / 6.0;
  //end_state(Q_W) = 1.0;
  way_pts_vec.push_back(end_state);
  period_vec.push_back(6.0);

  hydrusCmdTask task_descriptor;
  task_descriptor.period = 10000.0;
  task_descriptor.end_time = ros::Time::now().toSec() + 10000.0;
  hydrus_node->initSLQ(100.0, &period_vec, &way_pts_vec, task_descriptor);

  VectorXd x = VectorXd::Zero(12);
  VectorXd q = VectorXd::Zero(3);
  q(0) = 3.1416 / 6.0;
  // hydrus_node->updateHydrusLinks(&x, &q);

  ros::spinOnce();
  return 0;
}
