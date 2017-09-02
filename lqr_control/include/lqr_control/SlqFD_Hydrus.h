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

#ifndef SLQ_FINITE_DISCRETE_CONTROLLER_HYDRUS_H
#define SLQ_FINITE_DISCRETE_CONTROLLER_QUADROTOR_H

#include <lqr_control/LqrFD_Quadrotor.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/MultiArrayDimension.h>
#include <unistd.h>
#include <lqr_control/Dare.h>
#include <lqr_control/float64Array.h>

namespace lqr_discrete{
  class SlqFiniteDiscreteControlHydrus: public LqrDiscreteControlBase{
  public:
    SlqFiniteDiscreteControlHydrus(ros::NodeHandle nh, ros::NodeHandle nhp):
      LqrDiscreteControlBase(nh, nhp){};

    /* Hydrus */
    double link_length_;
    int n_links_;
    std::vector<double> link_weight_vec_;
    VectorXd *joint_ptr_;
    std::vector<VectorXd> joint_vec_;
    double hydrus_weight_;
    MatrixXd *I_ptr_;
    std::vector<std::vector<MatrixXd> > I_vec_;
    std::vector<std::vector<Vector3d> > link_center_pos_local_vec_;
    MatrixXd *M_para_ptr_;
    bool debug_;
    double uav_rotor_thrust_min_;
    double uav_rotor_thrust_max_;
    /* slq */
    MatrixXd *H_ptr_;
    MatrixXd *Riccati_P_ptr_;
    MatrixXd *P_ptr_;
    VectorXd *p_ptr_;
    MatrixXd *K_ptr_;
    MatrixXd *G_ptr_;
    VectorXd *g_ptr_;
    VectorXd *l_ptr_;
    VectorXd *r_ptr_;
    VectorXd *q_ptr_;
    MatrixXd *Q0_ptr_;
    double alpha_;
    int line_search_steps_;
    std::vector<Vector4d> u_fw_vec_;
    std::vector<Vector4d> u_fb_vec_;
    std::vector<MatrixXd> K_vec_;
    std::vector<VectorXd> *waypoints_ptr_;
    std::vector<double> *time_ptr_;
    LqrFiniteDiscreteControlQuadrotor *lqr_controller_ptr_;

    /* Ros service */
    ros::ServiceClient dare_client_;

    void initSLQ(double freq, std::vector<double> *time_ptr, std::vector<VectorXd> *waypoints_ptr);
    void updateMatrixA(VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr);
    void updateMatrixB(VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr);
    void updateMatrixAB(VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr);
    void getRiccatiH();
    void iterativeOptimization();
    void updateNewState(VectorXd *new_x_ptr, VectorXd *x_ptr, VectorXd *u_ptr, VectorXd *joint_ptr);
    bool feedforwardConverged();
    VectorXd stateAddition(VectorXd *x1_ptr, VectorXd *x2_ptr);
    VectorXd stateSubtraction(VectorXd *x1_ptr, VectorXd *x2_ptr);
    VectorXd getAbsoluteState(VectorXd *relative_x_ptr);
    VectorXd getRelativeState(VectorXd *absolute_x_ptr);
    void updateWaypointWeightMatrix(double time, double end_time, MatrixXd *W_ptr, bool goal_flag);
    void updateSLQEquations();
    void FDLQR();
    void checkControlInputFeasible(VectorXd *u);
    VectorXd getCurrentJoint(double time, int order = 0);
    void getHydrusLinksCenter(VectorXd *joint_ptr);
    void getHydrusInertialTensor(VectorXd *joint_ptr, int time_id);
    void printStateInfo(VectorXd *x, int id);
    void printControlInfo(VectorXd *u, int id);
    void printMatrixAB();
  };
}


#endif
