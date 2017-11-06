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
#include <geometry_msgs/Pose.h>
#include <tf/transform_broadcaster.h>

namespace lqr_discrete{
#define PI 3.141592653
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
    std::vector<VectorXd> joint_dt_vec_;
    std::vector<VectorXd> joint_ddt_vec_;
    double hydrus_weight_;
    VectorXd M_z_;
    Eigen::Matrix3d *I_ptr_;
    std::vector<std::vector<Eigen::Matrix3d> > I_vec_;
    std::vector<std::vector<Eigen::Matrix3d> > I_dt_vec_;
    std::vector<std::vector<Eigen::Vector3d> > link_center_pos_local_vec_;
    std::vector<std::vector<Eigen::Vector3d> > link_center_pos_local_dt_vec_;
    std::vector<Eigen::Vector3d> cog_pos_local_vec_;
    std::vector<Eigen::Vector3d> cog_pos_local_dt_vec_;
    bool debug_;
    double uav_rotor_thrust_min_;
    double uav_rotor_thrust_max_;
    bool transform_movement_flag_;
    /* slq */
    bool not_first_slq_flag_;
    MatrixXd *H_ptr_;
    MatrixXd *Riccati_P_ptr_;
    MatrixXd *IDlqr_F_ptr_;
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
    double alpha_candidate_;
    int line_search_steps_;
    std::vector<MatrixXd> lqr_F_vec_;
    std::vector<Vector4d> u_fw_vec_;
    std::vector<Vector4d> u_fb_vec_;
    std::vector<VectorXd> un_vec_;
    std::vector<MatrixXd> K_vec_;
    std::vector<VectorXd> *waypoints_ptr_;
    std::vector<double> *time_ptr_;
    LqrFiniteDiscreteControlQuadrotor *lqr_controller_ptr_;

    /* infinite state */
    VectorXd stable_u_last_;
    VectorXd xn_last_;
    bool infinite_feedback_update_flag_;

    /* Ros service */
    ros::ServiceClient dare_client_;

    void initHydrus();
    void initSLQ(double freq, std::vector<double> *time_ptr, std::vector<VectorXd> *waypoints_ptr);
    void updateMatrixA(int time_id);
    void updateMatrixB(int time_id);
    void updateMatrixAB(int time_id);
    void getRiccatiH();
    void iterativeOptimization();
    void updateNewState(VectorXd *new_x_ptr, VectorXd *x_ptr, VectorXd *u_ptr, int time_id);
    bool feedforwardConverged();
    VectorXd infiniteFeedbackControl(VectorXd *cur_real_x_ptr);
    VectorXd highFrequencyFeedbackControl(double relative_time, VectorXd *cur_real_x_ptr);
    VectorXd highFrequencyLQRFeedbackControl(double relative_time, VectorXd *cur_real_x_ptr);
    VectorXd stateAddition(VectorXd *x1_ptr, VectorXd *x2_ptr);
    VectorXd stateSubtraction(VectorXd *x1_ptr, VectorXd *x2_ptr);
    VectorXd getAbsoluteState(VectorXd *relative_x_ptr);
    VectorXd getRelativeState(VectorXd *absolute_x_ptr);
    void updateWaypointWeightMatrix(double time, double end_time, MatrixXd *W_ptr, bool goal_flag);
    void updateSLQEquations();
    void FDLQR();
    void checkControlInputFeasible(VectorXd *u, int time_id);
    VectorXd getCurrentJoint(double time, int order = 0);
    Eigen::Matrix3d getCurrentRotationMatrix(Eigen::Vector3d euler_angle, int order = 0);
    void getHydrusLinksCenter(VectorXd *joint_ptr);
    void getHydrusLinksCenterDerivative(VectorXd *joint_ptr, VectorXd *joint_dt_ptr);
    void getHydrusInertialTensor(VectorXd *joint_ptr, int time_id);
    void updateHydrusCogPosition(int time_id);
    void updateHydrusCogPositionDerivative(int time_id);
    MatrixXd getJacobianW(int id);
    VectorXd getStableThrust(int id);
    double calculateCostFunction();
    Eigen::Vector3d VectorXdTo3d(VectorXd vec);
    MatrixXd S_operation(VectorXd vec);
    void printStateInfo(VectorXd *x, int id);
    void printControlInfo(VectorXd *u, int id);
    void printMatrixAB();
    VectorXd getCurrentIdealPosition(double relative_time);
  };
}


#endif
