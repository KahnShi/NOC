// -*- mode: c++ -*-
/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2015, JSK Lab
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

#ifndef LQR_DISCRETE_CONTROLLER_BASE_H
#define LQR_DISCRETE_CONTROLLER_BASE_H

/* ros */
#include <ros/ros.h>

/* linear algebra */
#include <math.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Eigenvalues>


/* general header file */
#include <iostream>

using namespace Eigen;
namespace lqr_discrete
{
  class LqrDiscreteControlBase{
  public:
    LqrDiscreteControlBase(ros::NodeHandle nh, ros::NodeHandle nhp);
    ~LqrDiscreteControlBase();

    ros::NodeHandle nh_;
    ros::NodeHandle nhp_;

    ros::Timer  control_timer_;

    double control_freq_;
    double end_time_;
    int iteration_times_;
    int u_size_;
    int x_size_;
    MatrixXd *A_ptr_;
    MatrixXd *B_ptr_;
    VectorXd *input_ptr_;
    VectorXd *control_ptr_;
    MatrixXd *Q_ptr_;
    MatrixXd *R_ptr_;
    VectorXd *x0_ptr_;
    VectorXd *xn_ptr_;
    VectorXd *x_ptr_;
    VectorXd *u_ptr_;
    VectorXd *un_ptr_;
    std::vector<MatrixXd *> F_ptr_vec_;
    std::vector<VectorXd *> x_ptr_vec_;
    std::vector<VectorXd *> u_ptr_vec_;
    std::vector<Vector3d> x_vec_;

    void initLQR(double freq, double period, MatrixXd *A, MatrixXd *B, MatrixXd *Q, MatrixXd *R, VectorXd *s0);
  };
}
#endif
