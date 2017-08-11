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

#ifndef SLQ_FINITE_DISCRETE_CONTROLLER_QUADROTOR_H
#define SLQ_FINITE_DISCRETE_CONTROLLER_QUADROTOR_H

#include <lqr_control/LqrFD_Quadrotor.h>
#include <lqr_control/LqrID_Quadrotor.h>
namespace lqr_discrete{
  class SlqFiniteDiscreteControlQuadrotor: public LqrDiscreteControlBase{
  public:
    SlqFiniteDiscreteControlQuadrotor(ros::NodeHandle nh, ros::NodeHandle nhp):
      LqrDiscreteControlBase(nh, nhp){};
    MatrixXd *I_ptr_;
    MatrixXd *M_para_ptr_;
    double uav_mass_;
    bool debug_;
    /* slq */
    MatrixXd *H_ptr_;
    MatrixXd *P_ptr_;
    VectorXd *p_ptr_;
    MatrixXd *K_ptr_;
    MatrixXd *G_ptr_;
    VectorXd *g_ptr_;
    VectorXd *l_ptr_;
    VectorXd *r_ptr_;
    double alpha_;
    std::vector<MatrixXd> F_vec_;
    LqrFiniteDiscreteControlQuadrotor *lqr_controller_ptr_;
    void initSLQ(double freq, double period, VectorXd *x0, VectorXd *xn);
    void updateMatrixA();
    void updateMatrixB();
    void updateMatrixAB();
    void updateAll();
    void getRicattiH();
    void iterativeOptimization();
    void updateNewState(VectorXd *new_x_ptr);
  };
}


#endif
