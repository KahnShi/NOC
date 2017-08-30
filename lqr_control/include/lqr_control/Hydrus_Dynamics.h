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

#ifndef HYDRUS_DYNAMICS_H
#define HYDRUS_DYNAMICS_H

#include <iostream>
/* linear algebra */
#include <math.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Eigenvalues>

using namespace Eigen;

namespace hydrus_dynamics{
  class HydrusDynamics{
  public:
    HydrusDynamics(double l_length, double l_weight_1, double l_weight_2, double l_weight_3, double l_weight_4);
    ~HydrusDynamics();
    double link_length;
    double link_weight_1;
    double link_weight_2;
    double link_weight_3;
    double link_weight_4;
    /* state */
    double px;
    double py;
    double pz;
    double er;
    double ep;
    double ey;
    double d_px;
    double d_py;
    double d_pz;
    double d_er;
    double d_ep;
    double d_ey;
    /* control */
    double q1;
    double q2;
    double q3;
    double f1;
    double f2;
    double f3;
    double f4;
    /* matrix */
    VectorXd *x_ptr_;
    VectorXd *u_ptr_;
    MatrixXd *Ds_ptr_;
    MatrixXd *C_ptr_;
    VectorXd *Bs_ptr_;
    VectorXd *Csds_ptr_;
    VectorXd *gs_ptr_;
    MatrixXd *Ds3_ptr_;
    MatrixXd *Cs3_ptr_;
    VectorXd *dds2_ptr_;
    // Ds_x Bs_x Bs_u Csds_x Csds_dx gs_x Ds3_x Cs3_x Cs3_dx
    MatrixXd *Ds_x;
    MatrixXd *Bs_x;
    MatrixXd *Bs_u;
    MatrixXd *Csds_x;
    MatrixXd *Csds_dx;
    MatrixXd *gs_x;
    MatrixXd *Ds3_x;
    MatrixXd *Cs3_x;
    MatrixXd *Cs3_dx;

    voidcalculateSkewOperation(std::vector<MatrixXd> *S_vec_ptr);
  };
}


#endif
