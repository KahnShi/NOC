#!/usr/bin/env python

# Copyright (c) 2016, JSK(University of Tokyo)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Open Source Robotics Foundation, Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior
#       written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# Authors: Fan Shi
# Maintainer: Fan Shi <shifan@jsk.imi.i.u-tokyo.ac.jp>

from __future__ import print_function
import time
import sys
import math
import rospy
import tf
import numpy as np
import numpy.matlib
from std_msgs.msg import Float64, Float64MultiArray, MultiArrayDimension
import scipy
import scipy.linalg as linalg
from lqr_control.srv import *
from lqr_control.msg import *

class dareSolver:
    def init(self):
        rospy.init_node('dare_solver_sever' , anonymous=True)
        rospy.Service('dare_solver', Dare, self.__solve_riccati_equation)
        self.__debug = False

    def __solve_riccati_equation(self, req):
        rospy.loginfo("[Dare Solver]Receive and start to solve Riccati equation.")
        x_size = req.A.array.layout.dim[0].size
        u_size = req.B.array.layout.dim[1].size
        R = 50.0 * np.matlib.eye(u_size, dtype=float)
        Q = np.matlib.eye(x_size, dtype=float)
        for i in range(0, 6):
            Q[i, i] = 10.0
        A = np.matlib.zeros((x_size, x_size))
        B = np.matlib.zeros((x_size, u_size))

        for i in range(0, x_size):
            for j in range(0, x_size):
                A[i, j] = req.A.array.data[i * x_size + j]

        for i in range(0, x_size):
            print()
            for j in range(0, u_size):
                B[i, j] = req.B.array.data[i * u_size + j]
                print(B[i, j], ', ', end='')

        P = linalg.solve_discrete_are(A, B, Q, R)

        if self.__debug:
            print("Output Matrix A:")
            for i in range(0, x_size):
                for j in range(0, x_size):
                    print(A[i, j], ', ', end='')
                print('')
            print("Output Matrix B:")
            for i in range(0, x_size):
                for j in range(0, u_size):
                    print(B[i, j], ', ', end='')
                print('')
            print("Output Matrix Q:")
            for i in range(0, x_size):
                for j in range(0, x_size):
                    print(Q[i, j], ', ', end='')
                print('')
            print("Output Matrix R:")
            for i in range(0, u_size):
                for j in range(0, u_size):
                    print(R[i, j], ', ', end='')
                print('')
            print("Output Matrix P:")
            for i in range(0, x_size):
                for j in range(0, x_size):
                    print(P[i, j], ', ', end='')
                print('')

        ## publish P matrix
        dim0 = MultiArrayDimension()
        dim0.label = "height"
        dim0.size = x_size
        dim0.stride = x_size * x_size
        dim1 = MultiArrayDimension()
        dim1.label = "width"
        dim1.size = x_size
        dim1.stride = x_size
        P_array = float64Array()
        P_array.array.layout.dim.append(dim0)
        P_array.array.layout.dim.append(dim1)
        P_array.array.layout.data_offset = 0
        for i in range(0, x_size):
            for j in range(0, x_size):
                P_array.array.data.append(P[i, j])

        rospy.loginfo("[Dare Solver]Riccati matrix is sent.")
        return DareResponse(P_array)

if __name__ == '__main__':
    try:
        dare_solver = dareSolver()
        dare_solver.init()
        rospy.spin()
    except rospy.ROSInterruptException:
        pass
