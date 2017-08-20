#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *

def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        ### Your FK code here
        # Create symbols
        # tricky part: for q and d, it's 1 ~ i
        #              for a and alpha is't 0 ~ i-1
        # qs
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')

        # ds
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')

        # as
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')

        # alphas
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

    	# Create Modified DH parameters
        s = {
            alpha0:      0, a0:      0, d1:  (0.42 + 0.33), q1: q1,
            alpha1: -pi/2., a1:   0.35, d2:              0, q2: q2 - pi/2,
            alpha2:      0, a2:   1.25, d3:              0, q3: q3,
            alpha3: -pi/2., a3: -0.054, d4:  (0.96 + 0.54), q4: q4,
            alpha4:  pi/2., a4:      0, d5:              0, q5: q5,
            alpha5: -pi/2., a5:      0, d6:              0, q6: q6,
            alpha6:      0, a6:      0, d7: (0.193 + 0.11), q7: 0
        }

    	# Define Modified DH Transformation matrix
        T0_1 = GetHTFromDH(q1, d1, a0, alpha0)
        T0_1 = T0_1.subs(s)

        T1_2 = GetHTFromDH(q2, d2, a1, alpha1)
        T1_2 = T1_2.subs(s)

        T2_3 = GetHTFromDH(q3, d3, a2, alpha2)
        T2_3 = T2_3.subs(s)

        T3_4 = GetHTFromDH(q4, d4, a3, alpha3)
        T3_4 = T3_4.subs(s)

        T4_5 = GetHTFromDH(q5, d5, a4, alpha4)
        T4_5 = T4_5.subs(s)

        T5_6 = GetHTFromDH(q6, d6, a5, alpha5)
        T5_6 = T5_6.subs(s)

        T6_G = GetHTFromDH(q7, d7, a6, alpha6)
        T6_G = T6_G.subs(s)

    	# Create individual transformation matrices
        #print("T0_1 = ", GetTransEval(T0_1))

        #print("T0_2 = ", GetTransEval(T0_2))
        #T0_2 = T0_1*T1_2

        #print("T0_3 = ", GetTransEval(T0_3))
        #T0_3 = T0_2*T2_3

        #print("T0_4 = ", GetTransEval(T0_4))
        #T0_4 = T0_3*T3_4

        #print("T0_5 = ", GetTransEval(T0_5))
        #T0_5 = T0_4*T4_5

        #print("T0_6 = ", GetTransEval(T0_6))
        #T0_6 = T0_5*T5_6

        #print("T0_G = ", GetTransEval(T0_G))

        # tan(b0) = a3/d4
        ang_b0 = atan2(0.054, (0.96 + 0.54))

        # Compensate for rotation discrepancy between DH parameters and Gazebo
        # get transform for z-axis
        R_z = Matrix([
                    [cos(pi), -sin(pi), 0],
                    [sin(pi),  cos(pi), 0],
                    [      0,        0, 1]
                    ])
        # get transform for y-axis
        R_y = Matrix([
                    [ cos(-pi/2), 0, sin(-pi/2)],
                    [          0, 1,          0],
                    [-sin(-pi/2), 0, cos(-pi/2)]
                    ])
        # get the compensate transform
        R_comp = simplify(R_z * R_y)

        # Construct rot for ee symbolically
        r, p, y = symbols('r p y')

        R_ee_x = Matrix([
                      [1,         0,          0],
                      [0, cos(r), -sin(r)],
                      [0, sin(r),  cos(r)]
                      ])
        R_ee_y = Matrix([
                      [cos(p),  0, sin(p)],
                      [0,           1,          0],
                      [-sin(p), 0, cos(p)]
                      ])
        R_ee_z = Matrix([
                     [cos(y), -sin(y), 0],
                     [sin(y),  cos(y), 0],
                     [       0,         0, 1]
                     ])
        R_ee = R_ee_z * R_ee_y * R_ee_x
        R_ee = R_ee * R_comp

        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

	    # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            ### Your IK code here
            # Construct the rotation matrix of ee from the roll pitch yaw and pos data
            R_ee = R_ee.subs({'r': roll, 'p': pitch, 'y': yaw})

            # compute wc position
            EE = Matrix([
                        [px],
                        [py],
                        [pz]
                        ])
            # according to the formula from IK lesson
            WC = EE - d7*R_ee[:,2]
            # compute wc
            WC = WC.subs(s)

	        # Calculate joint angles using Geometric IK method
	        # calculate theta1 from wc pos directly
            theta1 = atan2(WC[1], WC[0])

            # calculate theta2
            A = (0.96 + 0.54)  # d4
            C = 1.25           # a2

            # calculate B
            wc_x = WC[0]
            wc_y = WC[1]
            wc_z = WC[2]

            b0 = sqrt(wc_x**2 + wc_y**2)
            b1 = b0 - 0.35 # consider the offset of a1
            b2 = wc_z - 0.75  # the offset of d1
            B = sqrt(b1**2 + b2**2)

            # got all the sides, compute the angles
            cos_ang_a = (B**2 + C**2 - A**2) / (2*B*C)
            cos_ang_b = (A**2 + C**2 - B**2) / (2*A*C)
            cos_ang_c = (A**2 + B**2 - C**2) / (2*A*B)

            ang_a = acos(cos_ang_a)
            ang_b = acos(cos_ang_b)
            ang_c = acos(cos_ang_c)

            # got all the angles, compute theta2 and theta3 from them
            # tan(a0) = b2/b1
            ang_a0 = atan2(b2, b1)
            theta2 = pi/2 - ang_a - ang_a0

            # consider the small slope from joint 3 to joint 4
            # ang_b0 is a constant that can be moved out of the loop to save computation time
            theta3 = pi/2 - (ang_b + ang_b0)

            # construct rot from joint 0 to 3
            R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]

            # evaluate it
            R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

            # get rot from joint 3 to 6
            R3_6 = R0_3.inv("LU") * R_ee

            # get theta4, theta5, theta6
            # same code taken from Lesson 2-8
            # define r11, r12, r13
            r11 = R3_6.row(0)[0]
            r12 = R3_6.row(0)[1]
            r13 = R3_6.row(0)[2]

            # define r21, r22, r23
            r21 = R3_6.row(1)[0]
            r22 = R3_6.row(1)[1]
            r23 = R3_6.row(1)[2]

            # define r31, r32, r33
            r31 = R3_6.row(2)[0]
            r32 = R3_6.row(2)[1]
            r33 = R3_6.row(2)[2]

            theta4 = atan2(r21, r11) # rotation about Z-axis
            theta5 = atan2(-r31, sqrt(r11*r11 + r21*r21)) # rotation about Y-axis
            theta6 = atan2(r32, r33) # rotation about X-axis

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)

# a helper function to generate homogeneous transforms from DH parameters
def GetHTFromDH(q, d, a, alpha):
    # q, d, a, alpha are DH parameters
    T = Matrix([
                [           cos(q),           -sin(q),           0,             a],
                [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                [                0,                 0,           0,             1]
               ])
    return T

def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
