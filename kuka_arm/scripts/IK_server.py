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
        print("T0_1 = ", GetTransEval(T0_1))

        T0_2 = simplify(T0_1*T1_2)
        print("T0_2 = ", GetTransEval(T0_2))

        T0_3 = simplify(T0_2*T2_3)
        print("T0_3 = ", GetTransEval(T0_3))

        T0_4 = simplify(T0_3*T3_4)
        print("T0_4 = ", GetTransEval(T0_4))

        T0_5 = simplify(T0_4*T4_5)
        print("T0_5 = ", GetTransEval(T0_5))

        T0_6 = simplify(T0_5*T5_6)
        print("T0_6 = ", GetTransEval(T0_6))

        T0_G = simplify(T0_6*T6_G)
        print("T0_G = ", GetTransEval(T0_G))

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
            R_ee_x = Matrix([
                          [1,         0,          0],
                          [0, cos(roll), -sin(roll)],
                          [0, sin(roll),  cos(roll)]
                          ])
            R_ee_y = Matrix([
                          [cos(pitch),  0, sin(pitch)],
                          [0,           1,          0],
                          [-sin(pitch), 0, cos(pitch)]
                          ])
            R_ee_z = Matrix([
                         [cos(yaw), -sin(yaw), 0],
                         [sin(yaw),  cos(yaw), 0],
                         [       0,         0, 1]
                         ])
            R_ee = R_ee_z * R_ee_y * R_ee_x

            # Compensate for rotation discrepancy between DH parameters and Gazebo
            # get homogeneous transform for z-axis
            R_z = Matrix([
                        [cos(np.pi), -sin(np.pi), 0, 0],
                        [sin(np.pi),  cos(np.pi), 0, 0],
                        [         0,           0, 1, 0],
                        [         0,           0, 0, 1]
                        ])
            # get homogeneous transform for y-axis
            R_y = Matrix([
                        [ cos(-np.pi/2), 0, sin(-np.pi/2), 0],
                        [             0, 1,             0, 0],
                        [-sin(-np.pi/2), 0, cos(-np.pi/2), 0],
                        [             0, 0,             0, 1]
                        ])
            # get the compensate transform
            R_comp = simplify(R_z * R_y)

            R_ee = R_ee * R_comp

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
	    #
            ###

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

# a helper function to numerically evaluate transforms when all input q is 0
def GetTransEval(T):
    return T.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0, q7: 0})

def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
