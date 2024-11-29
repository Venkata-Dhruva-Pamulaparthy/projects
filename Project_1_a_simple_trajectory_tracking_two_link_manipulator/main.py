# -*- coding: utf-8 -*-
"""

@author: V.D. Pamulaparthy
"""
#demo of trajectory controller that uses nonlinear feedback-linearization for a two link arm
from robot_arm2 import robot_dynamics 
from matplotlib import pyplot as plt
import numpy as np
dyn=robot_dynamics()
dyn.dynamics()
t=np.arange(0,100)
t1=np.arange(0,99)
plt.figure()
plt.title('desired trajectory tracking: joint angle 2 ')
plt.plot(t,dyn.qd1,label='desired trajectory')
plt.plot(t,dyn.q1,label='joint trajectory')
plt.ylabel('joint angle q2,desired joint angle qd2')
plt.xlabel('time t')
plt.grid()
plt.legend()
plt.figure()
plt.title('desired trajectory tracking: joint angle 1')
plt.plot(t,dyn.qd2,label='desired trajector')
plt.plot(t,dyn.q2,label='joint trajectory')
plt.legend()
plt.ylabel('joint angle q1,desired joint angle qd1')
plt.xlabel('time t')
plt.grid()
plt.legend()
plt.figure()
plt.grid()
plt.title('tracking errors for joint angle q1')
plt.plot(t1,dyn.e1,label='position tracking errors')
plt.plot(t1,dyn.e11,label='velocity tracking errors')
plt.ylabel('tracking errors')
plt.xlabel('time t')
plt.legend()
plt.figure()
plt.grid()
plt.title('tracking errors for joint angle q2')
plt.plot(t1,dyn.e2, label='position tracking errors')
plt.plot(t1,dyn.e12, label='velocity tracking errors')
plt.ylabel('tracking errors')
plt.xlabel('time t')
plt.legend()
