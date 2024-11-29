# -*- coding: utf-8 -*-
"""
Created on Tue May 15 06:17:50 2018

@author: V.D. Pamulaparthy
"""
#demo of trajectory controller that uses nonlinear feedback-linearization for a two link arm
from math import sin,cos
import numpy as np

class trajectory():
      def __init__(self,qd1=np.zeros(100),qd2=np.zeros(100),qd11=np.zeros(100),qd12=np.zeros(100),qd111=np.zeros(100),qd112=np.zeros(100),t=np.zeros(100)):
           self.qd1=qd1
           self.qd2=qd2
           self.qd11=qd11
           self.qd12=qd12
           self.qd111=qd111
           self.qd112=qd112
           self.t=t
           
      def desired_trajectory(self):
           self.t=np.arange(0,10,0.1)
           for j in range(100):
             self.qd1[j]=(0.1*sin(np.pi*self.t[j]))
             self.qd2[j]=(0.1*cos(np.pi*self.t[j]))
             self.qd11[j]=(0.1*np.pi*cos(np.pi*self.t[j]))
             self.qd12[j]=(0.1*np.pi*(-sin(np.pi*self.t[j])))
             self.qd111[j]=(0.1*((np.pi)**2))*(-sin(np.pi*self.t[j]))
             self.qd112[j]=(0.1*((np.pi)**2)*(-cos(np.pi*self.t[j])))
             
             

############################################################################################################################################

class robot_dynamics(trajectory):
      def __init__(self,dt=0.1,q01=0.1,q02=0,q011=0,q012=0,kp=100,kd=20,m1=1,m2=1,a1=1,a2=1,g=9.8):    
            trajectory. __init__(self,qd1=np.zeros(100),qd2=np.zeros(100),qd11=np.zeros(100),qd12=np.zeros(100),qd111=np.zeros(100),qd112=np.zeros(100),t=np.zeros(100))
            self.dt=dt
            self.q01=q01
            self.q02=q02
            self.q011=q011
            self.q012=q012
            self.kp=100
            self.kd=20
            self.m1=m1 
            self.m2=m2
            self.a1=a1
            self.a2=a2
            self.g=g
     
      
      def inertia(self):
        self.M11=((self.m1+self.m2)*(self.a1**2)+self.m2*((self.a2)**2) +2*self.m2*self.a1*self.a2*cos(self.q2[self.i]))
        self.M12=self.m2*((self.a2)**2)+self.m2*self.a1*self.a2*cos(self.q2[self.i])
        self.M22=self.m2*((self.a2)**2)
        self.M=np.matrix([[self.M11, self.M12],[self.M12,self.M22]])
        self.MI=(1/(self.M11*self.M22-(self.M12**2)))*np.matrix([[self.M22, -self.M12],[-self.M12,self.M11]])
    
      def nonlinearities(self):
         self.N1 = -self.m2*self.a1*self.a2*(2*self.q11[self.i]*self.q12[self.i]+(self.q12[self.i]**2))*sin(self.q2[self.i])
         self.N1 =self.N1+(self.m1+self.m2)*self.g*self.a1*cos(self.q1[self.i])+self.m2*self.g*self.a2*cos(self.q1[self.i]+self.q2[self.i])
         self.N2= self.m2*self.a1*self.a2*(self.q11[self.i]**2)*sin(self.q2[self.i])+self.m2*self.g*self.a2*cos(self.q1[self.i]+self.q2[self.i])
         self.N=np.matrix([[self.N1],[self.N2]])
           
         
      def errors(self):
            print(self.i)
            self.e1.append(-self.q1[self.i]+self.qd1[self.i])
            self.e2.append( -self.q2[self.i]+self.qd2[self.i])
            self.e11.append( -self.q11[self.i]+self.qd11[self.i])
            self.e12.append( -self.q12[self.i]+self.qd12[self.i])
            
            
      def dynamics(self):
            self.q1=np.zeros(100)
            self.q2=np.zeros(100)
            self.q11=np.zeros(100)
            self.q12=np.zeros(100)
            self.q1[0]=self.q01
            self.q2[0]=self.q02
            self.q11[0]=self.q011
            self.q12[0]=self.q012
            self.e1=[]
            self.e2=[]
            self.e11=[]
            self.e12=[]
            self.desired_trajectory()
            for self.i in range(99):
                  self.inertia()
                  self.nonlinearities()
                  self.controller()
                  self.q1[self.i+1] = self.q1[self.i]+ self.dt*self.q11[self.i]
                  self.q2[self.i+1] = self.q2[self.i]+ self.dt*self.q12[self.i]
                  self.q11[self.i+1] = self.q11[self.i]+self.dt*((self.MI[0,0]*(-self.N1+self.U[0,0]))+(self.MI[0,1]*(-self.N2+self.U[1,0])))
                  self.q12[self.i+1] = self.q12[self.i]+self.dt*((self.MI[1,0]*(-self.N1+self.U[0,0]))+(self.MI[1,1]*(-self.N2+self.U[1,0])))
            
      def controller(self):
             self.errors()
             self.inertia()
             self.u1=self.qd111[self.i]+self.kp*self.e1[self.i]+self.kd*self.e11[self.i]
             self.u2=self.qd112[self.i]+self.kp*self.e2[self.i]+self.kd*self.e12[self.i]
             self.u=np.matrix([[self.u1],[self.u2]])
             self.U=np.matmul(self.M,self.u)+self.N
             return self.u    
         
