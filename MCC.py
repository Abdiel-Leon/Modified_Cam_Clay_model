#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 17:31:54 2023

@author: abdiel
"""
import numpy as np
import matplotlib.pyplot as plt

class MCC:
    
    def __init__(self):
        phi = 30
        phi = phi*3.1416/180
        
        self.pc = 1e4 # preconsolidation pressure
       # p0_dot = Compute_p0_dot()

        
        self.M = 6*np.sin(phi)/(3 + np.sin(phi))
        
    def Draw_yield_surface(self):
        q = 0; p_eff=0
        p = []
        q_MCC= []
        
        while p_eff<self.pc:        
          #solve for q at f=0'
          q = np.sqrt(-self.M*self.M*p_eff*(p_eff-self.pc)) #MCC
          p.append(p_eff)
          q_MCC.append(q)
          p_eff = p_eff + 10
                            
        return p, q_MCC
    
    def Calculate_Elastic_state(self):
        q_sp1 = 0;
        # test deviatoric stress path
        p_sp1= 4000; f = 0
        q_sp1_l =[]; p_sp1_l =[];
        for i in range(5000):
            f = q_sp1*q_sp1 + self.M*self.M*p_sp1*(p_sp1-self.pc) #MCC
            q_sp1 = q_sp1 + 1
            
            p_sp1_l.append(p_sp1)
            q_sp1_l.append(q_sp1)
                    
                            
        return p_sp1_l, q_sp1_l    
    
    
test = MCC()
p,q_MCC=test.Draw_yield_surface()
#p_sp1,q_sp1 = test.Calculate_Elastic_state()
#print(p,q_)

plt.plot(p,q_MCC)
