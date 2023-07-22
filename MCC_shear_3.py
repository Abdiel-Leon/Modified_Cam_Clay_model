#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 17:31:54 2023
Need debugging and further improvements!
@author: abdiel
"""
import numpy as np
import matplotlib.pyplot as plt
import math

class MCC:
    
    def __init__(self):
        phi = 30
        phi = phi*3.1416/180
        
        self.pc = 2.0e5 # preconsolidation pressure
       # p0_dot = Compute_p0_dot()   
        self.M = 1.2#6*np.sin(phi)/(3 + np.sin(phi))
        #Material parameters
        v = 0.3        
        E = 20000e3*(2*(1+v))
        self.K = 36626762.87673086#E/(3*(1-2*v))
        self.G = 20000e3#E/(2*(1+v))
        

    def Draw_yield_surface(self,pc):
        q = 0; p_eff=0
        p_yield = []
        q_MCC= []
        while p_eff<pc:        
          #solve for q at f=0'
          q = np.sqrt(-self.M*self.M*p_eff*(p_eff-pc)) #MCC
          p_yield.append(p_eff)
          q_MCC.append(q)
          p_eff = p_eff + 10
                            
        return p_yield, q_MCC
    
    def Calculate_Elastic_trial_state(self):
        du_norm = 0
        dt = 0.01    
        eps1_dot = -0.01
        tol = 1e-6
        n=2000
        sig3 = -2e5
        N = 1.788
        lambda_e =0.077
        kappa = 0.0066
        
        # calculate residuals of non-linear equations in trial state

        
        #storage vectors
        eps1_vec = np.zeros(n); eps_p_p = np.zeros(n); eps_q_p = np.zeros(n)
        p_vec = np.zeros(n); q_vec = np.zeros(n); sig1_vec = np.zeros(n); 
        pc_vec = np.zeros(n); v_vec = np.zeros(n);     p_yield = []; q_yield=[]  
        K = np.zeros(n);
        
        #initialize sig1
        sig1_vec[0] = sig3
        d_evol_p = 0     
        d_es_p =0        
        pc_vec[0] =self.pc      
        
        # running drained Tx compression text
        
        for i in range(n-1):
            K[i] = self.K# - 2*eps1_vec[i] * eps1_vec[i])
            
        # TODO
        # ADD pressure dependent bulk modulus
        # decrease in pressure after first yield, check that!
 
            
        #rate form (elastic)

            # bulk moduou           
            sig1_dot = pow((1+3*K[i]/(self.G)),-1)*(9*K[i]*eps1_dot)    
            sig1_vec[i+1] =  sig1_vec[i] + dt*sig1_dot 

            eps3_dot = -sig1_dot/(2*self.G) + eps1_dot

            #trial invariants p-q            
            p_vec[i+1] = -1/3*(sig1_vec[i+1] + 2*sig3) -K[i]*eps_p_p[i]
            q_vec[i+1] = -(sig1_vec[i+1] -sig3) - 3*self.G*eps_q_p[i]
            
            v_vec[i+1] = N - lambda_e*math.log(p_vec[i+1])   
            
            print(v_vec[i+1]*p_vec[i+1]/0.0066)
            
            # initial yield surface
            pc_vec[i+1] = self.pc  
            
            # yield surface            
            f = pow(q_vec[i+1],2) + self.M*self.M*p_vec[i+1]*(p_vec[i+1]-pc_vec[i+1]) #MCC


            
            if f>0: # yielding
                  
                  df_dp =2*self.M*self.M*p_vec[i+1] - self.M*self.M*pc_vec[i+1]
                  df_dq = 2*q_vec[i+1]
                  df_dpc = -self.M*self.M*p_vec[i+1]
                  a = -v_vec[i+1] /(lambda_e-kappa) # negative sign to comply with signs of eps invariants: evol and es
                  
                  es_dot = -2/3*(eps1_dot - eps3_dot)
                  evol_dot = -(eps1_dot + 2*eps3_dot)
                 
                  lambda_dot =(df_dq*3*self.G*es_dot + df_dp*K[i]*evol_dot)/( df_dp*K[i]*df_dp+ 3*df_dq*self.G*df_dq  + df_dpc*a*df_dp*pc_vec[i])       
      
                  d_evol_p = dt*lambda_dot*df_dp
                  d_es_p = dt*lambda_dot*df_dq
                  
                  # update
                  eps_p_p[i+1] = eps_p_p[i] + d_evol_p
                  eps_q_p[i+1] =  eps_q_p[i]  + d_es_p



            
            # corrected p-q values
            p_vec[i+1] = p_vec[i+1] -   K[i]*d_evol_p 
            q_vec[i+1] = q_vec[i+1] -   3*self.G*d_es_p 
                           
            
            # expansion of yield surface (hardening)
            d_pc = (d_evol_p) *(v_vec[i+1]/(lambda_e-kappa))*pc_vec[i]   #v_vec[i+1]
            pc_vec[i+1] =pc_vec[i] + d_pc
            eps1_vec[i+1] = eps1_vec[i] + dt*eps1_dot
            
            
            self.Plots(eps1_vec,eps_p_p,pc_vec,q_vec,sig1_vec,p_vec,n,i)               
            
                  
            #NCL
                  

    def Plots(self,eps1_vec,eps_p_p,pc_vec,q_vec,sig1_vec,p_vec,n,i):

        if i == n-2:
            plt.show()                            
            #plt.figure(1)
            plt.subplot(411)        
            plt.plot(eps1_vec, q_vec)
            plt.subplot(412)        
            plt.plot(eps1_vec, p_vec)            
            plt.subplot(413)                
            plt.plot(p_vec,q_vec, '--',label='Stress path',c='b',linewidth=2.5)
            plt.plot(p_vec,1.2*p_vec, '--',label='CSL',c='r',linewidth=2.5)
            
            plt.subplot(414)                        
            plt.plot(eps1_vec, eps_p_p)
        if i % 400 == 0:
            plt.subplot(413)                        
            #draw yield surface                  
            p_yield, q_yield = self.Draw_yield_surface(pc_vec[i])        
            plt.plot(p_yield,q_yield)                
            plt.show()        

    

    

  
    
    
test = MCC()
test.Calculate_Elastic_trial_state()
#print(p,q_)

#plt.plot(p,q_MCC)
#plt.plot(p_sp1,q_sp1)
