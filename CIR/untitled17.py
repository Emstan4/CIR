# -*- coding: utf-8 -*-
"""
Created on Fri May 13 15:13:57 2016

@author: Top4
"""


from __future__ import division
import numpy as np
from matplotlib import pyplot as plot



#DMF
T_c = 562.2  #k
P_c = 48.9*100 #bar
w_c = 0.212
a   = 1.92E-2
V_liq2 = 89.37/1000 #m3/mol
b = 0
rho_d = 876
MM_d  = 78.11
V_d = (MM_d/rho_d)
 

# water
T_c1 =  553.5#k
P_c1 = 40.73*100 #bar
w_c1 = 0.211
V_liq1 = 108.04/1000 #m3/mol
a1 = 0#-1.09e-2
b1 = 0
R=8.31424
T = 25 + 273 # kPa
rho_w = 779
MM_w  = 80.74
V_w = (MM_w/rho_w)

# Antone constants

A_2 = 6.87947202
B_2 = 1196.76
C_2 = 219.161
A_1 = 6.85146
B_1 = 1206.47
C_1 = 223.136


P_sat1 =(10**((A_1)-((B_1)/(C_1+(T-273.15)))))*0.133333
P_sat2 =(10**((A_2)-(B_2/(C_2+(T-273.15)))))*0.133333
    
# Wilson Coefficient

A21 = 0.8731
A12 = 1.295


A=106.685953468227
B=106.595763947678
C=810.49994847524
A11=0.544750013862755*-100
B11=81.05
C11=892
B1= A - B*np.exp(C/T)
B2= A11 - B11*np.exp(C11/T)
def gamma_phi(P):
    xo = 0
    error = 1e-3
    diff = 2
    yo = 0
    
    while diff > error:
        x_2 = 1 - xo
        #Water

        #T_r1 =T/T_c1
#    #Water
#    
#        B_o =(0.1445)-(0.333/T_r1 )-(0.1385/(T_r1 **2))-(0.0121/(T_r1 **3))-(0.000607/(T_r1 **8))
#        B_2 =0.0637+(0.331/T_r1**2)-(0.423/T_r1**3)-(0.008/T_r1**8)
#        B_3 =(a1/T_r1 **6)-(b1/T_r1 **8)
#        B1 =((R*T_c1)/P_c1)*(B_o+(w_c1*B_2)+B_3)
       
        phi1 = np.exp(B1*P/(R*(T)))
        phi_sat1= np.exp(B1*P_sat1/(R*T))
        
        #DMF
#        T_r =T/T_c
#        B_o2 =(0.1445)-(0.333/T_r )-(0.1385/(T_r **2))-(0.0121/(T_r **3))-(0.000607/(T_r **8))
#        B_22 =0.0637+(0.331/T_r **2)-(0.423/T_r **3)-(0.008/T_r **8)
#        B_32 =(a/T_r **6)-(b/T_r **8)
#        B2 =((R*T_c)/P_c)*(B_o2+(w_c*B_22)+B_32)
        phi2 = np.exp(B2*P/(R*(T)))
        phi_sat2= np.exp(B2*P_sat2/(R*T))
        
        
        
        phi_1a=(phi1/phi_sat1)*np.exp(-V_liq1*(P-P_sat1)/(R*T))
        
        
        phi_2a=(phi2/phi_sat2)*np.exp(-V_liq2*(P-P_sat2)/(R*T))
        
        
        lngamma1=-np.log(xo + (x_2*A12)) + x_2*((A12/(xo + (x_2*A12))) - (A21/(x_2 + (xo*A21))))
        gamma1=np.exp(lngamma1)
        
        
        lngamma2=-np.log(x_2 + (xo*A21)) - xo*((A12/(xo + (x_2*A12)))-(A21/(x_2+(xo*A21))))
        gamma2=np.exp(lngamma2)
        
        F_p1 = gamma1
        F_p2 = gamma2
        F_pp1 = (P*phi_1a/P_sat1) 
        F_pp2 = (P*phi_2a/P_sat2) 
        diff = (abs(yo*F_pp1 - xo*F_p1) + abs((1 - yo)*F_pp2 - ((1 - xo)*F_p2)))
        x = (F_pp1*(F_pp2 - F_p2)/(F_p1*F_pp2 - F_pp1*F_p2))
        y = x*(F_p1/F_pp1)
        xo = x
        yo = y
    
    
    
    return [y, x] 
    
# Experimental Data

xe = [0, 0.0977, 0.1301, 0.3010, 0.5036, 0.6201, 0.6315, 0.7267, 0.80, 0.857, 0.8956, 0.9412, 0.9447, 0.9725]
ye = [0, 0.4133, 0.4624, 0.699, 0.8344, 0.8913, 0.8960, 0.9292, 0.9552, 0.9686, 0.9784, 0.9872, 0.9880, 0.9940]

s = 273
Te = [109.5+s, 98.15+s, 96.6+s, 86.7+s, 79.05+s, 75.85+s, 75.2+s, 72.35+s, 70.7+s, 69.35+s, 68.6+s, 67.6+s, 67.05+s, 67.35+s]



pplot = []
multi = 101.325/760
npoints = 1000
pseries = np.linspace(95.05*multi, 97*multi, npoints)
xi = []
yi = []
pe =  np.array([9.267461727,10.97176341,13.65783898,16.37613847,16.9867583,19.43398832,21.42718455,24.03466462,26.74904395,26.844185,27.2972197,29.00787015,31.04486903,34.14598065,35.62870018,35.99599519,37.41619857,39.15074369,42.03395549,43.44955037,44.87044898,46.98272585,50.28499502,53.2317943,55.19282956])
pe_mod = multi*pe
xe = [0,0.0441,0.1141,0.1853,0.2013,0.2652,0.3166,0.3821,0.4470,0.4492,0.4596,0.4977,0.5406,0.6011,0.6283,0.6349,0.6600,0.6900,0.7394,0.7640,0.7893,0.8286,0.8952,0.9583,1.0000]
ye = [0, 0.1989,0.3891,0.5439,0.5710,0.6526,0.7054,0.7606,0.8099,0.8113,0.8179,0.8415,0.8670,0.8963,0.9088,0.9116,0.9225,0.9343,0.9500,0.9572,0.9635,0.9731,0.9854,0.9949,1.0000]
for i in range(npoints):
    
    p = pseries[i]
    plotter = gamma_phi(p)
    xi.append((plotter[0]))
    yi.append(plotter[1])
plot.xlabel("$x$")  
plot.ylabel("$P$ $kPa$") 
plot.title("$Pxy$")   
plot.plot(yi,xi)#, xe, pe_mod, 'ro')
plot.show()