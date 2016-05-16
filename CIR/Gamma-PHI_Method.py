# -*- coding: utf-8 -*-
"""
Created on Wed May 20 18:53:25 2015

@author: Charles
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot




#DMF
T_c = 562.2  #K
P_c = 48.9*100 #kpa
w_c = 0.212
a   = 0#1.92E-2
rho_d = 876
MM_d  = 78.11
V_liq1 = (MM_d/rho_d)#77.44/1000 #m3/kmol or L/mol
b = 0

# water
T_c1 = 553.5 #k
P_c1 = 40.73*100 #bar
w_c1 = 0.211

rho_w = 779
MM_w  = 87.16
V_liq2 = MM_d/rho_w #m3/kmol
a1 = 0#-1.09e-2
b1 = 0
R=8.31424
T = 25 + 273 # K
# Antone constants

 
A_1 = 6.87947202
B_1 = 1196.76
C_1 = 219.161
A_2 = 6.85146
B_2 = 1206.47
C_2 = 223.136

def Psats(T):
    conversion_factor = 101.325/760 #from mmHg to kPa
    P_satA =(10**((A_1)-(B_1/(C_1+(T-273.15)))))*conversion_factor
    P_satB =(10**((A_2)-(B_2/(C_2+(T-273.15)))))*conversion_factor
    return [P_satA, P_satB]
# Wilson Coefficient   
A21 = 0.6
A12 = 0.96
# B models
def B_model(T):
    TrA =T/T_c
    TrB =T/T_c1
    #Component A
    B_o =(0.1445)-(0.333/TrA )-(0.1385/(TrA **2))-(0.0121/(TrA **3))-(0.000607/(TrA **8))
    B_2 =0.0637+(0.331/TrA**2)-(0.423/TrA **3)-(0.008/TrA **8)
    B_3 =(a/TrA**6)-(b/TrA**8)    
    B_comp_A =((R*T_c)/P_c)*(B_o+(w_c*B_2)+B_3)
    
    #Component B
    B_o =(0.1445)-(0.333/TrB )-(0.1385/(TrB **2))-(0.0121/(TrB **3))-(0.000607/(TrB **8))
    B_2 =0.0637+(0.331/TrB**2)-(0.423/TrB**3)-(0.008/TrB**8)
    B_3 =(a1/TrB **6)-(b1/TrB **8)
    B_comp_B =((R*T_c1)/P_c1)*(B_o+(w_c1*B_2)+B_3)
    
    return [B_comp_A, B_comp_B]

    
def fugacity(T, P):
    B_T = B_model(T)
    fug_A = np.exp(B_T[0]*P/(R*T))
    fug_B = np.exp(B_T[1]*P/(R*T))
    return [fug_A, fug_B]

#fugacity at saturated conditions
def fuga_sat(T):
    P_sat = Psats(T)
    B_T = B_model(T)
    fug_sat_A = np.exp(B_T[0]*P_sat[0]/(R*T))
    fug_sat_B = np.exp(B_T[1]*P_sat[1]/(R*T))
    return [fug_sat_A, fug_sat_B]
    
def phi(T, P):
    P_sat = Psats(T)
    fug = fugacity(T, P)
    fug_sat = fuga_sat(T)
    phi_A = (fug[0]/fug_sat[0])*np.exp(-V_liq1*(P-P_sat[0])/(R*T))
    phi_B = (fug[1]/fug_sat[1])*np.exp(-V_liq2*(P-P_sat[1])/(R*T))
    return [phi_A, phi_B]

def lngamma(x):
    x_2 = 1 - x
    lngamma_A = -np.log(x + (x_2*A12)) + x_2*((A12/(x + (x_2*A12))) - (A21/(x_2 + (x*A21)))) 
    lngamma_B = -np.log(x_2 + (x*A21)) - x*((A12/(x + (x_2*A12))) - (A21/(x_2 + (x*A21))))
    return [lngamma_A, lngamma_B]
    
def gamma(x):
    ln_gamma = lngamma(x)
    gamma_A = np.exp(ln_gamma[0])
    gamma_B = np.exp(ln_gamma[1])
    return [gamma_A, gamma_B]
    
    
def solver(T, y):
    error = 1e-5
    xo = 0.1
    y2 = 1 - y
    Po = 11.6
    difference = 1
    while difference > error:
        gam = gamma(xo)
        phi_list = phi(T, Po)
        psat = Psats(T)
        
        F_p1 = gam[0]
        F_p2 = gam[1]
        F_p3 = y*phi_list[0]/psat[0]
        F_p4 = y2*phi_list[1]/psat[1]
        
        mat1 = np.array([[-F_p1, F_p3],
                         [ F_p2, F_p4]])
        mat2 = np.array([[0],
                         [F_p2]])
         
        mat1 = np.linalg.inv(mat1)                
        x = np.dot(mat1,mat2)
        difference = abs(x[0,0] - xo)
        P = x[1,0]
        x = x[0,0]        
        xo = x
        Po = P
        
    return [x, P]   
    

npoints = 1000
yseries = np.linspace(0.01, 1, npoints)
xi = []
yi = []


xe = [
0.1035,
0.175,
0.2760,
0.377,
0.433,
0.509,
0.583,
0.694,
0.7945,
0.9005,
0.9500
]

ye = [
0.1375,
0.217,
0.313,
0.4015,
0.446,
0.505,
0.562,
0.6505,
0.741,
0.8565,
0.922]


pe = [102.05,
104.5,
106.75,
108.1,
108.45,
108.65,
108.3,
106.9,
104.5,
100.6,
98.15]

pe = (101/760)*np.array(pe)

print pe
for i in range(npoints):
    
    y = yseries[i]
    plotter = solver(T, y)
    xi.append(plotter[0])
    yi.append(plotter[1])
plot.subplot(2,1,1)
plot.plot(yseries, yi,xi, yi, ye, pe, 'ro', xe, pe, 'bo')

plot.subplot(2,1,2)
plot.plot(xe, ye, 'ro', xi, yseries, yseries, yseries)
plot.show()    
