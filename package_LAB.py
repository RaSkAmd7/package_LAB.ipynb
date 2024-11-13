import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

#-----------------------------------        
def Lead_Lag_RT(MV,Kp,T_lead,T_lag,Ts,PV,PVInit=0,method='EBD'):
    
    """
Lead_Lag_RT(MV,Kp,T_lead,T_lag,Ts,PV,PVInit=0,method='EBD')
    The function "Lead_Lag_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :T_lead: lead time constant [s]
    :T_lag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
    if (T_lag != 0):
        K = Ts/T_lag
        alpha = T_lead/Ts
        
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+alpha)*MV[-1])-alpha*MV[-2]))
            elif method == 'EFD':
                PV.append((1-K)*PV[-1] + K*Kp*(alpha*MV[-1]+(1-alpha)*MV[-2]))
            #elif method == 'TRAP':
               #PV.append((1/(2*T+Ts))*((2*T-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))            
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*MV[-1])
    else:
        PV.append(Kp*MV[-1])
#----------------------------------------------------------------------------------------------------------------------
def calculate_pid(E,T_FD,T_s,T_I,K_C,T_D,PV,SP,MV_I_Init=0,MV_D_Init=0):
    global MV_I, MV_D, MV_P
    
    E.append(PV[-1] - SP[-1]) 
    if len(E) < 2:
        E.append(PV[-1] - SP[-1])
       
    # Proportional term
    MV_P = K_C * E[-1]

    # Integral term (MV_I[k] = MV_I[k-1] + K_C * Ts / T_I * E[k])
    if len(MV_I) == 0:
        MV_I.append(MV_I_Init)
    MV_I.append(MV_I[-1] + (K_C * T_s / T_I) * E[-1])

    # Derivative term (MV_D[k] = ...)
    if len(MV_D) == 0:
        MV_D.append(MV_D_Init)
    MV_D.append((T_FD / (T_FD + T_s)) * MV_D[-1] + (K_C * T_D / (T_FD + T_s)) * (E[-1] - E[-2]))


    # Total output (MV[k] = MV_P[k] + MV_I[k] + MV_D[k])
    MV.append(MV_P[-1] + MV_I[-1] + MV_D[-1])
