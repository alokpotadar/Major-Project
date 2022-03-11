import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable, Table, Column
from astropy import units as u
import csv
import argparse
#from np import np.linalg.solve


#initializing constants
#----------------------

#MIDSPAN DATA
parser=argparse.ArgumentParser()
parser.add_argument("wavelength", type=str)
parser.add_argument("amplitude",type=str)
args=parser.parse_args()
print(args.wavelength)
s_upon_w = int(args.wavelength)

midspan_chord=0.07
midspan_Cl=6.80
midspan_incidence_arr=np.linspace(1,12,12)

#WINGTIP DATA
tip_chord=1.524
tip_Cl=5.8
tip_incidence=3.5

total_span=0.495

wing_area=(midspan_chord)*total_span
print("wing area",wing_area)
aspect_ratio=(np.power(total_span,2))/wing_area

V=25
n=[2500]
Cl_sine_arr2=[]
Cd_sine_arr2=[]
T_ARR=[]
TT0_ARR=[]
ZS_ARR=[]
THETA_ARR=[]
C_l_spanwise=[]

lamda=total_span*s_upon_w/(16*32)
print(lamda)
amplitude=int(args.amplitude)/100
for segment in n:
    
    span_divs=np.linspace(0,total_span,segment+1)[1:]
    y_s=span_divs - (span_divs//lamda)*lamda
    #print('y_s:::::::::::::',y_s)
    for midspan_incidence in midspan_incidence_arr:
        print("----------THIS IS FOR "+str(midspan_incidence)+" AoA ----------------");c_s=[]
        theta_arr=np.linspace(0,(np.pi)/2,segment+1)[1:]
        #print('length of theta array',theta_arr)
        print(midspan_Cl)
        d_arr=[]
        l_full_arr=[]
        for i,theta in enumerate(theta_arr):
            l_arr=[]
            k=0
            c=midspan_chord + amplitude*0.5*(np.sin(y_s[i]*2*np.pi/lamda));c_s.append(c);''';theta = np.arctan((c*np.tan(theta))/midspan_chord)'''
            a=midspan_Cl*(1 + ((midspan_Cl - tip_Cl)*k/ midspan_Cl))
            alpha_degree=midspan_incidence*(1+k*((midspan_incidence-tip_incidence)/midspan_incidence))
            alpha_radian=np.radians(alpha_degree)
            mu=(c*a*0.25)/total_span
            muxalpha=(mu*alpha_radian)
            #Multiplying by sin(theta)
            co_effs=[0 for i in range(segment)]
            for i in range(segment):
                # print(i)
                co_effs[i] = (np.sin(((2*i)+1)*theta))*((np.sin(theta)) + (((2*i)+1)*mu))
                l_arr.append(co_effs[i])
    
            # c1=np.sin(theta)*((np.sin(theta)) + mu)
            # c3=np.sin(3*theta)*(np.sin(theta) + (3*mu))
            # c5=np.sin(5*theta)*(np.sin(theta) + (5*mu))
            # c7=np.sin(7*theta)*(np.sin(theta) + (7*mu))
            muxalphaxsin=muxalpha*np.sin(theta)
            d_arr.append(muxalphaxsin)
            # l_arr.append(np.round(c1,7))
            # l_arr.append(np.round(c3,7))
            # l_arr.append(np.round(c5,7))
            # l_arr.append(np.round(c7,7))
            l_full_arr.append(l_arr)
    
        # print(d_arr)
        # print("D",np.shape(d_arr))
        # print(l_full_arr)
        # print("L",np.shape(l_full_arr))
        coeff_arr=np.array(l_full_arr)
        d_arr=np.array(d_arr)
        x=np.linalg.solve(coeff_arr,d_arr)
        # print(x)
        print(np.allclose(np.dot(coeff_arr,x),d_arr));plt.plot(theta_arr,c_s);plt.close() #This is to check if our solution is correct. If True, the solution is correct
        theta_arr=np.linspace(0,(np.pi)/2,segment+1)
        z_upon_s_arr=[]
        T_arr=[]
        T_upon_T0=[]
        for t in theta_arr:
            z=0.5*np.cos(t)
            T=0.0
            for i in range(segment):
                T=T+x[i]*np.sin(((2*i)+1)*t)
            T=4*total_span*0.5*V*T
            # print(z)
            # print(T)
            T_arr.append(T)
            z_upon_s_arr.append(z)
        T_upon_T0=np.round(T_arr/T_arr[-1],5)
        # print(T_arr)
        # print(z_upon_s_arr)
        # print(T_upon_T0)
    
        ###Aerodynamic characteristics
        delta=0.0
        for i in range(segment):
            delta=delta+(((2*i)+1)*((x[i])**2))
    
        delta=delta/((x[0])**2)
        delta=delta-1
        # print(np.round(delta,6))
    
        C_l=np.pi*aspect_ratio*x[0]
        C_d=(C_l**2)*(1+delta)/(np.pi*aspect_ratio)
    
        print("Cl",C_l)
        print(C_d)
    
        rho=1.257
        lift=C_l*0.5*rho*(V**2)*wing_area
        induced_drag=C_d*0.5*rho*(V**2)*wing_area
        Cl_sine_arr2.append(str(C_l))
        Cd_sine_arr2.append(str(C_d))
        C_l_spanwise.append(np.pi*aspect_ratio*np.array(x))
        print(lift)
        print(induced_drag)
        theta_arr=np.array(theta_arr)*180/np.pi
        T_ARR.append(str(T_arr))
        TT0_ARR.append(str(T_upon_T0))
        ZS_ARR.append(z_upon_s_arr)
        THETA_ARR.append(theta_arr)

################SAVING ALL THE NUMPY ARRAYS ###################


with open('Sine_A_'+str(amplitude/10)+'_w_'+args.wavelength+'new.csv','w',newline='\n') as file:
    writer = csv.writer(file,delimiter='|')
    writer.writerow(Cl_sine_arr2)
    writer.writerow(Cd_sine_arr2)
    writer.writerow(T_ARR)
    writer.writerow(TT0_ARR)
    writer.writerow(C_l_spanwise)
