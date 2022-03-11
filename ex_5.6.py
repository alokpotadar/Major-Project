import numpy as np
import matplotlib.pyplot as plt
#from np import np.linalg.solve


#initializing constants
#----------------------

#MIDSPAN DATA
midspan_chord=3.048
midspan_Cl=5.5
midspan_incidence=5.5

#WINGTIP DATA
tip_chord=1.524
tip_Cl=5.8
tip_incidence=3.5

total_span=12.192
V=89.4
#----------------------

#Calculation of basic parameters

wing_area=(midspan_chord+tip_chord)*total_span/2
print("wing area",wing_area)
aspect_ratio=(np.power(total_span,2))/wing_area


'''Generating equations of the airfoil characteristics
   ---------------------------------------------------

We know that:
z/s=-(cos(theta)/2)
let:
z/s=-(cos(theta)/2)=k

'''
n=[4,10,25,50,100,250,500,1000,2000]
Cl_arr=[]
Cd_arr=[]
T_ARR=[]
TT0_ARR=[]
ZS_ARR=[]
THETA_ARR=[]
for segment in n:
	print("----------THIS IS FOR "+str(segment)+" SEGMENTS ----------------")
	theta_arr=np.linspace(0,(np.pi)/2,segment+1)[1:]
	# print(theta_arr)
	d_arr=[]
	l_full_arr=[]
	for theta in theta_arr:
		l_arr=[]
		k=-1*(np.cos(theta))
		c=midspan_chord*(1 + ((midspan_chord - tip_chord)*k/midspan_chord))
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
	print(np.allclose(np.dot(coeff_arr,x),d_arr)) #This is to check if our solution is correct. If True, the solution is correct

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
	Cl_arr.append(C_l)
	Cd_arr.append(C_d)
	print(lift)
	print(induced_drag)
	theta_arr=np.array(theta_arr)*180/np.pi
	T_ARR.append(T_arr)
	TT0_ARR.append(T_upon_T0)
	ZS_ARR.append(z_upon_s_arr)
	THETA_ARR.append(theta_arr)

for i in range(len(n)):
	plt.plot(THETA_ARR[i],ZS_ARR[i],linestyle='dashed',label=str(n[i])+" segments")
plt.ylabel('z/s')
plt.xlabel('theta')
plt.title('z/s vs theta')
plt.legend()
plt.savefig('zs.png')
plt.close()
for i in range(len(n)):
	plt.plot(THETA_ARR[i],T_ARR[i],linestyle='dashed',label=str(n[i])+" segments")
plt.ylabel('T')
plt.xlabel('theta')
plt.title('T vs theta')
plt.legend()
plt.savefig('T.png')
plt.close()
for i in range(len(n)):
	plt.plot(THETA_ARR[i],TT0_ARR[i],linestyle='dashed',label=str(n[i])+" segments")
plt.ylabel('TT0')
plt.xlabel('theta')
plt.title('T/T0 vs theta')
plt.legend()
plt.savefig('TT0.png')
plt.close()

plt.plot(n,Cl_arr,linestyle='dashed',marker="x",label='Cl')
plt.ylabel("Cl")
plt.xlabel("Number of Segments")
plt.legend()
plt.savefig('Cl.png')
plt.close()
plt.plot(n,Cd_arr,linestyle='dashed',marker="x",label='Cd')
plt.ylabel("Cd")
plt.xlabel("Number of Segments")
plt.legend()
plt.savefig('Cd.png')
plt.close()

Cl_percent_change=(100*(Cl_arr-Cl_arr[-1]))/Cl_arr[-1]
Cd_percent_change=(100*(Cd_arr-Cd_arr[-1]))/Cd_arr[-1]
c2=np.array(Cl_arr)/np.array(Cd_arr)
c2_percent_change=(100*(c2-c2[-1]))/c2[-1]
plt.plot(n,Cl_percent_change,linestyle='dashed',marker='.',label='Cl percent change')
plt.plot(n,Cd_percent_change,linestyle='dashed',marker='x',label='Cd percent change')
plt.plot(n,c2_percent_change,linestyle='dashed',marker='^',label='Cl/Cd percent change')
plt.ylabel('Percent Change')
plt.xlabel('Number of Segments')
plt.legend()
plt.savefig('Cl_Cd_percent_change.png')
plt.show()