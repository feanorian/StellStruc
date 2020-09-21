
import scipy.constants as sp
import numpy as np
import matplotlib.pyplot as plt
import sys 

####Here will be the constants we will use throughout the simulation
#Mass of sun
M_Sun = 1.989e30
# Radius of sun
R_Sun=6.96e8  
#luminosity os sun
L_Sun=3.828e26   

pi = np.pi

#Boltzmann constant
k_b = 1.38e-23
#speed of light
c = sp.c

# Radiation constant
a=7.5657e-16     
# Gravitational constant
G = sp.G
#polytropic index
gamma = 5./3.
# Mass of hydrogen atom in Kg
m_H=1.67e-27   

#Model Parameters with elemental abundances and mean molecular weight
#M = raw_input('Enter the Mass of the star (in solar masses):  ')
M_Star =  M_Sun 
if M_Star < 1.66*M_Sun:
	R_Star=R_Sun*1.06*(M_Star/M_Sun)**0.945
else:
	R_Star=R_Sun*1.33*(M_Star/M_Sun)**0.555
X=0.70
Y=0.28
Z=0.02
mu=1./(2.*X+0.75*Y+0.5*Z)
# energy generation rates for PP, CNO and Kramer's opacity
epsilon_pp=2.6E-37*X**2
#epsilon_CNO=7.9E-118*X*Z
epsilon_CNO=0
kappa0=4.3E21*Z*(X+1)

#tolerance
tol = 1.E-3


#Here we establish the number of points to evaluate and numIterations the simulation we
dim=10000 # Number of grid points
numIter=200 # maximum number of numIterations

#measure the difference between
k_values=np.zeros(numIter)
R_diff=np.zeros(numIter)
P_diff=np.zeros(numIter)
T_diff=np.zeros(numIter)
L_diff=np.zeros(numIter)

#these arrays store values of each variable integrating outward from center
M=np.zeros(dim)
R=np.zeros(dim)
T=np.zeros(dim)
P=np.zeros(dim)
L=np.zeros(dim)
rho=np.zeros(dim)
kappa=np.zeros(dim)
epsilon=np.zeros(dim)

##these arrays will store the M, R, T,eT_c. variables integrating from the srface inward
MS=np.zeros(dim)
RS=np.zeros(dim)
TS=np.zeros(dim)
PS=np.zeros(dim)
LS=np.zeros(dim)
rhoS=np.zeros(dim)
kappaS=np.zeros(dim)
epsilonS=np.zeros(dim)

##these arrays store the M, R, T,eT_c. variables after perturbation of conditions at center integrating from the center
Mp=np.zeros(dim)
Rp=np.zeros(dim)
Tp=np.zeros(dim)
Pp=np.zeros(dim)
Lp=np.zeros(dim)
rhop=np.zeros(dim)
kappap=np.zeros(dim)
epsilonp=np.zeros(dim)


##these arrays store the M, R, T,eT_c. variables after perturbation of conditions at surface integrating from the surface
MSp=np.zeros(dim)
RSp=np.zeros(dim)
TSp=np.zeros(dim)
PSp=np.zeros(dim)
LSp=np.zeros(dim)
rhoSp=np.zeros(dim)
kappaSp=np.zeros(dim)
epsilonSp=np.zeros(dim)



#initial boundary conditions

M[0]=1.E-06*M_Star
# Mass of each mass shell
dM=(M_Star-M[0])/(dim-1) 
#Average Core Pressure
P_c=1.E15            
# Surface pressure
P_s=1e8

# Average Core Temperature
T_c=1.E7                
# Surface temperature
T_s=5600.                
				 
# Luminosity
L_s=L_Sun             
# Stellar radius
R_s=1.*R_Star     

for k in range (0,numIter):
	P[0]=P_c
	T[0]=T_c
	rho[0]=P[0]*mu*m_H/(k_b*T[0])
	R[0]=(3.*M[0]/(4.*pi*rho[0]))**(1./3.)
	epsilon[0]=epsilon_pp*rho[0]*T[0]**4.+ epsilon_CNO*rho[0]*T[0]**16
	L[0]=epsilon[0]*M[0]
	kappa[0]=kappa0*rho[0]*T[0]**(-3.5)
	
	#integration from center
	for i in range (1,dim/2):
		R[i] = R[i-1] + dM/(4.*pi*R[i-1]**2*rho[i-1])
		P[i] = P[i-1] - dM*G*M[i-1]/(4.*pi*R[i-1]**4)   
		L[i] = L[i-1] + dM*epsilon[i-1]
		nabla_rad = (3.*kappa[i-1]*L[i-1]*P[i-1]/(16.*pi*a*c*T[i-1]**4*G*M[i-1]))
		#this condition determines if T is in radiative or convection zone, and calculates proper temp
		if nabla_rad < (gamma-1)/gamma:
			T[i] = T[i-1] - (dM*3.*kappa[i-1]*L[i-1]/
				(16.*pi*a*c*R[i-1]**2*T[i-1]**3)/(4.*pi*R[i-1]**2))			
		else:
			T[i] = T[i-1] + (dM*(gamma-1.)/gamma*
				T[i-1]/P[i-1]*(P[i]-P[i-1])/dM)					
	
	#condition to ensure T and P are positive
		if T[i] <= 0. or P[i] <= 0.:
			print 'T[i] <= 0. or P[i] <= 0.',i,T[i],P[i]
			break
	
		M[i] = M[i-1] + dM
		rho[i] = P[i]*mu*m_H/(k_b*T[i])
		epsilon[i] = (epsilon_pp*rho[i]*T[i]**4.+ epsilon_CNO*rho[i]*T[i]**16)
		kappa[i] = kappa0*rho[i]*T[i]**(-3.5)
	
	#Boundary conditions on surface    
	MS[dim-1]=M_Star
	LS[dim-1]=L_s
	RS[dim-1]=R_s
	TS[dim-1]=T_s
	PS[dim-1]=P_s
	rhoS[dim-1]=PS[dim-1]*mu*m_H/(k_b*TS[dim-1])
	kappaS[dim-1]=kappa0*rhoS[dim-1]*TS[dim-1]**(-3.5)
	epsilonS[dim-1]=(epsilon_pp*rhoS[dim-1]*TS[dim-1]**4+
		epsilon_CNO*rhoS[dim-1]*TS[dim-1]**16)
	
	dr_R_Star=dM/(4.*pi*RS[dim-1]**3*rhoS[dim-1])
	if dr_R_Star/R_s > 1.E-02:
		print 'Ending run since radial step size is too large'
		print 'dr_R_Star/Rs=',dr_R_Star/R_s,dr_R_Star,R_s
		exit()
	
	#integration from surface
	for j in range (dim-2,dim/2-2,-1):
		RS[j] = RS[j+1] - dM/(4.*pi*RS[j+1]**2*rhoS[j+1])
		PS[j] = PS[j+1] + dM*G*MS[j+1]/(4.*pi*RS[j+1]**4)
		LS[j] = LS[j+1] - dM*epsilonS[j+1]
		nabla_rad = (3.*kappaS[j+1]*LS[j+1]*PS[j+1]/
			(16.*pi*a*c*TS[j+1]**4*G*MS[j+1]))
		if nabla_rad < (gamma-1)/gamma:
			TS[j] = TS[j+1] + (dM*3.*kappaS[j+1]*LS[j+1]/
				(16.*pi*a*c*RS[j+1]**2*TS[j+1]**3)/(4.*pi*RS[j+1]**2))
		else:
			TS[j] = TS[j+1] - (dM*(gamma-1.)/gamma * 
				TS[j+1]/PS[j+1]*(PS[j+1] - PS[j])/dM)
		
		
		
		if RS[j] <= 0. or LS[j] <= 0.:
			print 'RS[j] <= 0. or LS[j] <= 0.',j,RS[j],LS[j]
			break

		MS[j] = MS[j+1] - dM
		rhoS[j] = PS[j]*mu*m_H/(k_b*TS[j])
		epsilonS[j] = (epsilon_pp*rhoS[j]*TS[j]**4.+
			epsilon_CNO*rhoS[j]*TS[j]**16)
		kappaS[j] = kappa0*rhoS[j]*TS[j]**(-3.5)
	
	
	
	#Unperturbed Values at the midpoint. The integrations from the center and the surface will converge here
	Rmid0=R[dim/2-1]
	RSmid0=RS[dim/2-1]
	Pmid0=P[dim/2-1]
	PSmid0=PS[dim/2-1]
	Tmid0=T[dim/2-1]
	TSmid0=TS[dim/2-1]
	Lmid0=L[dim/2-1]
	LSmid0=LS[dim/2-1]
	k_values[k]=k
	R_diff[k] = Rmid0-RSmid0
	P_diff[k] = Pmid0-PSmid0
	T_diff[k] = Tmid0-TSmid0
	L_diff[k] = Lmid0-LSmid0	
	

	


	Rmidpoint=(R[dim/2-1]+RS[dim/2-1])/2.
	Pmidpoint=(P[dim/2-1]+PS[dim/2-1])/2.
	Tmidpoint=(T[dim/2-1]+TS[dim/2-1])/2.
	Lmidpoint=(L[dim/2-1]+LS[dim/2-1])/2.
	print '**************************************'
	print 'numIteration step k=',k
	print 'Convergence tolerance parameter=',tol
	print 'State of convergence shown below:'
	print 'R_diff[k])/Rmidpoint=',abs(R_diff[k])/Rmidpoint 
	print 'P_diff[k])/Pmidpoint=',abs(P_diff[k])/Pmidpoint 
	print 'T_diff[k])/Tmidpoint=',abs(T_diff[k])/Tmidpoint 
	print 'L_diff[k])/Lmidpoint=',abs(L_diff[k])/Lmidpoint 
	print '**************************************'
	print ' '





	#The halting condition if stellar structure equations can be solved given our initial conditions
	if abs(R_diff[k])/Rmidpoint < tol:
		if abs(P_diff[k])/Pmidpoint < tol:
			if abs(T_diff[k])/Tmidpoint < tol :
				if abs(L_diff[k])/Lmidpoint < tol:
					print 'Convergence achieved'
					break
# 1 percent perturbation of the pressure from the initial shooting solutions above
	Mp[0]=M[0]
	Pp[0]=1.01*P_c
	Tp[0]=T_c
	rhop[0]=Pp[0]*mu*m_H/(k_b*Tp[0])
	Rp[0]=(3.*Mp[0]/(4.*pi*rhop[0]))**(1./3.)
	epsilonp[0]=(epsilon_pp*rhop[0]*Tp[0]**4.+
		epsilon_CNO*rhop[0]*Tp[0]**16)
	Lp[0]=epsilonp[0]*Mp[0]
	kappap[0]=kappa0*rhop[0]*Tp[0]**(-3.5)
	#outward integration
	for i in range (1,dim/2):
		Rp[i] = Rp[i-1] + dM/(4.*pi*Rp[i-1]**2*rhop[i-1])
		Pp[i] = Pp[i-1] - dM*G*Mp[i-1]/(4.*pi*Rp[i-1]**4)
		Lp[i] = Lp[i-1] + dM*epsilonp[i-1]
		nabla_rad = (3.*kappap[i-1]*Lp[i-1]*Pp[i-1]/(16.*pi*a*c*Tp[i-1]**4*G*Mp[i-1]))
		if nabla_rad < (gamma-1)/gamma:
			Tp[i] = Tp[i-1] - (dM*3.*kappap[i-1]*Lp[i-1]/(16.*pi*a*c*Rp[i-1]**2*Tp[i-1]**3)/(4.*pi*Rp[i-1]**2))
		else:
			Tp[i] = Tp[i-1] + (dM*(gamma-1.)/gamma*Tp[i-1]/Pp[i-1]*(Pp[i]-Pp[i-1])/dM)
			
		if Tp[i] <= 0.  or Pp[i] <= 0.:
			print 'Tp[i] <= 0. or Pp[i] <= 0.',i,Tp[i],Pp[i]
			break
				
		Mp[i] = Mp[i-1] + dM
		rhop[i] = Pp[i]*mu*m_H/(k_b*Tp[i])
		epsilonp[i] = (epsilon_pp*rhop[i]*Tp[i]**4.+ epsilon_CNO*rhop[i]*Tp[i]**16)
		kappap[i] = kappa0*rhop[i]*Tp[i]**(-3.5)
			
			
	#Takes the midpoint from the perturbed solution and computes the differnce
	#of the midpointP from the initial solution
	Rmid1=Rp[dim/2-1]
	Pmid1=Pp[dim/2-1]
	Tmid1=Tp[dim/2-1]
	Lmid1=Lp[dim/2-1]
	Rdiff1 = Rmid1-RSmid0
	Pdiff1 = Pmid1-PSmid0
	Tdiff1 = Tmid1-TSmid0
	Ldiff1 = Lmid1-LSmid0
	
#Initial conditions	with 1 percent temp change
	Mp[0]=M[0] 
	Pp[0]=P_c
	Tp[0]=T_c*1.01
	rhop[0]=Pp[0]*mu*m_H/(k_b*Tp[0])
	Rp[0]=(3.*Mp[0]/(4.*pi*rhop[0]))**(1./3.)
	epsilonp[0]=(epsilon_pp*rhop[0]*Tp[0]**4.+ epsilon_CNO*rhop[0]*Tp[0]**16)
	Lp[0]=epsilonp[0]*Mp[0]
	kappap[0]=kappa0*rhop[0]*Tp[0]**(-3.5)      

	#integrate from the center with 1% Temp perturbation
	for i in range (1,dim/2):
		Rp[i] = Rp[i-1] + dM/(4.*pi*Rp[i-1]**2*rhop[i-1])
		Pp[i] = Pp[i-1] - dM*G*Mp[i-1]/(4.*pi*Rp[i-1]**4)
		Lp[i] = Lp[i-1] + dM*epsilonp[i-1]
		nabla_rad = (3.*kappap[i-1]*Lp[i-1]*Pp[i-1]/(16.*pi*a*c*Tp[i-1]**4*G*Mp[i-1]))

		#this condition determines if T is in radiative or convection zone, and calculates proper temp
		if nabla_rad < (gamma-1)/gamma:
			Tp[i] = Tp[i-1] - (dM*3.*kappap[i-1]*Lp[i-1]/(16.*pi*a*c*Rp[i-1]**2*Tp[i-1]**3)/(4.*pi*Rp[i-1]**2))
						
		else:
			Tp[i] = Tp[i-1] + (dM*(gamma-1.)/gamma*Tp[i-1]/Pp[i-1]*(Pp[i]-Pp[i-1])/dM)
			
		if Tp[i] <= 0. or Pp[i] <= 0.:
			print 'Tp[i] <= 0. or Pp[i] <= 0.',i,Tp[i],Pp[i]
			break

		Mp[i] = Mp[i-1] + dM
		rhop[i] = Pp[i]*mu*m_H/(k_b*Tp[i])
		epsilonp[i] = (epsilon_pp*rhop[i]*Tp[i]**4.+epsilon_CNO*rhop[i]*Tp[i]**16)                            
		kappap[i] = kappa0*rhop[i]*Tp[i]**(-3.5)                           



	Rmid2=Rp[dim/2-1]
	Pmid2=Pp[dim/2-1]
	Tmid2=Tp[dim/2-1]
	Lmid2=Lp[dim/2-1]
	Rdiff2 = Rmid2-RSmid0
	Pdiff2 = Pmid2-PSmid0
	Tdiff2 = Tmid2-TSmid0
	Ldiff2 = Lmid2-LSmid0


		
	MSp[dim-1]=M_Star
	LSp[dim-1]=L_s*1.01
	RSp[dim-1]=R_s
	TSp[dim-1]=T_s
	PSp[dim-1]=P_s
	rhoSp[dim-1]=PSp[dim-1]*mu*m_H/(k_b*TSp[dim-1])
	kappaSp[dim-1]=kappa0*rhoSp[dim-1]*TSp[dim-1]**(-3.5)
	epsilonSp[dim-1]=(epsilon_pp*rhoSp[dim-1]*TSp[dim-1]**4+epsilon_CNO*rhoSp[dim-1]*TSp[dim-1]**16)
		
		
	dr_R_Star=dM/(4.*pi*RS[dim-1]**3*rhoS[dim-1])
	if dr_R_Star/R_s > 1.E-02:
		print 'Ending run since radial step size is too large'
		print 'dr_R_Star/Rs=',dr_R_Star/R_s,dr_R_Star,R_s
		exit()
		
		
		
	for j in range (dim-2,dim/2-2,-1):
		RSp[j] = RSp[j+1] - dM/(4.*pi*RSp[j+1]**2*rhoSp[j+1])
		PSp[j] = PSp[j+1] + dM*G*MSp[j+1]/(4.*pi*RSp[j+1]**4)
		LSp[j] = LSp[j+1] - dM*epsilonSp[j+1]
		nabla_rad = (3.*kappaSp[j+1]*LSp[j+1]*PSp[j+1]/(16.*pi*a*c*TSp[j+1]**4*G*MSp[j+1]))
		if nabla_rad < (gamma-1)/gamma:
			TSp[j] = TSp[j+1] + (dM*3.*kappaSp[j+1]*LSp[j+1]/(16.*pi*a*c*RSp[j+1]**2*TSp[j+1]**3)/(4.*pi*RSp[j+1]**2))
		else:
			TSp[j] = TSp[j+1] - (dM*(gamma-1.)/gamma*TSp[j+1]/PSp[j+1]*(PSp[j+1] - PSp[j])/dM)	
				
		if RSp[j] <= 0. or LSp[j] <= 0.:
			print 'RSp[j] <= 0. or LSp[j] <= 0.',j,RSp[j],LSp[j]
			break

		MSp[j] = MSp[j+1] - dM
		rhoSp[j] = PSp[j]*mu*m_H/(k_b*TSp[j])
		epsilonSp[j] = (epsilon_pp*rhoSp[j]*TSp[j]**4.+epsilon_CNO*rhoSp[j]*TSp[j]**16)		
		kappaSp[j] = kappa0*rhoSp[j]*TSp[j]**(-3.5)		
		
		
	RSmid1=RSp[dim/2-1]
	PSmid1=PSp[dim/2-1]
	TSmid1=TSp[dim/2-1]
	LSmid1=LSp[dim/2-1]
	Rdiff3 = Rmid0-RSmid1
	Pdiff3 = Pmid0-PSmid1
	Tdiff3 = Tmid0-TSmid1
	Ldiff3 = Lmid0-LSmid1	
		
		
		

					
		
				
	#now we will perturb the radius 1% from initial solution        
	MSp[dim-1]=M_Star
	LSp[dim-1]=L_s
	RSp[dim-1]=R_s*1.01
	PSp[dim-1]=P_s
	TSp[dim-1]=T_s
	

	rhoSp[dim-1]=PSp[dim-1]*mu*m_H/(k_b*TSp[dim-1])
	kappaSp[dim-1]=kappa0*rhoSp[dim-1]*TSp[dim-1]**(-3.5)
	epsilonSp[dim-1]=(epsilon_pp*rhoSp[dim-1]*TSp[dim-1]**4+epsilon_CNO*rhoSp[dim-1]*TSp[dim-1]**16)
				   
	
	dr_R_Star=dM/(4.*pi*RS[dim-1]**3*rhoS[dim-1])
	if dr_R_Star/R_s > 1.E-02:
		print 'Ending run since radial step size is too large'
		print 'dr_R_Star/Rs=',dr_R_Star/R_s,dr_R_Star,R_s
	#integrates the perturbed radius from the surface
	for j in range(dim-2,dim/2-2,-1):
		RSp[j] = RSp[j+1] - dM/(4.*pi*RSp[j+1]**2*rhoSp[j+1])
		PSp[j] = PSp[j+1] + dM*G*MSp[j+1]/(4.*pi*RSp[j+1]**4)   
		


		LSp[j] = LSp[j+1] - dM*epsilonSp[j+1]
		nabla_rad = (3.*kappaSp[j+1]*LSp[j+1]*PSp[j+1]/(16.*pi*a*c*TSp[j+1]**4*G*MSp[j+1]))
								
		if nabla_rad < (gamma-1)/gamma:
			TSp[j] = TSp[j+1] + (dM*3.*kappaSp[j+1]*LSp[j+1]/(16.*pi*a*c*RSp[j+1]**2*TSp[j+1]**3)/(4.*pi*RSp[j+1]**2))
					  
		else:
			TSp[j] = TSp[j+1] - (dM*(gamma-1.)/gamma*TSp[j+1]/PSp[j+1]*(PSp[j+1] - PSp[j])/dM)

		#checking for overshoot of R and L
		if RSp[j] <= 0. or LSp[j] <= 0.:
			print 'RSp[j] <= 0. or LSp[j] <= 0.',j,RSp[j],LSp[j]
			break

		MSp[j] = MSp[j+1] - dM
		rhoSp[j] = PSp[j]*mu*m_H/(k_b*TSp[j])
		epsilonSp[j] = (epsilon_pp*rhoSp[j]*TSp[j]**4.+ epsilon_CNO*rhoSp[j]*TSp[j]**16)
		kappaSp[j] = kappa0*rhoSp[j]*TSp[j]**(-3.5)  
			
	#Determines the midpoint from the surface after perturbation
	RSmid2 = RSp[dim/2-1]
	PSmid2 = PSp[dim/2-1]
	TSmid2 = TSp[dim/2-1]
	LSmid2 = LSp[dim/2-1]
	Rdiff4 = Rmid0-RSmid2
	Pdiff4 = Pmid0-PSmid2
	Tdiff4 = Tmid0-TSmid2
	Ldiff4 = Lmid0-LSmid2
	
	
	#Here we construct a matrix whose elements are the derivatives of the changes at the midpoint relative to the boundaries
	#Thus we will solve all 4 stellar structure equations as an eigenvalue problem
	
	#Derivatives of radius differences
	d_deltaR_dP_c=(Rdiff1-R_diff[k])/(0.01*P_c)
	d_deltaR_dT_c=(Rdiff2-R_diff[k])/(0.01*T_c)
	d_deltaR_dLs=(Rdiff3-R_diff[k])/(0.01*L_s)
	d_deltaR_dRs=(Rdiff4-R_diff[k])/(0.01*R_s)

	#Derivatives of pressure differences
	d_deltaP_dP_c=(Pdiff1-P_diff[k])/(0.01*P_c)
	d_deltaP_dT_c=(Pdiff2-P_diff[k])/(0.01*T_c)
	d_deltaP_dLs=(Pdiff3-P_diff[k])/(0.01*L_s)
	d_deltaP_dRs=(Pdiff4-P_diff[k])/(0.01*R_s)

	#Derivatives of temperature differences
	d_deltaT_dP_c=(Tdiff1-T_diff[k])/(0.01*P_c)
	d_deltaT_dT_c=(Tdiff2-T_diff[k])/(0.01*T_c)
	d_deltaT_dLs=(Tdiff3-T_diff[k])/(0.01*L_s)
	d_deltaT_dRs=(Tdiff4-T_diff[k])/(0.01*R_s)

	#Derivatives of luminosity differences
	d_deltaL_dP_c=(Ldiff1-L_diff[k])/(0.01*P_c)
	d_deltaL_dT_c=(Ldiff2-L_diff[k])/(0.01*T_c)
	d_deltaL_dLs=(Ldiff3-L_diff[k])/(0.01*L_s)
	d_deltaL_dRs=(Ldiff4-L_diff[k])/(0.01*R_s)

	#The entries for the A matrix 
	A=np.zeros(16).reshape(4,4)
	A[0,0]=d_deltaR_dP_c
	A[0,1]=d_deltaR_dT_c
	A[0,2]=d_deltaR_dLs
	A[0,3]=d_deltaR_dRs
	A[1,0]=d_deltaP_dP_c
	A[1,1]=d_deltaP_dT_c
	A[1,2]=d_deltaP_dLs
	A[1,3]=d_deltaP_dRs
	A[2,0]=d_deltaT_dP_c
	A[2,1]=d_deltaT_dT_c
	A[2,2]=d_deltaT_dLs
	A[2,3]=d_deltaT_dRs
	A[3,0]=d_deltaL_dP_c
	A[3,1]=d_deltaL_dT_c
	A[3,2]=d_deltaL_dLs
	A[3,3]=d_deltaL_dRs
	

	y=np.zeros(4)
	y[0]=0.1*R_diff[k]
	y[1]=0.1*P_diff[k]
	y[2]=0.1*T_diff[k]
	y[3]=0.1*L_diff[k]
	
	x=np.linalg.solve(A,y)
	
	
	P_c=P_c-x[0]
	T_c=T_c-x[1]
	L_s=L_s-x[2]
	R_s=R_s-x[3]

print 'Number of iterations is %s' % int(max(k_values))

print 'Total Mass:  %s kg' % round(MS[dim-1],3)
print 'Radius: %s meters' %round(RS[dim-1],3)
print 'Surface Temp: %s K' % round(TS[dim-1],3)
print 'Surface Pressure: %s' %round(PS[dim-1],3)
print 'Surface Luminosity: %s Watts' %round(LS[dim-1],2)


plt.figure(1,figsize=(8.5,11))


plt.subplot(3,2,1)
plt.plot(R[0:dim/2-1],P[0:dim/2-1],'g-')
plt.plot(RS[dim/2:dim-1],PS[dim/2:dim-1],'r-')
plt.title('Pressure vs Radius')
plt.xlabel('Radius (m)')
plt.ylabel('Pressure (KPa)')


plt.subplot(3,2,2)
plt.plot(R[0:dim/2-1],L[0:dim/2-1],'g-')
plt.plot(RS[dim/2:dim-1],LS[dim/2:dim-1],'r-')
plt.title('Luminosity vs Radius')
plt.xlabel('Radius (m)')
plt.ylabel('Luminosity (Watts)')

plt.subplot(3,2,5)
plt.plot(R[0:dim/2-1],M[0:dim/2-1],'g-')
plt.plot(RS[dim/2:dim-1],MS[dim/2:dim-1],'r-')
plt.title('Mass vs Radius')
plt.xlabel('Radius (m)')
plt.ylabel('Mass (kg)')

plt.subplot(3,2,6)
plt.plot(R[0:dim/2-1],T[0:dim/2-1],'g-')
plt.plot(RS[dim/2:dim-1],TS[dim/2:dim-1],'r-')
plt.title('Temperature vs Radius')
plt.xlabel('Radius (m)')
plt.ylabel('Temperature (K)')

plt.show()


