#solves the EDQNM equation for active fluids in 2D https://arxiv.org/abs/1711.10865
#April 25, 2018 (MJ)
 
def export_and_plot(spectrum,wavenumber):#function to plot energy spectrum
 ax.plot(wavenumber,spectrum)
 plt.pause(0.001)
 ax.clear()
 
def find_eta(m,spectrum):
 sum0=0.0
 for i in range(m):
  sum0=sum0+i*i*spectrum[i]
 return lambda0*gamma*np.sqrt(sum0*scale*scale*scale)+np.fabs((1-wavenumber[m]**2)**2+alpha)
 
def find_x(k,p,q):
 return (1.0*p*p+q*q-k*k)/(2*p*q)

def find_y(k,p,q):
 return (1.0*k*k+q*q-p*p)/(2*k*q)

def find_z(k,p,q):
 return (1.0*k*k+p*p-q*q)/(2*k*p)

def find_transfer_term(spectrum,transfer_term,wavenumber):#function to evaluate the transfer term using EDQNM closure theory
 eta=np.zeros(N)
 start = time.time()
 eta[0]=0.001
 for k in range(1,N):
  eta[k]=find_eta(k,spectrum)
 end = time.time()
 print "Took %f ms for eta" % ((end - start) * 1000.0) 
 E=spectrum
 start = time.time() 
 T=one*E.reshape(1,1,N)*(C*E.reshape(1,N,1)-D*E.reshape(N,1,1))/(one*eta.reshape(N,1,1)+one*eta.reshape(1,N,1)+one*eta.reshape(1,1,N))
 end = time.time()
 print "Took %f ms for eval" % ((end - start) * 1000.0)
 start = time.time() 
 transfer_term=T.sum(axis=(1,2))*scale*scale*(-4/np.pi)*lambda0*lambda0
 end = time.time()
 print "Took %f ms for truncation" % ((end - start) * 1000.0)
 return transfer_term
 
####function MAIN()####
import numpy as np
import matplotlib.pyplot as plt
import time

global scale#spacing between wavenumbers dk
global N #number of grid points
global end_time 
global dt
global alpha
global beta
global lambda0
global gamma
global wavenumber
global C
global D

scale=0.025
N=256
end_time=1000
dt=0.025
alpha=-0.8
beta=0.01
lambda0=3.5
gamma=0.55

fig,ax=plt.subplots()
plt.ion()

wavenumber=np.arange(N)*scale
wavenumber[0]=0.001

one=np.ones((N,N,N))

X=np.zeros((N,N,N)) 
Y=np.zeros((N,N,N)) 
Z=np.zeros((N,N,N))  

for k in range(1,N):
  for p in range(1,N):
   for q in range(1,N):
    if p+q>k and q+k>p and p+k>q:
     X[k,p,q]=find_x(k,p,q)
     Y[k,p,q]=find_y(k,p,q)
     Z[k,p,q]=find_z(k,p,q)
     
C=one*wavenumber.reshape(N,1,1)*(X*Y-Z+2*Z*Z*Z)/(one*wavenumber.reshape(1,1,N)*np.sqrt(one-X*X))
D=C*wavenumber.reshape(1,N,1)
C=C*wavenumber.reshape(N,1,1)

time0=0

spectrum=wavenumber**2*np.exp(-wavenumber**2)#initialize spectrum
spectrum[0]=0.0

spectrum0=wavenumber**2*np.exp(-wavenumber**2)#initialize spectrum
spectrum0[0]=0.0

transfer_term=0*np.random.rand(N)/10#initialize transfer term
transfer_term[0]=0.0

linear_factor=2*(-(1-wavenumber**2)**2-alpha)#linear factor in the evolution equation for energy speectrum
linear_factor=np.exp(linear_factor*dt/2)

f=open('energy_opt_new2.dat','w')

while time0<end_time:#time evolution using RK2 scheme
 start0 = time.time()
 print (time0)
 #export_and_plot(transfer_term,wavenumber)
 energy=np.sum(spectrum)*scale
 transfer_term=find_transfer_term(spectrum,transfer_term,wavenumber)
 spectrum0=linear_factor*np.exp(-8*beta*energy*dt/2)*(spectrum-0.5*dt*transfer_term)
 energy=np.sum(spectrum0)*scale
 transfer_term=find_transfer_term(spectrum0,transfer_term,wavenumber)
 spectrum=linear_factor*linear_factor*np.exp(-8*beta*energy*dt)*spectrum-dt*linear_factor*np.exp(-8*beta*energy*dt/2)*transfer_term
 f.write("%E\n"%energy)
 f.flush()
 np.savetxt("spectrum_dt_opt_new2.dat",spectrum)
 np.savetxt("transfer_term_dt_opt_new2.dat",transfer_term)
 time0=time0+dt
 end0 = time.time()
 print "Took %f ms for tot" % ((end0 - start0) * 1000.0)
####end of program####
