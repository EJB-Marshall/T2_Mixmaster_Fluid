import numpy as np
from matplotlib import pyplot as plt
import h5py


def Dx(f,dx):
    """ This function calculates the first derivative
        using a 2nd or 4th order central finite difference stencil"""
    h = dx
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    # # #2nd Order Central Finite Difference
    df[1:n-1] = (-f[0:n-2]+f[2:n])/(2*h)

    #Periodic Boundary Condition
    df[0] = (-f[n-1]+f[1])/(2*h)
    df[n-1] = (-f[n-2]+f[0])/(2*h)

    # #4th Order Central Finite Difference
    # df[2:n-2] = (f[0:n-4]-8*f[1:n-3]+8*f[3:n-1]-f[4:n])/(12*h)

    # #Periodic Boundary Condition
    # df[1] = (f[n-1]-8*f[0]+8*f[2]-f[3])/(12*h)
    # df[0] = (f[n-2]-8*f[n-1]+8*f[1]-f[2])/(12*h)
    # df[n-2] = (f[n-4]-8*f[n-3]+8*f[n-1]-f[0])/(12*h)
    # df[n-1] = (f[n-3]-8*f[n-2]+8*f[0]-f[1])/(12*h)

    return df



def conserved_to_primitive(cons ,K):
        
        tau, S1, S2, S3 = cons

        Q = ((S1**2 + S2**2 + S3**2))/((K+1)**2*tau**2)

        Gamma2 = (1 - 2*K*(1+K)*Q + np.sqrt(1-4*K*Q))\
                /(2*(1 - (1+K)**2*Q))
        
        mu = tau/((K+1)*Gamma2 - K)

        nu1 = S1/((K+1)*Gamma2*mu)

        nu2 = S2/((K+1)*Gamma2*mu)

        nu3 = S3/((K+1)*Gamma2*mu)
        

        return np.array([nu1,nu2,nu3,mu])


Output_1 = h5py.File('/Users/elliotmarshall/Desktop/T2_Fluid_Fortran/T2_Primitive_Vars/HDF_Files/500_K05.hdf5', 'r')

t = Output_1['Time'][:]
grid = Output_1['x_coordinates'][:]
K = Output_1['K'][()]

dx = np.abs(grid[1]-grid[0])



# Load in solution data
Sigma_Minus = Output_1['Gravitational/Sigma_Minus'][:]
Sigma_Times = Output_1['Gravitational/Sigma_Times'][:]
N_Minus = Output_1['Gravitational/N_Minus'][:]
N_Times = Output_1['Gravitational/N_Times'][:]
E11 = Output_1['Gravitational/E11'][:]
Sigma2 = Output_1['Gravitational/Sigma2'][:]
Sigma3 = Output_1['Gravitational/Sigma3'][:]
LambdaTilde = Output_1['Gravitational/LambdaTilde'][:]
tau_log = Output_1['Fluid/tau_log'][:]
nu1 = Output_1['Fluid/nu1'][:]
nu2 = Output_1['Fluid/nu2'][:]
nu3 = Output_1['Fluid/nu3'][:]


nu_norm = nu1**(2) + nu2**(2)  + nu3**(2)
Gamma2 = 1/(1 - nu_norm)
T_00 = np.exp(tau_log)
T_01 = (K+1)/(1+K*nu_norm)*T_00*nu1
T_02 = (K+1)/(1+K*nu_norm)*T_00*nu2
T_03 = (K+1)/(1+K*nu_norm)*T_00*nu3
T_11 = (K+1)/(1+K*nu_norm)*T_00*nu1*nu1 + K*T_00/(Gamma2*(1+K*nu_norm))
T_12 = (K+1)/(1+K*nu_norm)*T_00*nu1*nu2 
T_13 = (K+1)/(1+K*nu_norm)*T_00*nu1*nu3 
T_22 = (K+1)/(1+K*nu_norm)*T_00*nu2*nu2 + K*T_00/(Gamma2*(1+K*nu_norm))
T_23 = (K+1)/(1+K*nu_norm)*T_00*nu2*nu3 
T_33 = (K+1)/(1+K*nu_norm)*T_00*nu3*nu3 + K*T_00/(Gamma2*(1+K*nu_norm))


Sigma_Plus = 0.5*(1 - Sigma_Times**(2) - Sigma_Minus**(2) - Sigma2**(2) - Sigma3**(2) \
                    - N_Minus**(2) - N_Times**(2) - T_00 - LambdaTilde)


mu_log = tau_log - np.log((K+1)*Gamma2-K)


# print(Sigma_Minus[0,:] - 0.1*np.sin(grid))


### Define Hubble-Normalised Variables and Generalised (?) Kasner Exponents
Sigma_Plus_H = Sigma_Plus/(1-Sigma_Plus)
Sigma_Minus_H = Sigma_Minus/(1-Sigma_Plus)
N_Minus_H = N_Minus/(1-Sigma_Plus)
N_Times_H = N_Times/(1-Sigma_Plus)
Sigma3_H = Sigma3/(1-Sigma_Plus)
Sigma2_H = Sigma2/(1-Sigma_Plus)
Sigma_H_1 = -2*Sigma_Plus_H
Sigma_H_2 = Sigma_Plus_H + np.sqrt(3)*Sigma_Minus_H
Sigma_H_3 =  Sigma_Plus_H - np.sqrt(3)*Sigma_Minus_H
p1 = 1/3*(1-2*Sigma_Plus_H)
p2 = 1/3*(1 + Sigma_Plus_H + np.sqrt(3)*Sigma_Minus_H)
p3 = 1/3*(1 + Sigma_Plus_H - np.sqrt(3)*Sigma_Minus_H)


plt.rcParams.update({"text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern",
    "savefig.bbox": "tight",
    "savefig.format": "pdf"})
plt.rc('font', size=16)


# print(t[100])

# # plt.plot(grid,Dx(mu_log[29000,:],dx),label=r'$\frac{\partial_{x}\rho}{\rho}$')
# plt.tight_layout()
# plt.plot(t,np.max(tau_log,axis=1),label=r'$\log(|\tilde{T}^{00}|)$')
# plt.plot(t,np.max(np.log(np.abs(Sigma2)),axis=1),label=r'$\log(|E^{1}_{1}|)$')
# plt.xlabel(r'$t$')
# # plt.ylabel(r'$\nu^{1}$')
# plt.legend()
# plt.show()



# for i in range(0,np.shape(t)[0],50):

#     # plt.plot(grid,nu1[i,:],label='nu1')
#     # plt.plot(grid[500:],nu1[i,500:],label='nu1')
#     # plt.plot(grid[500:],np.flip(-nu1[i,:500]),label='nu1')
#     # plt.plot(grid,nu2[i,:],label='nu2')
#     # plt.plot(grid,nu3[i,:],label='nu3')
#     plt.plot(grid,np.sqrt(nu1[i,:]**2 + nu2[i,:]**2 + nu3[i,:]**2),label='nu1')
#     # plt.plot(grid,Dx(tau_log[i,:],dx),label='log(T^00)')
#     # plt.plot(grid,Dx(mu_log[i,:],dx),label='log(T^00)')
#     # plt.plot(grid, Sigma_Minus_H[i,:])
#     # plt.xlabel('x')
#     # plt.ylabel(r'$|\nu|$')
#     # plt.ylabel(r'Density Gradient')
# #     plt.plot(grid,Sigma_Minus[i,:])
# #     plt.legend()
#     # print(np.max(nu1[i,:]**2 + nu2[i,:]**2 + nu3[i,:]**2))
#     # if np.any(np.max(nu1[i,:]**2 + nu2[i,:]**2 + nu3[i,:]**2)>=1.0):
#     #     print("Unphysical velocity!")
#     #     print(np.max(nu1[i,:]**2 + nu2[i,:]**2 + nu3[i,:]**2))

#     # plt.ylim([-1.1,1.1])
#     plt.ylim([-0.1,1.1])
#     # plt.tight_layout()
#     plt.title(t[i])
#     plt.draw()
#     plt.pause(0.05)
#     plt.cla()



point = 291
# for i in range(0,np.shape(t)[0],100):

#     plt.plot(t[:i],nu1[:i,point]**2 + nu2[:i,point]**2 + nu3[:i,point]**2,label='|nu|')
#     plt.legend()
#     # print(np.max(nu1[i,:]**2 + nu2[i,:]**2 + nu3[i,:]**2))
#     # if np.any(np.max(nu1[i,:]**2 + nu2[i,:]**2 + nu3[i,:]**2)>=1.0):
#     #     print("Unphysical velocity!")
#     #     print(np.max(nu1[i,:]**2 + nu2[i,:]**2 + nu3[i,:]**2))

#     # plt.ylim([-1.1,1.1])
#     plt.ylim([-0.1,1.1])
#     plt.xlabel('t')
# #     plt.ylabel(r'$|\nu|$')
#     plt.title(t[i])
#     plt.legend()
#     plt.draw()
#     plt.pause(0.05)
#     plt.cla()




fig, ax = plt.subplots(2,1)




### Fluid Variables
# plt.plot(t[:],np.sqrt(nu1[:,point]**2 + nu2[:,point]**2 + nu3[:point]**2),label=r'$|\nu|$')
# ax[0].plot(t[:],1-np.sqrt(nu1[:,point]**2 + nu2[:,point]**2 + nu3[:,point]**2),label='1-|nu|')
ax[0].plot(t[:],np.abs(nu1[:,point]),label=r'$|\nu^{1}|$')
# ax[0].plot(t[:],np.abs(nu2[:,point]),label=r'$|\nu^{2}|$')
# ax[0].plot(t[:],np.abs(nu3[:,point]),label=r'$|\nu^{3}|$')
# plt.plot(t[:],Sigma_Plus[:,point],label='Sigma_plus')


### Kasner Exponents and Sums
# plt.plot(t[:],np.minimum(p1[:,point],np.minimum(p2[:,point],p3[:,point])),label='p_max')
# plt.plot(t[:],0*t[:] + K,label='K')
ax[1].plot(t[:],p1[:,point]-K,label=r'$P_{1}-K$')
# ax[1].plot(t[:],p2[:,point]-K,label=r'$P_{2}-K$')
# ax[1].plot(t[:],p3[:,point]-K,label=r'$P_{3}-K$')
# plt.plot(t[:],p1[:,point]+p2[:,point]+p3[:,point], label='p_sum')
# plt.plot(t[:],p1[:,point]**2 + p2[:,point]**2 + p3[:,point]**2, label='p_sum_square')
# ax[1].plot(t[:],t[:]*0 + (3*K-1) + 3*(1-K)*Sigma_Plus[:,point],label='Kasner-Arc Critical Value?')


### Gravitational Variables (Hubble-Normalised)
# plt.plot(t[:],t[:]*0 + -0.5*(3*(K+1)-4),label='Kasner-Arc Critical Value?')
# ax[1].plot(t[:],t[:]*0 + (3*K-1) + 2*Sigma_Plus_H[:,point],label='Kasner-Arc Critical Value')
# ax[1].plot(t[:],t[:]*0 + (np.sqrt(3)*Sigma_Minus_H[:,point]-3*Sigma_Plus_H[:,point]),label='Kasner-Arc Critical v3')
# plt.plot(t[:],Sigma_Plus_H[:,point],label='Sigma_Plus, Hubble')

### Gravitational Variables (Beta-Normalised)
# plt.plot(t[:],np.max(np.abs(E11),axis=1),label=r'$E^{1}_{1}$')
# plt.plot(t[:],np.max(np.abs(T_00),axis=1),label=r'$T^{00}$')
# plt.plot(t[:],N_Minus[:,point],label='N_minus')
# plt.plot(t[:],N_Times[:,point],label='N_x')
# plt.plot(t[:],Sigma_Minus[:,point],label='Sigma_minus')
# plt.plot(t[:],Sigma_Plus[:,point],label='Sigma_Plus')
# plt.plot(t[:],Sigma_Times[:,point],label='Sigma_x')
# plt.plot(t[:],Sigma2[:,point],label='Sigma_2')
# plt.plot(t[:],Sigma3[:,point],label='Sigma_3')

ax[0].set_yscale("log")
# # ax[0].set_ylim([10**(-300),10**()])
# # ax[0].set_yticks([10**0,10**(-40),10**(-80),10**(-120),10**(-160),10**(-200),10**(-240)], minor=False)
# ax[0].set_yticks([10**0,10**(-54),10**(-118)], minor=False)
ax[0].legend( prop={'size': 14})
ax[1].legend( prop={'size': 14})
ax[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
# # plt.legend()
# # plt.yscale('log')
plt.xlabel(r'$t$')
plt.ylabel(r'$|\nu|$')
# # plt.tight_layout()
plt.show()

