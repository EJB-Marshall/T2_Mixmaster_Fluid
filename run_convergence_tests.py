import numpy as np
import h5py
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.interpolate import CubicSpline

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



def DxF(f,dx):
    """ This function calculates the first derivative
        using a 2nd order forward finite difference stencil and backward stencil at right boundary"""
    h = dx
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    # #2nd Order
    # df[0:n-2] = (-3/2*f[0:n-2]+2*f[1:n-1]-0.5*f[2:n])/h
    # df[0] = (-3/2*f[0]+2*f[1]-0.5*f[2])/h
    # df[n-2] = (3/2*f[n-2]-2*f[n-3]+0.5*f[n-4])/h
    # df[n-1] = (3/2*f[n-1]-2*f[n-2]+0.5*f[n-3])/h
    

    # #4th Order 
    df[0] = (-25/12*f[0]+4*f[1]-3*f[2]+4/3*f[3]-1/4*f[4])/h
    df[1] = (-25/12*f[1]+4*f[2]-3*f[3]+4/3*f[4]-1/4*f[5])/h
    df[2:n-2] = (1/12*f[0:n-4]-2/3*f[1:n-3]+2/3*f[3:n-1]-1/12*f[4:n])/h
    df[n-2] = (25/12*f[n-2]-4*f[n-3]+3*f[n-4]-4/3*f[n-5]+1/4*f[n-6])/h
    df[n-1] = (25/12*f[n-1]-4*f[n-2]+3*f[n-3]-4/3*f[n-4]+1/4*f[n-5])/h

    return df

def Dx3(f,dx):
    """ This function calculates the third derivative
        using a 4th order central finite difference stencil"""
    h = dx
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    
    #4th Order Central Finite Difference
    df[3:n-3] = (1/8*f[0:n-6]-1*f[1:n-5]+13/8*f[2:n-4]\
                 -13/8*f[4:n-2]+1*f[5:n-1]-1/8*f[6:n])/(h**3)

    #Periodic Boundary Condition
    df[2] = (1/8*f[-1]-1*f[0]+13/8*f[1]\
                 -13/8*f[3]+1*f[4]-1/8*f[5])/(h**3)
    df[1] = (1/8*f[-2]-1*f[-1]+13/8*f[0]\
                 -13/8*f[2]+1*f[3]-1/8*f[4])/(h**3)
    df[0] = (1/8*f[-3]-1*f[-2]+13/8*f[-1]\
                 -13/8*f[1]+1*f[2]-1/8*f[3])/(h**3)
    
    df[n-3] = (1/8*f[n-6]-1*f[n-5]+13/8*f[n-4]\
                 -13/8*f[n-2]+1*f[n-1]-1/8*f[0])/(h**3)
    df[n-2] = (1/8*f[n-5]-1*f[n-4]+13/8*f[n-3]\
                 -13/8*f[n-1]+1*f[0]-1/8*f[1])/(h**3)
    df[n-1] = (1/8*f[n-4]-1*f[n-3]+13/8*f[n-2]\
                 -13/8*f[0]+1*f[1]-1/8*f[2])/(h**3)
    
    return df

def Dx4(f,dx):
    """ Computes periodic 4th derivative for artificial viscosity terms
      Inputs: f - function to be differentiated
                r - Grid instance"""
    
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    h = dx
    df[2:n-2] = (1*f[0:n-4]-4*f[1:n-3]+6*f[2:n-2]-4*f[3:n-1]+1*f[4:n])/h**4

    df[1] = (1*f[n-1]-4*f[0]+6*f[1]-4*f[2]+1*f[3])/h**4
    df[0] = (1*f[n-2]-4*f[n-1]+6*f[0]-4*f[1]+1*f[2])/h**4
    df[n-2] = (1*f[n-4]-4*f[n-3]+6*f[n-2]-4*f[n-1]+1*f[0])/h**4
    df[n-1] = (1*f[n-3]-4*f[n-2]+6*f[n-1]-4*f[0]+1*f[1])/h**4

    return df

def DxF2(f,dx):
    """ This function calculates the first derivative
        using a 2nd order forward finite difference stencil and backward stencil at right boundary"""
    h = dx
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    # #2nd Order
    # df[0:n-2] = (-3/2*f[0:n-2]+2*f[1:n-1]-0.5*f[2:n])/h
    # df[0] = (-3/2*f[0]+2*f[1]-0.5*f[2])/h
    # df[n-2] = (3/2*f[n-2]-2*f[n-3]+0.5*f[n-4])/h
    # df[n-1] = (3/2*f[n-1]-2*f[n-2]+0.5*f[n-3])/h
    
    # # #4th Order 
    df[0:n-4] = (-25/12*f[0:n-4]+4*f[1:n-3]-3*f[2:n-2]+4/3*f[3:n-1]-1/4*f[4:n])/h
    df[n-4] = (25/12*f[n-4]-4*f[n-5]+3*f[n-6]-4/3*f[n-7]+1/4*f[n-8])/h
    df[n-3] = (25/12*f[n-3]-4*f[n-4]+3*f[n-5]-4/3*f[n-6]+1/4*f[n-7])/h
    df[n-2] = (25/12*f[n-2]-4*f[n-3]+3*f[n-4]-4/3*f[n-5]+1/4*f[n-6])/h
    df[n-1] = (25/12*f[n-1]-4*f[n-2]+3*f[n-3]-4/3*f[n-4]+1/4*f[n-5])/h

    return df


############################################################################
#Evaluate constraint quantity and calculate L2 norm functions
############################################################################

def computel2norm(f,dx):
    """This function outputs a vector whose 
    ith entry is the L2 norm at timestep i"""
    l2norm = np.zeros_like(f[:,1])
    for i in range(0,np.shape(f)[0]):
        l2norm[i] = np.sqrt(dx*np.sum(f[i,:]**2))
    return l2norm

def computel1norm(f,dx):
    """This function outputs a vector whose 
    ith entry is the L2 norm at timestep i"""
    l1norm = np.zeros_like(f[:,1])
    for i in range(0,np.shape(f)[0]):
        l1norm[i] = dx*np.sum(np.abs(f[i,:]))
    return l1norm

def compute_integral_abs(f,dx):
    """This function outputs a vector whose 
    ith entry is the absolute value of spatial integral at timestep i"""
    int_eval = np.zeros_like(f[:,1])
    for i in range(0,np.shape(f)[0]):
        int_eval[i] = np.abs(dx*np.sum(f[i,:]))
    return int_eval

def computeH3norm(f,dx):
    """This function outputs a vector whose 
    ith entry is the H3 norm at timestep i"""
    l2norm = np.zeros_like(f[:,1])
    for i in range(0,np.shape(f)[0]):
        l2norm[i] = np.sqrt(dx*np.sum(Dx3(f[i,:],dx)**2))
    return l2norm



def constrainteval(soln,dx,t):
    """This function evaluates the constraints"""

    Sigma_Minus, Sigma_Times, N_Minus, N_Times, \
    E11, Sigma2, Sigma3, LambdaTilde,\
    tau, S1, S2, S3 = soln
    
    CM1 = np.zeros_like(Sigma_Minus)
    CM2 = np.zeros_like(Sigma_Minus)
    CBeta = np.zeros_like(Sigma_Minus)

    for i in range(0,np.shape(t)[0]):

        r = 3*Sigma_Times[i,:]*N_Minus[i,:] - 3*N_Times[i,:]*Sigma_Minus[i,:] - 3/2*S1[i,:]

        CM1[i,:] = E11[i,:]*Dx(Sigma3[i,:],dx) - r*Sigma3[i,:] - np.sqrt(3)*Sigma3[i,:]*N_Times[i,:] + np.sqrt(3)*S2[i,:]

        CM2[i,:] = E11[i,:]*Dx(Sigma2[i,:],dx) - r*Sigma2[i,:] - np.sqrt(3)*Sigma2[i,:]*N_Times[i,:] + np.sqrt(3)*S3[i,:] \
                    + 2*np.sqrt(3)*Sigma3[i,:]*N_Minus[i,:]

        CBeta[i,:] = E11[i,:]*Dx(LambdaTilde[i,:],dx) - 2*LambdaTilde[i,:]*r
        
      
    return np.array([CM1,CM2,CBeta])


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



def restrict_grid(soln,num_grids):
    """ In a cell-centred scheme, the cell centres will not be aligned
        for different resolutions. Thus, we must average the cell-centred values from
         high resolution runs to align with the low resolution cells. """
    
    n = soln.shape[0] # Number of grid points for high res run
    averaged_runs = []

    averaged_runs.append(soln)

    for k in range(1,num_grids+1):

        averaged_runs.append(np.zeros((int(n/(2**(k))))))

        for parent_index in range(int(n/(2**(k)))):
            child_index = 2 * parent_index

            ### Conservative update
            averaged_runs[k][parent_index] = 0.5*(averaged_runs[k-1][child_index] + averaged_runs[k-1][child_index+1])

    return averaged_runs



############################################################################
# Recover Primitive Variables
############################################################################
constraint_list = []
t_list = []
grid_list = []
convergence_list = []
norm_list = []
test_list = []
Output_1 = h5py.File('/Users/elliotmarshall/Desktop/T2_Fluid_Fortran/T2_Primitive_Vars/HDF_Files/500_K05_Convergence2.hdf5', 'r')
t_original = Output_1['Time'][:]


# file_names = ['200_FutureFluid.hdf5','400_FutureFluid.hdf5','800_FutureFluid.hdf5','1600_FutureFluid.hdf5','3200_FutureFluid.hdf5']#,'6400_FutureFluid.hdf5']
# file_names = ['400_FutureFluid.hdf5','800_FutureFluid.hdf5','1600_FutureFluid.hdf5','3200_FutureFluid.hdf5']
# file_names = ['200_Exact_OT.hdf5','400_Exact_OT.hdf5','800_Exact_OT.hdf5','1600_Exact_OT.hdf5','3200_Exact_OT.hdf5','6400_Exact_OT.hdf5']
# file_names = ['400_K05_Convergence.hdf5','800_K05_Convergence.hdf5','1600_K05_Convergence.hdf5','3200_K05_Convergence.hdf5']
file_names = ['500_K05_Convergence2.hdf5','1000_K05_Convergence2.hdf5','2000_K05_Convergence2.hdf5','4000_K05_Convergence2.hdf5']

for k in range(len(file_names)):
    Output_1 = h5py.File('/Users/elliotmarshall/Desktop/T2_Fluid_Fortran/T2_Primitive_Vars/HDF_Files/'+file_names[k], 'r')
    t = Output_1['Time'][:]
    grid = Output_1['x_coordinates'][:]
    K = Output_1['K'][()]

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


    dx = np.abs(grid[1]-grid[0])

    time = 8000
    print(t_original[time])

    # print(Sigma2[0,:])


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


    soln = np.array([Sigma_Minus, Sigma_Times, N_Minus, N_Times, \
    E11, Sigma2, Sigma3, LambdaTilde,\
    T_00, T_01, T_02, T_03])

    CM1,CM2,CBeta = constrainteval(soln,dx,t)        
    converge_time = np.zeros_like(Sigma_Times[0,:])
    # for j in range(np.shape(N_Times)[1]):
    #     if k>0:
    #         spl = CubicSpline(t,nu1[:,j])
    #         converge_time[j] = spl(t_original[time])

    for j in range(np.shape(N_Times)[1]):
        if k>0:
            spl = CubicSpline(np.flip(t),np.flip(Sigma_Minus[:,j]))
            converge_time[j] = spl(t_original[time])




    # totalH3 = computeH3norm(nu1,dx)**2 + computeH3norm(mu,dx)**2  
    # totall2 = computel2norm(nu1,dx)**2 + computel2norm(mu,dx)**2 
    # norm_ratio = totalH3/totall2

    l2norm = computel2norm(CM1,dx)
    # norm_spl = CubicSpline(t,norm_ratio)
    t_list.append(t)
    grid_list.append(grid)
    test_list.append(Sigma_Minus)
    convergence_list.append(converge_time)
    norm_list.append(l2norm)




### Compute restricted high resolution run
avg_high_res = restrict_grid(convergence_list[3],3)




plt.rcParams.update({"text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern",
    "savefig.bbox": "tight",
    "savefig.format": "pdf"})
plt.rc('font', size=16)





# plt.figure()
plt.plot(t_list[0],np.log2(norm_list[0]),label=r'$500$')
plt.plot(t_list[1],np.log2(norm_list[1]),label=r'$1000$')
plt.plot(t_list[2],np.log2(norm_list[2]),label=r'$2000$')
plt.plot(t_list[3],np.log2(norm_list[3]),label=r'$4000$')
# plt.plot(t_list[3],np.max(np.log(np.abs(E11)),axis=1),label='E11')
# plt.plot(t_list[3],np.max(tau_log,axis=1),label='T^00')
# plt.plot(t_list[4],np.log2(norm_list[4]),label=r'$3200$')
plt.xlabel(r'$t$')
plt.ylabel(r'$\log_{2}\|C\|_{2}$')
plt.legend()
# plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
# plt.yscale("log")
# plt.ylim([10**(-15),1])
# print(t_original[time])
plt.show()


# plt.figure()
# plt.plot(grid_list[0],np.log2(np.abs(test_list[0][time,:])),label='400')
# plt.plot(grid_list[1],np.log2(np.abs(convergence_list[1][:])),label='800')
# plt.plot(grid_list[2],np.log2(np.abs(convergence_list[2][:])),label='1600')
# # plt.plot(grid_list[3],np.log2(np.abs(convergence_list[3][:])),label='1600')
# # plt.plot(grid_list[4],np.log2(np.abs(convergence_list[4][:])),label='3200')
# # plt.plot(grid_list[5],convergence_list[5][:],label='3200')
# plt.xlabel(r'$x$')
# plt.ylabel(r'$\log_{2}|\Delta|$')
# plt.legend()
# plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
# plt.tight_layout()
# plt.show()



plt.figure(1)
plt.plot(grid_list[0],np.log2(np.abs(test_list[0][time,:]-avg_high_res[-1][:])),label='500')
plt.plot(grid_list[1],np.log2(np.abs(convergence_list[1][:]-avg_high_res[-2][:])),label='1000')
plt.plot(grid_list[2],np.log2(np.abs(convergence_list[2][:]-avg_high_res[-3][:])),label='2000')
# plt.plot(grid_list[3],np.log2(np.abs(convergence_list[3][:]-avg_high_res[-4][:])),label='1600')
# plt.plot(grid_list[4],np.log2(np.abs(convergence_list[4][:]-avg_high_res[-5][:])),label='3200')
# plt.plot(grid_list[0],np.log2(np.abs(test_list[0][time,:]-test_list[5][time*32,::32])),label='200')
# plt.plot(grid_list[1],np.log2((np.abs(test_list[1][time*2,:]-test_list[5][time*32,::16]))),label='400')
# plt.plot(grid_list[2],np.log2((np.abs(test_list[2][time*4,:]-test_list[5][time*32,::8]))),label='800')
# plt.plot(grid_list[3],np.log2((np.abs(test_list[3][time*8,:]-test_list[5][time*32,::4]))),label='1600')
# plt.plot(grid_list[4],np.log2((np.abs(test_list[4][time*16,:]-test_list[5][time*32,::2]))),label='3200')
# plt.xlabel(r'$x$')
plt.ylabel(r'$\log_{2}|\Delta|$')
plt.legend()
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
# print(t_original[time])
plt.show()



### Exact OT Solution Convergence Plot

# plt.figure(2)
# plt.plot(grid_list[0],np.log2(np.abs(test_list[0][time,:]-1/(np.sqrt(3))*np.exp(t_list[0][time])*np.sin(np.exp(2*t_list[0][time])-2*grid_list[0]))),label='200')
# plt.plot(grid_list[1],np.log2(np.abs(convergence_list[1][:]-1/(np.sqrt(3))*np.exp(t_list[0][time])*np.sin(np.exp(2*t_list[0][time])-2*grid_list[1]))),label='400')
# plt.plot(grid_list[2],np.log2(np.abs(convergence_list[2][:]-1/(np.sqrt(3))*np.exp(t_list[0][time])*np.sin(np.exp(2*t_list[0][time])-2*grid_list[2]))),label='800')
# plt.plot(grid_list[3],np.log2(np.abs(convergence_list[3][:]-1/(np.sqrt(3))*np.exp(t_list[0][time])*np.sin(np.exp(2*t_list[0][time])-2*grid_list[3]))),label='1600')
# plt.plot(grid_list[4],np.log2(np.abs(convergence_list[4][:]-1/(np.sqrt(3))*np.exp(t_list[0][time])*np.sin(np.exp(2*t_list[0][time])-2*grid_list[4]))),label='3200')
# plt.plot(grid_list[5],np.log2(np.abs(convergence_list[5][:]-1/(np.sqrt(3))*np.exp(t_list[0][time])*np.sin(np.exp(2*t_list[0][time])-2*grid_list[5]))),label='6400')
# # # plt.plot(grid_list[0],np.log2(np.abs(test_list[0][time,:]-test_list[5][time*32,::32])),label='200')
# # # plt.plot(grid_list[1],np.log2((np.abs(test_list[1][time*2,:]-test_list[5][time*32,::16]))),label='400')
# # # plt.plot(grid_list[2],np.log2((np.abs(test_list[2][time*4,:]-test_list[5][time*32,::8]))),label='800')
# # # plt.plot(grid_list[3],np.log2((np.abs(test_list[3][time*8,:]-test_list[5][time*32,::4]))),label='1600')
# # # plt.plot(grid_list[4],np.log2((np.abs(test_list[4][time*16,:]-test_list[5][time*32,::2]))),label='3200')
# plt.xlabel(r'$x$')
# plt.ylabel(r'$\log_{2}|\Delta|$')
# plt.legend()
# plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
# plt.tight_layout()
# print(t_original[time])
# plt.show()
