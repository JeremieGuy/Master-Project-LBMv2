# Author : Guy Jérémie
# Date : 
# inspired by Jonas Latt

from numpy import *
from functions import *
import time

########################### Flow & Topology Definition #############################################

# Max iterations (dt =1)
maxIter = 900                 

# System topology definition
nx, ny = 260, 200               # Number of lattice nodes (dx = 1)
tubeSize = 21                   # Diameters of the tubes in the system

# Fluid variables
viscosity = 0.01                # Kinematic viscosity
omega = 1 / (3*viscosity+0.5);  # Relaxation parameter
rho_initial = 2.5               # Inital density of the system

# Forces
F_initial = [0,-0.0001]         # Accelerating force F[nx,ny]
K_initial = [0.001,0.001]       # Initial resisting force of porous region K[2,nx,ny]

# Clot
clotSize = 20                   # Size of the clot lenghtwise in a tube section
clotCoord = [nx//2-clotSize//2, 
             nx//2+clotSize//2] # Clot coordinates
d = 2                           # Clot dissolving rate

########################## Lattice Constants ###########################################

# Velocity directions vectors
v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
            [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])

# Directionnal weights
w = array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

# Sound velocity adapted to lattice units   
cs2 = 1/3                       

#################### Main Function Definitions ######################################

# Macroscopic variable computation
# def macroscopic(fin):
#     # density
#     rho = zeros((nx,ny))
#     rho[openPath] = sum(fin[:,openPath], axis=0)
#     # velocity
#     u = zeros((2, nx, ny))
#     for i in range(9):
#         u[0,openPath] += v[i,0] * fin[i,openPath]
#         u[1,openPath] += v[i,1] * fin[i,openPath]
#     u[:,openPath] /= rho[openPath]
#     return rho, u

# # Equilibrium distribution function (rho = array)
# def equilibrium(rho, u):          
#     usqr = 3/2 * (u[0,openPath]**2 + u[1,openPath]**2)
#     feq = zeros((9,nx,ny))
#     for i in range(9):
#         cu = 3 * (v[i,0]*u[0,openPath] + v[i,1]*u[1,openPath])
#         feq[i,openPath] = rho[openPath]*w[i] * (1 + cu + 0.5*cu**2 - usqr)
#     return feq

# Macroscopic variable computation
def macroscopic(fin):
    rho = sum(fin, axis=0)
    u = zeros((2, nx, ny))
    for i in range(9):
        u[0,:,:] += v[i,0] * fin[i,:,:]
        u[1,:,:] += v[i,1] * fin[i,:,:]
    u /= rho
    return rho, u

# Equilibrium distribution function
def equilibrium(rho, u): 
    usqr = 3/2 * (u[0]**2 + u[1]**2)
    feq = zeros((9,nx,ny))
    for i in range(9):
        cu = 3 * (v[i,0]*u[0,:,:] + v[i,1]*u[1,:,:])
        feq[i,:,:] = rho*w[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq

# Acceleration field force definition and porous region counter-field
def addResistingClotForce(rho, u, F, K):
    FF = zeros((9,nx,ny))
    for i in range(9):
        FF[i,:,:] = rho*(v[i,0]*(F[0,:,:] - K[0,:,:]*u[0,:,:]) + v[i,1]*F[1,:,:])
        FF[i,:,:] = FF[i,:,:] * (w[i] / cs2)
    return FF

# Dissolving the clot proportionnaly to the flow
def dissolveClot(rho, u, K, clot):

    # Determining norm of the velocity
    u_norm = sqrt(u[0]**2+u[1]**2)

    # Computing flow
    flow = rho*u_norm
    
    # Computing amount of dissolution 
    dissolution1 = d*K[0,clot]*flow[clot]
    dissolution2 = d*K[1,clot]*flow[clot]

    # Calculating the resulting K values
    K_tmp1 = K[0,clot] - dissolution1
    K_tmp2 = K[0,clot] - dissolution2

    # Assessing if dissolution is > 0 -> if not set value to 0
    condition1 = K_tmp1 > 0  # Create a boolean mask
    condition2 = K_tmp2 > 0
    filtered_K_tmp1 = where(condition1, K_tmp1, 0) # Check if value validates condition
    filtered_K_tmp2 = where(condition2, K_tmp2, 0) # else replace by 0

    # Updating the dissolved K
    K[0, clot] =  filtered_K_tmp1
    K[1, clot] =  filtered_K_tmp2

    return K

################################## Masks ####################################

# Bounceback nodes mask
bounceback = full((nx,ny),False)
bounceback[:,0] = True                              # Top border
bounceback[:,ny-1] = True                           # Bottom border
bounceback[0,:] = True                              # Left border
bounceback[nx-1,:] = True                           # Right border
bounceback[tubeSize+1:(nx-1-tubeSize),              
           tubeSize+1:(ny-1-tubeSize)] = True       # Obstacle
bounceback[tubeSize+1:(nx-1-tubeSize),
           1+2*tubeSize+1:1+3*tubeSize+1] = False   # Branch

# Open path mask
openPath = invert(bounceback)

# Clot mask for clot in the upper tube
clot = full((nx,ny),False)
clot[clotCoord[0]:clotCoord[1]+1,1:1+tubeSize] = True

# Force array resistance mask for porous region
K = zeros((2,nx,ny))
K[0,clot] = K_initial[0]
K[1,clot] = K_initial[1]

# Pulsing field for fluid aceleration in the lower left tube section
pulseField = full((nx,ny),False)
pulseField[1:tubeSize+1,(ny//2-10):(ny//2+11)] = True

# Accelerating force mask
F = zeros((2,nx,ny))
F[0,pulseField] = F_initial[0]
F[1,pulseField] = F_initial[1]

##################### Initialising Output Monitoring Functions #####################

# Generating working directories
mainDirectory, clotDirectory, velocityDirectory, clotForceDirectory = createRepositories(maxIter,nx,ny,viscosity,rho_initial,F_initial,K_initial)

############################# System Initliaization #################################

# Velocity initialisation
vel = zeros((2,nx,ny))

# Density initialisation
rho = full((nx,ny),rho_initial)

# Initialisation of the populations at equilibrium with the given density & velocity.
fin = equilibrium(rho, vel)
fout = equilibrium(rho, vel)

# Loading already converged variables for faster execution time
fin, fout, _, _ = getVariables(nx, ny, viscosity, rho_initial, F_initial, K_initial, 60000)

################################# Main time loop ######################################

# Monitoring execution time
start_time = time.time()

# main loop
for execTime in range(maxIter):

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Compute equilibrium.
    feq = equilibrium(rho, u)

    # BGK collision step for open path
    fout[:,openPath] = fin[:,openPath] - omega * (fin[:,openPath] - feq[:,openPath])

    # Bounce-back condition
    for i in range(9):
        fout[i, bounceback] = fin[8-i, bounceback]

    # Force field injection
    fout += addResistingClotForce(rho, u, F, K)
    
    # Streaming step for every direction i
    fin[0,:,:] = roll(roll(fout[0,:,:],1,axis=0),1,axis=1)      # i = 0
    fin[1,:,:] = roll(fout[1,:,:],1,axis=0)                     # i = 1
    fin[2,:,:] = roll(roll(fout[2,:,:],1,axis=0),-1,axis=1)     # i = 2
    fin[3,:,:] = roll(fout[3,:,:],1,axis=1)                     # i = 3
    fin[4,:,:] = fout[4,:,:]                                    # i = 4
    fin[5,:,:] = roll(fout[5,:,:],-1,axis=1)                    # i = 5
    fin[6,:,:] = roll(roll(fout[6,:,:],-1,axis=0),1,axis=1)     # i = 6
    fin[7,:,:] = roll(fout[7,:,:],-1,axis=0)                    # i = 7
    fin[8,:,:] = roll(roll(fout[8,:,:],-1,axis=0),-1,axis=1)    # i = 8

    # Velocity visualisation
    if (execTime%10==0):
        # visualise(u)
        showClotForce(clotForceDirectory, K_initial, K, clot, execTime)
    
    # Clot dissolution step
    K = dissolveClot(rho, u, K, clot)
    
    # Clot velocity and velocity profiles graphics generation
    # if(execTime%1000==0):
    #     plotResults(clotDirectory, nx, tubeSize, clotCoord, u, rho, execTime)
    #     plotVelocityProfiles(velocityDirectory, nx, ny, tubeSize, u, execTime)
    
    # Displaying current progress
    print("iteration : " + str(execTime) + "/" + str(maxIter), end="\r")

# Final execution time
end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

######################## Final Iteration Monitoring ########################## 

# Final Clot velocity and velocity profiles graphics generation
plotResults(clotDirectory, nx, tubeSize, clotCoord, u, rho, maxIter)
plotVelocityProfiles(velocityDirectory, nx, ny, tubeSize, u, maxIter)


########################### Converged System Saving ############################# 

# Saving converged system to load directly at next run
# saveVariables(nx, ny, viscosity, rho_initial, F_initial, K_initial, maxIter, fin, fout, rho, u)