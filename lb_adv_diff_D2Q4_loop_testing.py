# Author : Guy Jérémie
# Date : 
# inspired by Jonas Latt

from numpy import *
from functions import *
import time

########################### Flow & Topology Definition #############################################

# Max iterations (dt =1)
maxIter = 35000                 

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
d = 1                           # Clot dissolving rate
rhoTPA_initial = 1              # tPA added concentration

########################## Lattice Constants ###########################################

# Velocity directions vectors, D2Q9
v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
            [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])

# Directionnal weights, D2Q9
w = array([1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

# Velocity directions vector for tPA, D2Q4
v_tPA = array([[ 1, 0], [ 0, 1], [ 0, -1], [ -1, 0]])

# Directionnal weighta for tPA, D2Q4
w_tPA = array([1/4, 1/4, 1/4, 1/4])

# Sound velocity adapted to lattice units   
cs2 = 1/3                       

#################### Main Function Definitions ######################################

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

# Macroscopic variable computation for tPA (only rho)
def macroscopicTPA(tPAin):
    rhoTPA = sum(tPAin, axis=0)
    return rhoTPA

# Equilibrium distribution function
def equilibriumTPA(rhoTPA, u):
    tPAeq = zeros((4,nx,ny))
    for i in range(4):
        vu = v_tPA[i,0]*u[0,:,:] + v_tPA[i,1]*u[1,:,:]
        tPAeq[i,:,:] = w_tPA[i]*rhoTPA*(1 + (1/cs2)*vu)
    return tPAeq

# Acceleration field force definition and porous region counter-field
def addResistingClotForce(rho, u, F, K):
    FF = zeros((9,nx,ny))
    for i in range(9):
        FF[i,:,:] = rho*(v[i,0]*(F[0,:,:] - K[0,:,:]*u[0,:,:]) + v[i,1]*(F[1,:,:] - K[1,:,:]*u[1,:,:]))
        FF[i,:,:] = FF[i,:,:] * (w[i] / cs2)
    return FF

# Dissolving the clot proportionnaly to tPA amount
def bindAndDissolve(tPAin, K):
    #Compute the total amount of population in tPAin
    sumTPAin = sum(tPAin, axis=0)

    # Compute the dissolution amount (not : k[0] and K[1] change in the same way)
    dissolutionAmount = sumTPAin*K[0,:,:]*d
    
    # Check if each site has had dissolution ...
    hasDissolved = dissolutionAmount == 0

    # TESTVALUE
    # if 49500<execTime<50500 :
    #     saveDissolutionAmount(tPAClotFile, displayTEST, dissolutionAmount)

    # ... and remove tPAin accordingly
    for i in range(4):
        tPAin[i,:,:] = where(hasDissolved, tPAin[i,:,:], 0)

    # Dissolving the clot : 
    # 1. Start by applying dissolution
    K_tmp0 = K[0,:,:] - dissolutionAmount
    K_tmp1 = K[1,:,:] - dissolutionAmount

    # 2. Check if dissolution amount exceeds K (resulting in a K < 0)
    isNotNegative0 = K_tmp0 >= 0
    isNotNegative1 = K_tmp1 >= 0

    # 3. If it does, set K to 0 (= no clot, all dissolved) 
    #    and if not, keep K value
    positive_K_tmp0 = where(isNotNegative0, K_tmp0, 0)
    positive_K_tmp1 = where(isNotNegative1, K_tmp1, 0)

    # 4. Check if K<1e-7, we consider it = 0
    isTooSmall0 = positive_K_tmp0 > 1e-7 
    isTooSmall1 = positive_K_tmp1 > 1e-7

    # 5. Adjust accordingly
    K[0,:,:] = where(isTooSmall0, positive_K_tmp0, 0)
    K[1,:,:] = where(isTooSmall1, positive_K_tmp1, 0)

    # Return new tPA population and new clot value
    return tPAin, K

################################## Masks ####################################

# Bounceback nodes mask
bounceback = full((nx,ny),False)
bounceback[:,0] = True                              # Top border
bounceback[:,ny-1] = True                           # Bottom border
bounceback[0,:] = True                              # Left border
bounceback[nx-1,:] = True                           # Right border
bounceback[tubeSize+1:(nx-1-tubeSize),              
           tubeSize+1:(ny-1-tubeSize)] = True       # Obstacle
# bounceback[tubeSize+1:(nx-1-tubeSize),
#            1+2*tubeSize+1:1+3*tubeSize+1] = False   # Branch

# Open path mask
openPath = invert(bounceback)

# Clot mask for clot in the upper tube
clot = full((nx,ny),False)
clot[clotCoord[0]:clotCoord[1]+1,1:1+tubeSize] = True

displayTEST = full((nx,ny),False)
displayTEST[clotCoord[0]-3:clotCoord[1]+1-15,1:1+tubeSize] = True

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
mainDirectory, clotDirectory, velocityDirectory, clotForceDirectory, tPADensityDirectory, tPAFlowDirectory = createRepositories(maxIter,nx,ny,viscosity,rho_initial, rhoTPA_initial, d, F_initial,K_initial)

tPAClotFile = open(mainDirectory + "/monitor_clot_tpa" + ".txt", 'w')

############################# System Initliaization #################################

# Velocity initialisation
vel = zeros((2,nx,ny))

# Density initialisation
rho = full((nx,ny),rho_initial)

# Initialisation of the populations at equilibrium with the given density & velocity.
fin = equilibrium(rho, vel)
fout = equilibrium(rho, vel)

# Loading already converged variables for faster execution time
fin, fout, _, u = getVariables("loop_Kxy", nx, ny, viscosity, rho_initial, F_initial, K_initial, 65000)
# _, _, _, u = getVariables("loop", nx, ny, viscosity, rho_initial, F_initial, K_initial, 65000)

# tPA density initialisation
rhoTPA = zeros((nx,ny))
rhoTPA[1:tubeSize+1,ny//2] = rhoTPA_initial

# tPA population initialisation
tPAin = equilibriumTPA(rhoTPA, u)
tPAout = equilibriumTPA(rhoTPA, u)

################################# Main time loop ######################################

# Monitoring execution time
start_time = time.time()

# main loop
for execTime in range(maxIter):

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)       # fluid 
    rhoTPA = macroscopicTPA(tPAin)  # tPA 

    # injecting tPA constantly
    rhoTPA[1:tubeSize+1,ny//2] = rhoTPA_initial

    # Compute equilibrium.
    feq = equilibrium(rho, u)           # fluid
    tPAeq = equilibriumTPA(rhoTPA, u)   # tPA

    # BGK collision step for open path
    fout[:,openPath] = fin[:,openPath] - omega * (fin[:,openPath] - feq[:,openPath])            # fluid
    tPAout[:,openPath] = tPAin[:,openPath] - omega * (tPAin[:,openPath] - tPAeq[:,openPath])    # tPA

    # Bounce-back condition 
    for i in range(9):                                  # fluid
        fout[i, bounceback] = fin[8-i, bounceback]
    for i in range(4):                                  # tPA
        tPAout[i, bounceback] = tPAin[3-i, bounceback]
    

    # Force field injection
    fout += addResistingClotForce(rho, u, F, K)
    
    # Streaming step for fluid in every direction i=0:8
    fin[0,:,:] = roll(roll(fout[0,:,:],1,axis=0),1,axis=1)      # i = 0
    fin[1,:,:] = roll(fout[1,:,:],1,axis=0)                     # i = 1
    fin[2,:,:] = roll(roll(fout[2,:,:],1,axis=0),-1,axis=1)     # i = 2
    fin[3,:,:] = roll(fout[3,:,:],1,axis=1)                     # i = 3
    fin[4,:,:] = fout[4,:,:]                                    # i = 4
    fin[5,:,:] = roll(fout[5,:,:],-1,axis=1)                    # i = 5
    fin[6,:,:] = roll(roll(fout[6,:,:],-1,axis=0),1,axis=1)     # i = 6
    fin[7,:,:] = roll(fout[7,:,:],-1,axis=0)                    # i = 7
    fin[8,:,:] = roll(roll(fout[8,:,:],-1,axis=0),-1,axis=1)    # i = 8

    # Streaming step for tPA in every direction i=0:4
    tPAin[0,:,:] = roll(tPAout[0,:,:],1,axis=0)                     # i = 0
    tPAin[1,:,:] = roll(tPAout[1,:,:],1,axis=1)                     # i = 1
    tPAin[2,:,:] = roll(tPAout[2,:,:],-1,axis=1)                    # i = 2
    tPAin[3,:,:] = roll(tPAout[3,:,:],-1,axis=0)                    # i = 3


    
    # if 49500<execTime<50500 :
    #     SaveTPAValues(tPAClotFile, True, tPAin, tPAout, K, displayTEST, execTime)
    
    
    # Dissolve clot
    tPAin, K = bindAndDissolve(tPAin, K)

    # if 49500<execTime<50500 :
    #     SaveTPAValues(tPAClotFile, False, tPAin, tPAout, K, displayTEST, execTime)

    if (execTime%250==0):
        saveTPADensity(tPADensityDirectory, rhoTPA, execTime)
        showClotForce(clotForceDirectory, K_initial, K, clot, execTime)
        saveTPAFlow(tPAFlowDirectory, rhoTPA, u, execTime)

    # Velocity visualisation
    # if (execTime%10==0):
        # visualise(u)
        # visualiseTPA(rhoTPA)
        # showClotForce(clotForceDirectory, K_initial, K, clot, execTime)
    
    # Clot velocity and velocity profiles graphics generation
    # if(execTime%500==0):
        
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
tPAClotFile.close()


########################### Converged System Saving ############################# 

# Saving converged system to load directly at next run
# saveVariables("loop_Kxy", nx, ny, viscosity, rho_initial, F_initial, K_initial, maxIter, fin, fout, rho, u)