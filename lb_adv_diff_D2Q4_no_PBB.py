# Author : Guy Jérémie
# Date : 
# inspired by Jonas Latt

from numpy import *
from functionsLB import *
from functionsMonitoring import *
import time

########################### Flow & Topology Definition #############################################
# System topology definition
class Lattice:
    maxIter = 46000              # Max iterations (dt =1)            
    nx, ny = 260, 200               # Number of lattice nodes (dx = 1)
    tubeSize = 21                   # Diameters of the tubes in the system

# Fluid variables
class Fluid:
    viscosity = 0.01                # Kinematic viscosity
    omega = 1 / (3*viscosity+0.5);  # Relaxation parameter
    rho_initial = 2.5               # Inital density of the system
    F_initial = [0,-0.0001]         # Accelerating force F[nx,ny]
    
# Forces
class Clot:
    K_initial = [0.001,0.001]       # Initial resisting force of porous region K[2,nx,ny]
    clotSize = 30                   # Size of the clot lenghtwise in a tube section
    coord = [Lattice.nx//2-clotSize//2, Lattice.nx//2+clotSize//2] # Clot coordinates
    d = 1                           # Clot dissolving rate

# tPA 
class TPA:
    rho_initial = 1              # tPA added concentration

########################## Lattice Constants ###########################################

class D2Q9:
    # Velocity directions vectors, D2Q9
    v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
        [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ]) 
    # Directionnal weights, D2Q9
    w = array([1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36]) 
    # Fluid sound velocity adapted to lattice units
    cs2 = 1/3                          

class D2Q4:
    # Velocity directions vector for tPA, D2Q4
    v = array([[ 1, 0], [ 0, 1], [ 0, -1], [ -1, 0]])
    # Directionnal weighta for tPA, D2Q4
    w = array([1/4, 1/4, 1/4, 1/4])
    # tPA sound velocity adapted to lattice units
    cs2 = 1/2                                    

################################## Masks ####################################

# Bounceback nodes mask (Loop = False, Branch = True)
bounceback = generateBouncebackMask(Lattice, False)

# Open path mask
openPath = invert(bounceback)

# Clot mask for clot in the upper tube
clotMask = generateClotMask(Lattice, Clot)

# Force array resistance mask for porous region
K = generateK(Lattice, Clot, clotMask)

# Pulsing field for fluid aceleration in the lower left tube section
pulseField = generatePulseFieldMask(Lattice)

# Accelerating force values
F = zeros((2,Lattice.nx, Lattice.ny))
F[0,pulseField] = Fluid.F_initial[0]
F[1,pulseField] = Fluid.F_initial[1]

##################### Initialising Output Monitoring Functions #####################
class Directory:
    clotVel = False
    vel = False
    clot = True
    tpaRho = False
    tpaFlow = True
    tpaConcentration = True
    clotLeft = True

# Generating working directories
mainDir, _, _, clotDir, _, tPAFlowDir, tPAConcDir, _ = createRepositories(Lattice, Fluid, Clot, TPA, Directory)

KLeftMost = []

############################# System Initliaization #################################

# Velocity initialisation
vel = zeros((2,Lattice.nx, Lattice.ny))

# Density initialisation
rho = full((Lattice.nx, Lattice.ny), Fluid.rho_initial)

# Initialisation of the populations at equilibrium with the given density & velocity.
fin = equilibrium(rho, vel, Lattice, D2Q9)
fout = equilibrium(rho, vel, Lattice, D2Q9)

# Loading already converged variables for faster execution time
fin, fout, _, u = getVariables("loop_Kxy", Lattice, Fluid, Clot, 65000)

# tPA density initialisation
rhoTPA = zeros((Lattice.nx, Lattice.ny))
rhoTPA[1:Lattice.tubeSize+1, Lattice.ny//2] = TPA.rho_initial

# tPA population initialisation
tPAin = equilibriumTPA(rhoTPA, u, Lattice, D2Q4)
tPAout = equilibriumTPA(rhoTPA, u, Lattice, D2Q4)

################################# Main time loop ######################################

# Monitoring execution time
start_time = time.time()

# main loop
for execTime in range(Lattice.maxIter):

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin, Lattice, D2Q9)       # fluid 
    rhoTPA = macroscopicTPA(tPAin)  # tPA 

    # injecting tPA constantly
    rhoTPA[1:Lattice.tubeSize+1, Lattice.ny//2] = Fluid.rho_initial

    # Compute equilibrium.
    feq = equilibrium(rho, u, Lattice, D2Q9)           # fluid
    tPAeq = equilibriumTPA(rhoTPA, u, Lattice, D2Q4)   # tPA

    # BGK collision step for open path
    fout[:,openPath] = fin[:,openPath] - Fluid.omega * (fin[:,openPath] - feq[:,openPath])            # fluid
    tPAout[:,openPath] = tPAin[:,openPath] - Fluid.omega * (tPAin[:,openPath] - tPAeq[:,openPath])    # tPA

    # Bounce-back condition 
    for i in range(9):                                  # fluid
        fout[i, bounceback] = fin[8-i, bounceback]
    for i in range(4):                                  # tPA
        tPAout[i, bounceback] = tPAin[3-i, bounceback]
    
    # Force field injection
    fout += addResistingClotForce(rho, u, F, K, Lattice, D2Q9)
    
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
    tPAin[0,:,:] = roll(tPAout[0,:,:],1,axis=0)                 # i = 0
    tPAin[1,:,:] = roll(tPAout[1,:,:],1,axis=1)                 # i = 1
    tPAin[2,:,:] = roll(tPAout[2,:,:],-1,axis=1)                # i = 2
    tPAin[3,:,:] = roll(tPAout[3,:,:],-1,axis=0)                # i = 3

    # Dissolve clot
    tPAin, K = bindAndDissolve(tPAin, K, Clot)

    if (execTime%250==0 and 13000<execTime<45000):
        # saveTPADensity(tPARhoDir, rhoTPA, execTime)
        showClotForce(clotDir, Clot, K, clotMask, execTime)
        saveTPAFlow(tPAFlowDir, rhoTPA, u, execTime)
    
    if (execTime%1000==0):
        plotTPAConcentration(tPAConcDir, Lattice, Clot, tPAin, K, bounceback, execTime)
    
    if(execTime%20==0 and 13000<execTime<45000):
        KLeftMost.append(leftMostK(K, Clot, clotMask))

    # Displaying current progress
    print("iteration : " + str(execTime) + "/" + str(Lattice.maxIter), end="\r")
    

# Final execution time
end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

######################## Final Iteration Monitoring ########################## 

# Final graphics generation
plotKLeftMost(mainDir, KLeftMost, Lattice, Clot)
plotTPAConcentration(tPAConcDir, Lattice, Clot, tPAin, K, bounceback, execTime+1)
# plotResults(clotVelDir, Lattice, Clot, u, rho, Lattice.maxIter)
# plotVelocityProfiles(velDir, Lattice, u, Lattice.maxIter)
# plotVelocityProfiles(velDir, Lattice, u, Lattice.maxIter)

########################### Converged System Saving ############################# 

# Saving converged system to load directly at next run
# saveVariables("loop_Kxy", Lattice, Fluid, Clot, fin, fout, rho, u)