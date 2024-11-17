# Author : Guy Jérémie
# inspired by Jonas Latt

from numpy import *
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import os
import time

########################### Flow definition #############################################

maxIter = 8001  # Total number of time iterations.
Re = 10.0         # Reynolds number.
nx, ny = 20, 11 # Number of lattice nodes.
R = ny//2       # rayon of tube section
# nulb    = uLB*R*2/Re;             # Viscoscity in lattice units. velocity*characteristic length (= H)/Raynolds
nulb = 0.01
omega = 1 / (3*nulb+0.5);    # Relaxation parameter.
velocity = 0.04             # Velocity in lattice units.
cs2 = 1/3                # sound veocity adapted to lattice units     
initialDensity = 2.5 # inital density of the system
F = [0.0001,0] # pushing force F[x,y]
K_initial = [0.001,0] # initial resisting force of porous region K[2,x,y]
clotCoord = [8,12] # clot coordinates

########################## Lattice Constants ###########################################

v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
            [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ])
t = array([ 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])

col1 = array([0, 1, 2])
col2 = array([3, 4, 5])
col3 = array([6, 7, 8])

#################### Main Function Definitions ######################################

# macroscopic variable computation
def macroscopic(fin):
    rho = zeros((nx,ny))
    rho[openPath] = sum(fin[:,openPath], axis=0)
    
    u = zeros((2, nx, ny))
    for i in range(9):
        u[0,openPath] += v[i,0] * fin[i,openPath]
        u[1,openPath] += v[i,1] * fin[i,openPath]
    u[:,openPath] /= rho[openPath]
    return rho, u

# Equilibrium distribution function (rho = array)
def equilibrium(rho, u):          
    usqr = 3/2 * (u[0,openPath]**2 + u[1,openPath]**2)
    feq = zeros((9,nx,ny))
    for i in range(9):
        cu = 3 * (v[i,0]*u[0,openPath] + v[i,1]*u[1,openPath])
        feq[i,openPath] = rho[openPath]*t[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq


# setting relevant BB nodes to 0 on the 31x23 system
def setBBNodeToZero():
    # top border
    fin[:,:,ny-1] = 0
    # bottom border
    fin[:,:,0] = 0
    # right border
    # fin[:,nx-1,:] = 0

    # obstacle
    # fin[:,0:249,52:99] = 0

     # top border
    fout[:,:,ny-1] = 0
    # bottom border
    fout[:,:,0] = 0
    # right border
    # fout[:,nx-1,:] = 0

    # obstacle
    # fout[:,0:249,52:99] = 0

    rho[:,ny-1] = 0
    rho[:,0] = 0

# Adding a force to u
def addResistingClotForce():

    # computing a resisting force on the porous region
    FF = zeros((9,nx,ny))
    for i in range(9):
        FF[i,:,:] = rho*(v[i,0]*(F[0] - K[0,:,:]*u[0,:,:]) + v[i,1]*F[1])
        FF[i,:,:] = FF[i,:,:] * (t[i] / cs2)

    return FF


############################ Monitoring functions #################################

def plotSystem():
    norm = plt.Normalize(vmin=0,vmax=0.7)
    flagsname = ["open path", "bounceback", "clot","inlet", "outlet"]
    plt.figure(figsize=(7.9,4))
    plt.title("Flags")
    values = unique(flags_plot.ravel())
    
    im = plt.imshow(flags_plot.transpose())
    colors = [ im.cmap(im.norm(value)) for value in values]
    patches = [ mpatches.Patch(color=colors[i], label=flagsname[i] ) for i in range(len(values)) ]
    plt.title("Flags")
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    # plt.plot([10,10],[100,100], color="red", linewidth=3)
    plt.savefig(new_dir_monitoring + "/system.png",bbox_inches='tight')
    # plt.show()
    plt.close()

def plotResults():

    plt.clf()
    x_full = arange(0,nx,1)

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, 1, figsize=(8, 10))

    fig.suptitle("System state at iteration " + str(execTime), y=0.93)

    # full system visual
    im0 = ax0.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    ax0.set_title("Full system (Clot = [" + str(clotCoord[0]) + "," + str(clotCoord[1]) + "])")
    ax0.set_ylabel("Y coordinates")
    ax0.set_aspect('auto')
    ax0.axvline(x=clotCoord[0], color='red', linestyle='--', linewidth=1)
    ax0.axvline(x=clotCoord[1], color='red', linestyle='--', linewidth=1)

    fig.colorbar(im0, ax=ax0, orientation='vertical', label="Velocity [m/s]", pad=0.02)

    # velocities
    u_x = mean(u[0].transpose(),axis=0)
    ax1.plot(x_full,u_x)
    ax1.set_ylabel("Velocity u(x) mean [m/s]")
    ax1.axvline(x=clotCoord[0], color='red', linestyle='--', linewidth=1)
    ax1.axvline(x=clotCoord[1], color='red', linestyle='--', linewidth=1)

    # densities
    rho_x = mean(rho,axis=1)
    ax2.plot(x_full,rho_x)
    ax2.set_ylabel("Density rho(x) mean")
    ax2.axvline(x=clotCoord[0], color='red', linestyle='--', linewidth=1)
    ax2.axvline(x=clotCoord[1], color='red', linestyle='--', linewidth=1)

    # flow = rho*u
    ax3.plot(x_full,u_x*rho_x)
    ax3.set_ylabel("Flow u(x)*rho(x) [m/s]")
    ax3.set_xlabel("X Coordinates")
    ax3.axvline(x=clotCoord[0], color='red', linestyle='--', linewidth=1)
    ax3.axvline(x=clotCoord[1], color='red', linestyle='--', linewidth=1)
    # ax3.set_ylim([0, 0.03])

    name = new_dir_clot_velocity + "/" + "sanity_check_" + str(execTime)
    plt.savefig(name, bbox_inches='tight')

    plt.close()

def initializePopFile(populationFile):

    populationFile.write("MaxIter : " + str(maxIter) + ", system size : [" + str(nx) + "x"+ str(ny) + "]\n\n")
    # populationFile.write("System : \n\n" + str(transpose(flags_plot)) +"\n\n0 = open path, 1 = BB, 2 = inlet, 3 = outlet\n\n")

    ini = "Initial parameters :"
    ini += "\nMaxIter = " + str(maxIter) 
    ini += "\nsystem size = [" + str(nx) + "x"+ str(ny) + "]"
    ini += "\nViscosity = " + str(nulb)
    ini += "\nInitial Rho = " + str(initialDensity)
    ini += "\nForce F = " + str(F)
    ini += "\nRestiting porous region Force K = " + str(K)

    populationFile.write("Directions for respectively fin & fout : \n")
    y_range = [0,1,2,3,4,5,6,7,8,9,10]
    x_range = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]

    for y in y_range:
        # line1 = "y = " + str(y) + "\n"
        line0 = ""
        line1 = "" 
        line2 = ""
        line3 = ""

        for x in x_range:
            line0 += "(" + str(x) + "," + str(y) + ")     "
            line1 += "2" + " " + "5" + " " + "8" + "        "
            line2 += "1" + " " + "4" + " " + "7" + "        "
            line3 += "0" + " " + "3" + " " + "6" + "        "
        
        line0 += "|   "
        line1 += "|   "
        line2 += "|   "
        line3 += "|   "

        for x in x_range:
            line0 += "(" + str(x) + "," + str(y) + ")     "
            line1 += "6" + " " + "3" + " " + "0" + "        "
            line2 += "7" + " " + "4" + " " + "1" + "        "
            line3 += "8" + " " + "5" + " " + "2" + "        "
        
        populationFile.write("\n" + line0 + "\n" + line1 + "\n" + line2 + "\n" + line3 + "\n")

    populationFile.write("\nHorizontal Velocity                                                      | Vertical Velocity\n\n")

    for y in y_range:
        line = ""


        for x in x_range:
            line += "u[0," + str(x) + "," + str(y) + "]    "
        
        line += "|   "

        for x in x_range:
            line += "u[1," + str(x) + "," + str(y) + "]    "
        
        populationFile.write("\n" + line + "\n" )

    lineRho = "\nRho\n\n"

    for y in y_range:
        for x in x_range:
            lineRho += "rho[" + str(x) + "," + str(y) + "]   "
        
        lineRho += "\n"

    populationFile.write(lineRho)

    populationFile.write(dotline + "\n")

def drawoutput(populationFile):
    populationFile.write("Iteration : " + str(execTime) + "\n")
    # y_range = arange(0,ny,1)[::-1] 
    y_range = [0,1,2,3,4,5,6,7,8,9,10]
    x_range = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    # print(y_range)
    for y in y_range:
        line1 = ""
        line2 = ""
        line3 = ""

        for x in x_range:
            line1 += str('%.5f'%(fin[2,x,y])) + " " + str('%.5f'%(fin[5,x,y])) + " " + str('%.5f'%(fin[8,x,y]))  + "    "
            line2 += str('%.5f'%(fin[1,x,y])) + " " + str('%.5f'%(fin[4,x,y])) + " " + str('%.5f'%(fin[7,x,y])) + "    "
            line3 += str('%.5f'%(fin[0,x,y])) + " " + str('%.5f'%(fin[3,x,y])) + " " + str('%.5f'%(fin[6,x,y])) + "    "
        
        line1 += "|   "
        line2 += "|   "
        line3 += "|   "

        for x in x_range:
            line1 += str('%.5f'%(fout[6,x,y])) + " " + str('%.5f'%(fout[3,x,y])) + " " + str('%.5f'%(fout[0,x,y])) + "    "
            line2 += str('%.5f'%(fout[7,x,y])) + " " + str('%.5f'%(fout[4,x,y])) + " " + str('%.5f'%(fout[1,x,y])) + "    "
            line3 += str('%.5f'%(fout[8,x,y])) + " " + str('%.5f'%(fout[5,x,y])) + " " + str('%.5f'%(fout[2,x,y])) + "    "
        
        populationFile.write("\n" + line1 + "\n" + line2 + "\n" + line3 + "\n\n")
    
    populationFile.write("Horizontal velocity                                                           | Vertical velocity\n")

    for y in y_range:
        line = ""
        
        for x in x_range: line += str('%.5f'%(u[0,x,y])) + "    "

        line += "|   "

        for x in x_range: line += str('%.5f'%(u[1,x,y])) + "    "
        
        populationFile.write("\n" + line + "\n")

    lineRho = "\n Rho\n\n"
    for y in y_range:
        for x in x_range:
            lineRho += str(rho[x,y]) + "   "
        
        lineRho += "\n"

    populationFile.write(lineRho)

    populationFile.write(dotline + "\n")

##################### DEFINING CONTROL VARIABLES #####################


new_dir_monitoring = "./Monitoring/SmallSystem_FF_viscosity=" + str(nulb) + "_Rho=" + str(initialDensity) 
new_dir_monitoring += "_F=" + str(F) + "_K=" + str(K_initial)
new_dir_monitoring += "_it=" + str(maxIter)

if not os.path.exists(new_dir_monitoring):
    os.mkdir(new_dir_monitoring)
    print("Made new monitoring directory : " + new_dir_monitoring)

new_dir_clot_velocity = new_dir_monitoring + "/Clot_Velocity"
if not os.path.exists(new_dir_clot_velocity):
    os.mkdir(new_dir_clot_velocity)
    print("Made new clot velocity directory : " + new_dir_clot_velocity)

# new_dir_velocity = new_dir_monitoring + "/Velocity_profiles"
# if not os.path.exists(new_dir_velocity):
#     os.mkdir(new_dir_velocity)
#     print("Made new velocity profile directory : " + new_dir_velocity)

populationFile = open(new_dir_monitoring + "/monitor_clot_nodes_it=" + str(maxIter) + ".txt", 'w')

dotline = "\n"
for i in range(121): dotline+="-"
dotline += "\n"

########################################### PLOTTING VARIABLES & FLAGS ############################

#  Bounceback nodes mask
bounceback = full((nx,ny),False)
bounceback[:,0] = True
bounceback[:,ny-1] = True

# open path flags
openPath = invert(bounceback)


# Force array resistance for porous region
K = zeros((2,nx,ny))
K[0,bounceback] = 1
# adding a resisting force in x direction on the porous region
K[0,clotCoord[0]:clotCoord[1]+1,1:ny-1] = K_initial[0]

#### draw system
# open path = 0
flags_plot = zeros((nx,ny))

# bounceback = 1
flags_plot[:,0] = 1
flags_plot[:,ny-1] = 1

# clot = 2
flags_plot[clotCoord[0]:clotCoord[1]+1,1:ny-1] = 2

plotSystem()

initializePopFile(populationFile)

################################### System Initliaization ##########################################

# initial velocity
# vel = iniVel()
vel = zeros((2,nx,ny))

rho = full((nx,ny),initialDensity)
# rho[:,0] = 0
# rho[:,ny-1] = 0

# Initialization of the populations at equilibrium with the given density & velocity.
fin = equilibrium(rho, vel)
fout = equilibrium(rho, vel)

#### set BB nodes PDFs and density to 0 
# setBBNodeToZero()

# initial macrovariables for file output
rho2, u = macroscopic(fin)

# rho = full((nx,ny),initialDensity)
# rho[:,0] = 0
# rho[:,ny-1] = 0

# setBBNodeToZero()
#################################### Main time loop #################################################
start_time = time.time()

for execTime in range(maxIter):

    # right wall: outflow condition.
    # fin[col3,-1,:] = fin[col3,-2,:]
    # drawoutput(populationFile)

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Adding a horizontal force on the system
    u[0,openPath] += F[0]
    u[1,openPath] += F[1]

    # Left wall: inflow condition.
    # u[:,0,:] = vel[:,0,:]
    # rho[0,:] = 1/(1-u[0,0,:]) * ( sum(fin[col2,0,:], axis=0) + 2*sum(fin[col3,0,:], axis=0) )

    # Compute equilibrium.
    feq = equilibrium(rho, u)
    # fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]

    # Collision step for open path
    # fout = fin - omega * (fin - feq)
    fout[:,openPath] = fin[:,openPath] - omega * (fin[:,openPath] - feq[:,openPath])
    # fout[:,openPath] = fin[:,openPath]


    FF = addResistingClotForce()
    fout[:,openPath] += FF[:,openPath] 

    # Bounce-back condition
    for i in range(9):
        fout[i, bounceback] = fin[8-i, bounceback]


    
    # streaming step
    # i = 0
    fin[0,:,:] = roll(roll(fout[0,:,:],1,axis=0),1,axis=1)
    # i = 1
    fin[1,:,:] = roll(fout[1,:,:],1,axis=0)
    # i = 2
    fin[2,:,:] = roll(roll(fout[2,:,:],1,axis=0),-1,axis=1)
    # i = 3
    fin[3,:,:] = roll(fout[3,:,:],1,axis=1)
    # i = 4
    fin[4,:,:] = fout[4,:,:]
    # i = 5
    fin[5,:,:] = roll(fout[5,:,:],-1,axis=1)
    # i = 6
    fin[6,:,:] = roll(roll(fout[6,:,:],-1,axis=0),1,axis=1)
    # i = 7
    fin[7,:,:] = roll(fout[7,:,:],-1,axis=0)
    # i = 8
    fin[8,:,:] = roll(roll(fout[8,:,:],-1,axis=0),-1,axis=1)


    if (execTime%10==0):
        plt.clf()
        # plot velocities
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
        # plt.savefig("vel.{0:03d}.png".format(time//100))

        plt.pause(.01)
        plt.cla()
    
    if(execTime%1000==0):
        plotResults()
    
    

    print("iteration : " + str(execTime) + "/" + str(maxIter), end="\r")

    # latticePopulation.append(sum(fin))


end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

#################### SYSTEM CHECKING ###################### 

# output result file
plotResults()

####################### COMMENTS & QUESTIONS #################################
# When changing dimensions, Change : flags, Mask, Bounceback nodes, inlet, outlet, velocity profiles