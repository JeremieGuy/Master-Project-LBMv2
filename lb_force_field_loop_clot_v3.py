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
# iterations
maxIter = 50000# Total number of time iterations.

# system Size
nx, ny = 260, 200 # Number of lattice nodes.
R = ny//2       # rayon of tube section
tubeSize = 21 # diameters of the tubes in the system
Rtube = tubeSize//2 #radius of smaller tube sections

# system vairbales
Re = 10.0         # Reynolds number.
# nulb    = uLB*R*2/Re;             # Viscoscity in lattice unitss. velocity*characteristic length (= H)/Raynolds
nulb = 0.01
omega = 1 / (3*nulb+0.5);    # Relaxation parameter.
velocity = 0.04             # Velocity in lattice units.
cs2 = 1/3                # sound veocity adapted to lattice units     
initialDensity = 2.5 # inital density of the system

# Force
F = [0,-0.0001] # pushing force F[x,y]
K_initial = [0.001,0] # initial resisting force of porous region K[2,x,y]

# Clot
clotSize = 20 # size of the clot lenghtwise in a tube section
clotCoord = [nx//2-clotSize//2,nx//2+clotSize//2] # clot coordinates


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

# setting relevant BB nodes to 0 on the system
def setBBNodeToZero():
    fin[:,bounceback] = 0
    fout[:,bounceback] = 0
    rho[bounceback] = 0

# Adding a force to u
def addResistingClotForce():

    # computing a resisting force on the porous region
    FF = zeros((9,nx,ny))
    for i in range(9):
        FF[i,:,:] = rho*(v[i,0]*(F[0] - K[0,:,:]*u[0,:,:]) + v[i,1]*F[1])
        FF[i,:,:] = FF[i,:,:] * (t[i] / cs2)

    # computing a resisting force on the porous region
    # FF = zeros((9,nx,ny))
    # for i in range(9):
    #     FF[i,clot] = rho[clot]*(v[i,0]*(F[0] - K[0,clot]*u[0,clot]) + v[i,1]*F[1])
    #     FF[i,clot] = FF[i,clot] * (t[i] / cs2)

    return FF


############################ Monitoring functions #################################

def plotSystem():
    norm = plt.Normalize(vmin=0,vmax=0.7)
    # flagsname = ["Open Path", "Bounceback", "Clot","Inlet", "Outlet"]
    flagsname = ["Open Path", "Bounceback","Clot","Acceleration Field","Inlet", "Outlet"]
    plt.figure(figsize=(7.9,4))
    plt.title("Flags")
    values = unique(flags_plot.ravel())
    
    im = plt.imshow(flags_plot.transpose())
    colors = [ im.cmap(im.norm(value)) for value in values]
    patches = [ mpatches.Patch(color=colors[i], label=flagsname[i] ) for i in range(len(values)) ]
    plt.title("Flags")
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    
    # velocity profiles lines
    # plt.plot([nx//2, nx//2], [1, 1+tubeSize], linestyle='--', color='red', linewidth=1) # top
    plt.plot([(nx-1-tubeSize), nx-2], [ny//2, ny//2], linestyle='--', color='red', linewidth=1) # right
    plt.plot([nx//2, nx//2], [(ny-1-tubeSize), ny-2], linestyle='--', color='red', linewidth=1) # bottom
    plt.savefig(new_dir_monitoring + "/system.png",bbox_inches='tight')
    # plt.show()
    
    plt.close()


def plotResults():

    plt.clf()
    plotBiais = int(tubeSize*1.5)
    x_full = arange(1+plotBiais,nx-1-plotBiais,1)

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, 1, figsize=(8, 10))

    fig.suptitle("System state at iteration " + str(execTime), y=0.93)

    # full system visual
    im0 = ax0.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    ax0.set_title("Full system (Clot = [" + str(clotCoord[0]) + "," + str(clotCoord[1]) + "])")
    ax0.set_ylabel("Y coordinates")

    # graphic display zone
    ax0.plot([1+plotBiais,nx-1-plotBiais],[1,1], color="red", linestyle="--", linewidth=1)
    ax0.plot([1+plotBiais,nx-1-plotBiais],[1+tubeSize,1+tubeSize], color="red", linestyle="--", linewidth=1)
    ax0.plot([1+plotBiais,1+plotBiais],[1,1+tubeSize], color="red", linestyle="--", linewidth=1)
    ax0.plot([nx-1-plotBiais,nx-1-plotBiais],[1,1+tubeSize], color="red", linestyle="--", linewidth=1)

    # clot
    ax0.plot([clotCoord[0],clotCoord[1]],[1,1], color="blue", linestyle="--", linewidth=1)
    ax0.plot([clotCoord[0],clotCoord[1]],[1+tubeSize,1+tubeSize], color="blue", linestyle="--", linewidth=1)
    ax0.plot([clotCoord[0],clotCoord[0]],[1,1+tubeSize], color="blue", linestyle="--", linewidth=1)
    ax0.plot([clotCoord[1],clotCoord[1]],[1,1+tubeSize], color="blue", linestyle="--", linewidth=1)
    fig.colorbar(im0, ax=ax0, orientation='vertical', label="Velocity [m/s]", pad=0.02)

    # velocities
    u_x = mean(u[0,(1+plotBiais):(nx-1-plotBiais),1:(1+tubeSize+1)].transpose(),axis=0)
    ax1.plot(x_full,u_x)
    ax1.set_ylabel("Velocity u(x) mean [m/s]")
    ax1.axvline(x=clotCoord[0], color='blue', linestyle='--', linewidth=1)
    ax1.axvline(x=clotCoord[1], color='blue', linestyle='--', linewidth=1)

    # densities
    rho_x = mean(rho[1+plotBiais:nx-1-plotBiais,1:1+tubeSize+1],axis=1)
    ax2.plot(x_full,rho_x)
    ax2.set_ylabel("Density rho(x) mean")
    ax2.axvline(x=clotCoord[0], color='blue', linestyle='--', linewidth=1)
    ax2.axvline(x=clotCoord[1], color='blue', linestyle='--', linewidth=1)

    # flow = rho*u
    ax3.plot(x_full,u_x*rho_x)
    ax3.set_ylabel("Flow u(x)*rho(x) [m/s]")
    ax3.set_xlabel("X Coordinates")
    ax3.axvline(x=clotCoord[0], color='blue', linestyle='--', linewidth=1)
    ax3.axvline(x=clotCoord[1], color='blue', linestyle='--', linewidth=1)

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
    y_range = [9,10,11]
    x_range = [58,59,60,61,62]

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
    y_range = [9,10,11]
    x_range = [58,59,60,61,62]
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

def plotVelocityProfiles():

    #### Final velocities
    # Top cylindier
    # uTop = abs(u[0,nx//2,1:tubeSize+1])
    # umaxTop = max(uTop)
    # Right cylinder
    uRight = abs(u[1,(nx-1-tubeSize):nx-1,ny//2])
    umaxMid = max(uRight)
    # Bottom cylinder
    uBot = abs(u[0,nx//2,(ny-1-tubeSize):ny-1])
    umaxBot = max(uBot)

    # velocity Plotting variables
    r = abs(arange((-tubeSize//2)+1,(tubeSize//2)+1,1))
    
    x = arange(0,tubeSize,1)

    # expected velocities
    # expectedUTop = [umaxTop*(1-(i/Rtube)**2) for i in r]
    expecteduRight = [umaxMid*(1-(i/Rtube)**2) for i in r]
    expectedUBot = [umaxBot*(1-(i/Rtube)**2) for i in r]

    # plot
    plt.clf()

    fig, (ax2, ax3) = plt.subplots(1, 2,figsize=(11, 4))

    fig.suptitle('Velocity Profiles')

    # ax1.plot(x,uTop, label = "Real Profile")
    # ax1.plot(x,expectedUTop, label = "Expected Profile")
    # ax1.set_title('Top cylinder')
    # ax1.set_xlabel("Tube width coordinates")
    # ax1.set_ylabel("Velocity")
    # ax1.legend()
    # ax1.set_ylim([-0.01,umaxTop+0.01])

    ax2.plot(x,uRight, label = "Real Profile")
    ax2.plot(x,expecteduRight, label = "Expected Profile")
    ax2.set_title('Right cylinder')
    ax2.set_xlabel("Tube width coordinates")
    # ax2.set_ylabel("Velocity")
    ax2.legend()
    ax2.set_ylim([-0.01,umaxMid+0.01])

    ax3.plot(x,uBot, label = "Real Profile")
    ax3.plot(x,expectedUBot, label = "Expected Profile")
    ax3.set_title('Bottom cylinder')
    ax3.set_xlabel("Tube width coordinates")
    # ax3.set_ylabel("Velocity")
    ax3.legend()
    ax3.set_ylim([-0.01,umaxMid+0.01])

    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    
    name = new_dir_velocity + "/" + "velocity_profiles_" + str(execTime)
    plt.savefig(name, bbox_inches='tight')
    # plt.show()
    plt.close()
    plt.clf()

##################### DEFINING CONTROL VARIABLES #####################

new_dir_monitoring = "./Monitoring/FF_loop_clotv3_FULLFF_"+str(nx)+"x"+str(ny)+"_viscosity=" + str(nulb) + "_Rho=" + str(initialDensity) 
new_dir_monitoring += "_F=" + str(F) + "_K=" + str(K_initial)
new_dir_monitoring += "_it=" + str(maxIter)

if not os.path.exists(new_dir_monitoring):
    os.mkdir(new_dir_monitoring)
    print("Made new monitoring directory : " + new_dir_monitoring)

new_dir_clot_velocity = new_dir_monitoring + "/Clot_Velocity"
if not os.path.exists(new_dir_clot_velocity):
    os.mkdir(new_dir_clot_velocity)
    print("Made new clot velocity directory : " + new_dir_clot_velocity)

new_dir_velocity = new_dir_monitoring + "/Velocity_profiles"
if not os.path.exists(new_dir_velocity):
    os.mkdir(new_dir_velocity)
    print("Made new velocity profile directory : " + new_dir_velocity)

populationFile = open(new_dir_monitoring + "/monitor_clot_nodes_it=" + str(maxIter) + ".txt", 'w')

dotline = "\n"
for i in range(121): dotline+="-"
dotline += "\n"

########################################### PLOTTING VARIABLES & FLAGS ############################

#  Bounceback nodes mask
bounceback = full((nx,ny),False)
bounceback[:,0] = True
bounceback[:,ny-1] = True
bounceback[0,:] = True
bounceback[nx-1,:] = True

# obstacle
# bounceback[22:128,22:78] = True
bounceback[tubeSize+1:(nx-1-tubeSize),tubeSize+1:(ny-1-tubeSize)] = True

# open path flags
openPath = invert(bounceback)

# Pulsing field for fluid aceleration
pulseField = full((nx,ny),False)
# pulseField[1:22,40:61] = True
pulseField[1:tubeSize+1,(ny//2-10):(ny//2+11)] = True


# clot
clot = full((nx,ny),False)
# clot[(nx-1-tubeSize):nx-1,(ny//2-10):(ny//2+11)] = True
clot[clotCoord[0]:clotCoord[1]+1,1:1+tubeSize+1] = True

# Force array resistance for porous region
K = zeros((2,nx,ny))
K[0,clot] = K_initial[0]
K[1,clot] = K_initial[1]


# K[0,bounceback] = 1
# adding a resisting force in x direction on the porous region
# K[0,clotCoord[0]:clotCoord[1]+1,1:ny-1] = K_initial[0]

#### draw system
# open path = 0
flags_plot = zeros((nx,ny))

flags_plot[openPath] = 0

# bounceback = 1
flags_plot[bounceback] = 1

# clot = 2
flags_plot[clot] = 2

# aceleration pulse field = 3
flags_plot[pulseField] = 3

# flags_plot[nx//2,ny//2] = 4

plotSystem()

# initializePopFile(populationFile)

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

setBBNodeToZero()
#################################### Main time loop #################################################
start_time = time.time()

for execTime in range(maxIter):

    # right wall: outflow condition.
    # fin[col3,-1,:] = fin[col3,-2,:]
    # drawoutput(populationFile)

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin)

    # Adding a pulsing horizontal force on the system
    # if (execTime%100<50):
    #     u[0,pulseField] += F[0]
    #     u[1,pulseField] += F[1]

    # Adding a constant horizontal force on the system
    u[0,pulseField] += F[0]
    u[1,pulseField] += F[1]

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

    # FF = addResistingClotForce()
    # fout[:,openPath] += FF[:,openPath] 
    

    # Bounce-back condition
    for i in range(9):
        fout[i, bounceback] = fin[8-i, bounceback]

    fout += addResistingClotForce()
    
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


    # if (execTime%10==0):
    #     plt.clf()
    #     # plot velocities
    #     plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    #     # plt.savefig("vel.{0:03d}.png".format(time//100))

    #     plt.pause(.01)
    #     plt.cla()
    
    if(execTime%1000==0):
        plotResults()
        plotVelocityProfiles()
    
    

    print("iteration : " + str(execTime) + "/" + str(maxIter), end="\r")

    # latticePopulation.append(sum(fin))


end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

#################### SYSTEM CHECKING ###################### 

# output result file
plotResults()
plotVelocityProfiles()

####################### COMMENTS & QUESTIONS #################################
# When changing dimensions, Change : flags, Mask, Bounceback nodes, inlet, outlet, velocity profiles