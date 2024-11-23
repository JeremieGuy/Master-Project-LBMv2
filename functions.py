from numpy import *
# import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import os

# initialising the directories to save the output
def createRepositories(maxIter, nx, ny, viscosity, initialDensity, F_initial, K_initial):

    # Root monitoring directory
    root = "./Monitoring"
    if not os.path.exists(root):
        os.mkdir(root)
        print("Made new root monitoring directory : " + root)

    # Main working directory for current execution
    new_dir_monitoring = "./Monitoring/FF_branch_clean_"+str(nx)+"x"+str(ny)+"_viscosity=" + str(viscosity) + "_Rho=" + str(initialDensity) 
    new_dir_monitoring += "_F=" + str(F_initial) + "_K=" + str(K_initial)
    new_dir_monitoring += "_it=" + str(maxIter)
    if not os.path.exists(new_dir_monitoring):
        os.mkdir(new_dir_monitoring)
        print("Made new main monitoring directory : " + new_dir_monitoring)

    # Directory for clot velocity visualisation
    new_dir_clot_velocity = new_dir_monitoring + "/Clot_Velocity"
    if not os.path.exists(new_dir_clot_velocity):
        os.mkdir(new_dir_clot_velocity)
        print("Made new clot velocity directory : " + new_dir_clot_velocity)

    # Directory for velocity profiles visualisation
    new_dir_velocity = new_dir_monitoring + "/Velocity_profiles"
    if not os.path.exists(new_dir_velocity):
        os.mkdir(new_dir_velocity)
        print("Made new velocity profile directory : " + new_dir_velocity)

    return new_dir_monitoring, new_dir_clot_velocity, new_dir_velocity

# Drawing a figure of the system, with velocitiy profiles sites
def plotSystem(main_directory, nx, ny, tubeSize, bounceback, openPath, clot, pulseField):
    # Defining the plotting image 
    # open path = 0
    flags_plot = zeros((nx,ny))
    flags_plot[openPath] = 0

    # bounceback = 1
    flags_plot[bounceback] = 1

    # clot = 2
    flags_plot[clot] = 2

    # aceleration pulse field = 3
    flags_plot[pulseField] = 3

    # Generating a figure with labels
    norm = plt.Normalize(vmin=0,vmax=0.7)
    flagsname = ["Open Path", "Bounceback","Clot","Acceleration Field","Inlet", "Outlet"]
    plt.figure(figsize=(7.9,4))
    plt.title("Masks")
    values = unique(flags_plot.ravel())
    im = plt.imshow(flags_plot.transpose())
    colors = [ im.cmap(im.norm(value)) for value in values]
    patches = [ mpatches.Patch(color=colors[i], label=flagsname[i] ) for i in range(len(values)) ]
    plt.title("Flags")
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    
    # Adding velocity profiles lines
    rightProfile = 1 +3*tubeSize + ((ny-1-1-4*tubeSize)//2) # coordinates of the right tube

    # plt.plot([nx//2, nx//2], [1, 1+tubeSize], linestyle='--', color='red', linewidth=1) # top
    plt.plot([nx//2,nx//2],[1+2*tubeSize+1,1+3*tubeSize], linestyle='--', color='red', linewidth=1) # middle tube
    plt.plot([(nx-1-tubeSize), nx-2], [rightProfile, rightProfile], linestyle='--', color='red', linewidth=1) # right
    plt.plot([nx//2, nx//2], [(ny-1-tubeSize), ny-2], linestyle='--', color='red', linewidth=1) # bottom

    # Saving
    plt.savefig(main_directory + "/system.png",bbox_inches='tight')
    
    # Cleanup
    plt.show() 
    plt.close()


# Displaying the results
def plotResults(clot_directory, nx, tubeSize, clotCoord, u, rho, execTime):

    # Plotting variables
    plt.clf()
    plotBiais = int(tubeSize*1.5)
    x_full = arange(1+plotBiais,nx-1-plotBiais,1)

    # Plot
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, 1, figsize=(8, 10))
    fig.suptitle("System state at iteration " + str(execTime), y=0.93)

    # Full system visual
    im0 = ax0.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    ax0.set_title("Full system (Clot = [" + str(clotCoord[0]) + "," + str(clotCoord[1]) + "])")
    ax0.set_ylabel("Y coordinates")

    # Graphic display zone
    ax0.plot([1+plotBiais,nx-1-plotBiais],[1,1], color="red", linestyle="--", linewidth=1)
    ax0.plot([1+plotBiais,nx-1-plotBiais],[1+tubeSize,1+tubeSize], color="red", linestyle="--", linewidth=1)
    ax0.plot([1+plotBiais,1+plotBiais],[1,1+tubeSize], color="red", linestyle="--", linewidth=1)
    ax0.plot([nx-1-plotBiais,nx-1-plotBiais],[1,1+tubeSize], color="red", linestyle="--", linewidth=1)

    # Clot
    ax0.plot([clotCoord[0],clotCoord[1]],[1,1], color="blue", linestyle="--", linewidth=1)
    ax0.plot([clotCoord[0],clotCoord[1]],[1+tubeSize,1+tubeSize], color="blue", linestyle="--", linewidth=1)
    ax0.plot([clotCoord[0],clotCoord[0]],[1,1+tubeSize], color="blue", linestyle="--", linewidth=1)
    ax0.plot([clotCoord[1],clotCoord[1]],[1,1+tubeSize], color="blue", linestyle="--", linewidth=1)
    fig.colorbar(im0, ax=ax0, orientation='vertical', label="Velocity [m/s]", pad=0.02)

    # Velocities
    u_x = mean(u[0,(1+plotBiais):(nx-1-plotBiais),1:(1+tubeSize+1)].transpose(),axis=0)
    ax1.plot(x_full,u_x)
    ax1.set_ylabel("Velocity u(x) mean [m/s]")
    ax1.axvline(x=clotCoord[0], color='blue', linestyle='--', linewidth=1)
    ax1.axvline(x=clotCoord[1], color='blue', linestyle='--', linewidth=1)

    # Densities
    rho_x = mean(rho[1+plotBiais:nx-1-plotBiais,1:1+tubeSize+1],axis=1)
    ax2.plot(x_full,rho_x)
    ax2.set_ylabel("Density rho(x) mean")
    ax2.axvline(x=clotCoord[0], color='blue', linestyle='--', linewidth=1)
    ax2.axvline(x=clotCoord[1], color='blue', linestyle='--', linewidth=1)

    # Flow = rho*u
    ax3.plot(x_full,u_x*rho_x)
    ax3.set_ylabel("Flow u(x)*rho(x) [m/s]")
    ax3.set_xlabel("X Coordinates")
    ax3.axvline(x=clotCoord[0], color='blue', linestyle='--', linewidth=1)
    ax3.axvline(x=clotCoord[1], color='blue', linestyle='--', linewidth=1)

    # Saving
    name = clot_directory + "/" + "sanity_check_" + str(execTime)
    plt.savefig(name, bbox_inches='tight')
    plt.close()


# Plotting and saving the velocity profiles 
def plotVelocityProfiles(velocity_directory, nx, ny, tubeSize, u, execTime):

    # Top cylindier
    uTop = abs(u[0,nx//2,1+2*tubeSize+1:1+3*tubeSize+1])
    umaxTop = max(uTop)

    # Right cylinder
    rightProfile = 1 +3*tubeSize + ((ny-1-1-4*tubeSize)//2)
    uRight = abs(u[1,(nx-1-tubeSize):nx-1,rightProfile])
    umaxMid = max(uRight)

    # Bottom cylinder
    uBot = abs(u[0,nx//2,(ny-1-tubeSize):ny-1])
    umaxBot = max(uBot)

    # Velocity Plotting variables
    Rtube = tubeSize//2 # radius of tube sections
    r = abs(arange((-tubeSize//2)+1,(tubeSize//2)+1,1))
    x = arange(0,tubeSize,1)

    # Expected velocities (= parabola computed form u_max)
    expectedUTop = [umaxTop*(1-(i/Rtube)**2) for i in r]
    expecteduRight = [umaxMid*(1-(i/Rtube)**2) for i in r]
    expectedUBot = [umaxBot*(1-(i/Rtube)**2) for i in r]

    # Figure plot
    plt.clf()
    fig, (ax1,ax2, ax3) = plt.subplots(1, 3,figsize=(11, 4))
    fig.suptitle('Velocity Profiles')

    # Top tube
    ax1.plot(x,uTop, label = "Real Profile")
    ax1.plot(x,expectedUTop, label = "Expected Profile")
    ax1.set_title('Branch cylinder')
    ax1.set_xlabel("Tube width coordinates")
    ax1.set_ylabel("Velocity")
    ax1.legend()
    ax1.set_ylim([-0.01,umaxTop+0.01])

    # Right tube
    ax2.plot(x,uRight, label = "Real Profile")
    ax2.plot(x,expecteduRight, label = "Expected Profile")
    ax2.set_title('Right cylinder')
    ax2.set_xlabel("Tube width coordinates")
    ax2.legend()
    ax2.set_ylim([-0.01,umaxMid+0.01])

    # Bottom tube
    ax3.plot(x,uBot, label = "Real Profile")
    ax3.plot(x,expectedUBot, label = "Expected Profile")
    ax3.set_title('Bottom cylinder')
    ax3.set_xlabel("Tube width coordinates")
    ax3.legend()
    ax3.set_ylim([-0.01,umaxMid+0.01])

    # Plotting
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    
    # Saving 
    name = velocity_directory + "/" + "velocity_profiles_" + str(execTime)
    plt.savefig(name, bbox_inches='tight')

    # Cleanup
    # plt.show()
    plt.close()
    plt.clf()

# Visualing the current system velocity norms
def visualise(u):
        plt.clf()
        plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
        plt.pause(.01)
        plt.cla()
