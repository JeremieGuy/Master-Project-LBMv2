from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import os
import pickle

# Function call examples 
"""
# Creating monitoring directories
-> Has to be called first to generate directories to save results for the other monitoring functions
mainDirectory, clotDirectory, velocityDirectory, clotForceDirectory = createRepositories(maxIter,nx,ny,viscosity,rho_initial,F_initial,K_initial)
-> if only some directories are needed : 
mainDirectory, clotDirectory, _ , _= createRepositories(maxIter,nx,ny,viscosity,rho_initial,F_initial,K_initial)

# Displaying the system topology
plotSystem(mainDirectory, nx, ny, tubeSize, bounceback, openPath, clot, pulseField)

# Generating the velocity profiles at execTime iterations
-> if called before the main loop, set execTime to 0
plotVelocityProfiles(velocityDirectory, nx, ny, tubeSize, u, execTime)
    
# Generating the clot velocity, density and flow graphics at execTime iterations
-> if called before the main loop, set execTime to 0
plotResults(clotDirectory, ny, tubeSize, clotCoord, u, rho, execTime)

# Generating the image for the porous region dissolving rates
-> if called before the main loop, set execTime to 0
showClotForce(clotForceDirectory, K_initial, K, clot, execTime)

# Saving the converged system to reuse later
saveVariables(nx, ny, viscosity, rho_initial, F_initial, K_initial, maxIter, fin, fout, rho, u)

# Getting the saved already converged system 
fin, fout, _, _ = getVariables(nx, ny, viscosity, rho_initial, F_initial, K_initial, 60000)
-> the iterations are given directly to always take the same value we calculated
-> rho and u are not needed to run the simulation so we need to load only fin and fout
"""

# initialising the directories to save the output
def createRepositories(maxIter, nx, ny, viscosity, rho_initial, F_initial, K_initial):

    # Root monitoring directory
    root = "./Monitoring"
    if not os.path.exists(root):
        os.mkdir(root)
        print("Made new root monitoring directory : " + root)

    # Main working directory for current execution
    new_dir_monitoring = root + "/FULLMACRO_DISSOLVE_FF_branch_clean_"+str(nx)+"x"+str(ny)+"_viscosity=" + str(viscosity) + "_Rho=" + str(rho_initial) 
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

    # Directory for clot force rates
    new_dir_clot_force = new_dir_monitoring + "/Clot_force"
    if not os.path.exists(new_dir_clot_force):
        os.mkdir(new_dir_clot_force)
        print("Made new clot force directory : ", new_dir_clot_force)

    return new_dir_monitoring, new_dir_clot_velocity, new_dir_velocity, new_dir_clot_force

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

# Visualisation of the resisting forces inside the clot
def showClotForce(directory, K_initial, K, clot, execTime):
    # Blank state
    plt.clf()

    # Showing the porous media resistance value
    # Note : both x & y decrease simutanously so displaying only one of them is sufficient
    rows, cols = where(clot)
    row_start, row_end = rows.min(), rows.max() + 1
    col_start, col_end = cols.min(), cols.max() + 1

    result = K[0, row_start:row_end, col_start:col_end]

    minVal = 0
    maxVal = K_initial[0]

    plt.imshow(result.transpose(), vmin=minVal, vmax=maxVal)
    plt.colorbar(label='K Value')
    plt.title("Porous Site Resistance, it = " + str(execTime))
    name = directory + "/clot_FF_it=" + str(execTime)

    # Saving
    plt.savefig(name, bbox_inches="tight")

    # Cleanup
    # plt.show()
    plt.close()

# Saving variables to run simulations with an already converged system
def saveVariables(nx , ny, viscosity, rho_initial, F_initial, K_initial, maxIter, fin, fout, rho, u):
    # Creating variable storing directory
    varFolder = "./Variables"
    if not os.path.exists(varFolder):
        os.mkdir(varFolder)
        print("Made new variables storing directory : " + varFolder)

    # File for dumping objects containing all the fluid variables for reference
    filename = varFolder + "/branch_"+str(nx)+"x"+str(ny)+"_viscosity=" + str(viscosity) + "_Rho=" + str(rho_initial) 
    filename += "_F=" + str(F_initial) + "_K=" + str(K_initial)
    filename += "_it=" + str(maxIter)

    # Generating the file
    with open(filename + ".pkl", 'wb') as f:
        pickle.dump([fin, fout, rho, u], f)

# Getting saved already converged variables to start the system
def getVariables(nx, ny, viscosity, rho_initial, F_initial, K_initial, maxIter):
    # Get correct filename
    varFolder = "./Variables"
    filename = varFolder + "/branch_"+str(nx)+"x"+str(ny)+"_viscosity=" + str(viscosity) + "_Rho=" + str(rho_initial) 
    filename += "_F=" + str(F_initial) + "_K=" + str(K_initial)
    filename += "_it=" + str(maxIter)

    # Recovering variables
    with open(filename + ".pkl", "rb") as f:  # Python 3: open(..., 'rb')
        fin, fout, rho, u = pickle.load(f)
    
    # Returning variables
    return fin, fout, rho, u


