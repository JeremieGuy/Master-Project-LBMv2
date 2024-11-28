from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.colors import ListedColormap
import os
import pickle


# initialising the directories to save the output
def createRepositories(lattice, fluid, clot, tpa, Dir):

    # Root monitoring directory
    root = "./Monitoring"
    if not os.path.exists(root):
        os.mkdir(root)
        print("Made new root monitoring directory : " + root)

    # Main working directory for current execution
    mainDir = root + "/FF"
    if lattice.branch:
        txt = "_branch"
    else :
        txt = "_loop"
    mainDir += txt
    mainDir += "_PBB_"
    mainDir += "_viscosity=" + str(fluid.viscosity) + "_Rho=" + str(fluid.rho_initial) 
    mainDir += "_rhoTPA=" + str(tpa.rho_initial)
    mainDir += "_d=" + str(clot.d)
    mainDir += "_g=" + str(clot.gamma)
    mainDir += "_F=" + str(fluid.F_initial) + "_K=" + str(clot.K_initial)
    mainDir += "_it=" + str(lattice.maxIter)
    if not os.path.exists(mainDir):
        os.mkdir(mainDir)
        print("Made new main monitoring directory : " + mainDir)

    # Directory initialisation
    clotVelDir = ""
    velDir = ""
    clotDir = ""
    tpaRhoDir = ""
    tpaFlowDir = ""
    tpaConcentrationDir = ""
    clotLeftDir = ""

    # Directory for clot velocity visualisation
    if Dir.clotVel:
        clotVelDir = mainDir + "/Clot_Velocity"
        if not os.path.exists(clotVelDir):
            os.mkdir(clotVelDir)
            print("Made new clot velocity directory : " + clotVelDir)

    # Directory for velocity profiles visualisation
    if Dir.vel:
        velDir = mainDir + "/Velocity_profiles"
        if not os.path.exists(velDir):
            os.mkdir(velDir)
            print("Made new velocity profile directory : " + velDir)

    # Directory for clot force rates
    if Dir.clot:
        clotDir = mainDir + "/Clot_force"
        if not os.path.exists(clotDir):
            os.mkdir(clotDir)
            print("Made new clot force directory : ", clotDir)
    
    # Directory for tPA density
    if Dir.tpaRho:
        tpaRhoDir = mainDir + "/tPA_Densities"
        if not os.path.exists(tpaRhoDir):
            os.mkdir(tpaRhoDir)
            print("Made new clot force directory : ", tpaRhoDir)

    # Directory for tPA flow
    if Dir.tpaFlow:
        tpaFlowDir = mainDir + "/tPA_flow"
        if not os.path.exists(tpaFlowDir):
            os.mkdir(tpaFlowDir)
            print("Made new clot force directory : ", tpaFlowDir)
    
    # Directory for tPA Concentration
    if Dir.tpaConcentration:
        tpaConcentrationDir = mainDir + "/tPA_concentration"
        if not os.path.exists(tpaConcentrationDir):
            os.mkdir(tpaConcentrationDir)
            print("Made new clot force directory : ", tpaConcentrationDir)

    # Directory for Clot leftmost value
    if Dir.clotLeft:
        clotLeftDir = mainDir + "/clot_leftmost"
        if not os.path.exists(clotLeftDir):
            os.mkdir(clotLeftDir)
            print("Made new clot force directory : ", clotLeftDir)

    return mainDir, clotVelDir, velDir, clotDir, tpaRhoDir, tpaFlowDir, tpaConcentrationDir, clotLeftDir

# Drawing a figure of the system, with velocitiy profiles sites
def plotSystem(main_directory, lattice, bounceback, openPath, clot, pulseField):
    # Defining the plotting image 
    nx = lattice.nx
    ny = lattice.ny
    tubeSize = lattice.tubeSize

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
def plotResults(clot_directory, lattice, clot, u, rho, bounceback, execTime):
    nx = lattice.nx
    tubeSize = lattice.tubeSize
    clotCoord = clot.coord

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

    cmap2 = ListedColormap([(1, 1, 1, 0),  # Fully transparent for 0 (RGBA)
                       (0, 0, 0, 1)]) # Black for 1 (RGBA)
    borders = zeros((lattice.nx,lattice.ny))
    borders[bounceback] = 1

    ax0.imshow(borders.transpose(), cmap=cmap2)

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
def plotVelocityProfiles(velocity_directory, lattice, u, execTime):
    nx = lattice.nx
    ny = lattice.ny
    tubeSize = lattice.tubeSize

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
    

    # print(f"Does the path exist? {os.path.exists(velocity_directory)}")
    # Saving 
    name = velocity_directory
    name += "/vel_prof_" + str(execTime)+ ".png"
    plt.savefig(name , bbox_inches='tight')

    # Cleanup
    # plt.show()
    plt.close()
    plt.clf()

# Visualising the current system velocity norms
def visualise(u):
    plt.clf()
    plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    plt.pause(.01)
    plt.cla()

# Visualising tPA density
def visualiseTPA(rho):
    plt.clf()
    plt.imshow(rho.transpose(), cmap=cm.Reds)
    plt.pause(.01)
    # plt.show()
    plt.cla()

# Saving the visual representation of tPA density
def saveTPADensity(directory, rho, it):
    plt.clf()
    plt.imshow(rho.transpose(), cmap=cm.Reds)
    plt.title("iteration : " + str(it))
    name = directory + "/tPA_density_it=" + str(it)
    plt.savefig(name, bbox_inches='tight')
    plt.close()
    plt.clf()

# Saving the visual representation of tPA flow
def saveTPAFlow(directory, rho, u, it):
    plt.clf()
    plt.imshow((rho*sqrt(u[0]**2+u[1]**2)).transpose(), cmap=cm.Reds)
    plt.title("iteration : " + str(it))
    name = directory + "/tPA_flow_it=" + str(it)
    plt.savefig(name, bbox_inches='tight')
    plt.close()
    plt.clf()

# Visualisation of the resisting forces inside the clot
def showClotForce(directory, clot, K, clotMask, execTime):
    # Blank state
    plt.clf()

    # Showing the porous media resistance value
    # Note : both x & y decrease simutanously so displaying only one of them is sufficient
    rows, cols = where(clotMask)
    row_start, row_end = rows.min(), rows.max() + 1
    col_start, col_end = cols.min(), cols.max() + 1

    result = K[0, row_start:row_end, col_start:col_end]

    minVal = 0
    maxVal = clot.K_initial[0]

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
def saveVariables(type, lattice, fluid, clot, fin, fout, rho, u):

    # Creating variable storing directory
    varFolder = "./Variables"
    if not os.path.exists(varFolder):
        os.mkdir(varFolder)
        print("Made new variables storing directory : " + varFolder)

    # File for dumping objects containing all the fluid variables for reference
    filename = varFolder + "/" + type + "_"+str(lattice.nx)+"x"+str(lattice.ny)+"_viscosity="
    filename += str(fluid.viscosity) + "_Rho=" + str(fluid.rho_initial) 
    filename += "_F=" + str(fluid.F_initial) + "_K=" + str(clot.K_initial)
    filename += "_it=" + str(lattice.maxIter)

    # Generating the file
    with open(filename + ".pkl", 'wb') as f:
        pickle.dump([fin, fout, rho, u], f)
    
    # Closing the file
    f.close()

# Getting saved already converged variables to start the system
def getVariables(type, lattice, fluid, clot, it):


    # Get correct filename
    varFolder = "./Variables"
    filename = varFolder + "/" + type + "_"+str(lattice.nx)+"x"+str(lattice.ny)+"_viscosity=" 
    filename+= str(fluid.viscosity) + "_Rho=" + str(fluid.rho_initial) 
    filename += "_F=" + str(fluid.F_initial) + "_K=" + str(clot.K_initial)
    filename += "_it=" + str(it)

    # Recovering variables
    with open(filename + ".pkl", "rb") as f:  # Python 3: open(..., 'rb')
        fin, fout, rho, u = pickle.load(f)

    # Closing the file
    f.close()
    
    # Returning variables
    return fin, fout, rho, u
 
def saveDissolutionAmount(file, clot, dissolution):

    rows, cols = where(clot)
    row_start, row_end = rows.min(), rows.max() + 1
    col_start, col_end = cols.min(), cols.max() + 1


    file.write("\nDissolution Amount\n\n")
    
    for row in dissolution[row_start:row_end, col_start:col_end].transpose():
        file.write(" ".join(map(str, row.tolist())) + "\n")

def SaveTPAValues(file, before, tPAin, tPAout, K, clot, execTime):
    rows, cols = where(clot)
    row_start, row_end = rows.min(), rows.max() + 1
    col_start, col_end = cols.min(), cols.max() + 1

    line = "Coord : [" + str(row_start)
    line += ":" 
    line += str(row_end)
    line += "," + str(col_start)
    line += ":" + str(col_end) + "]\n"

    tmpin = sum(tPAin, axis=0)
    sumTPAin = tmpin[row_start:row_end, col_start:col_end]

    tmpout= sum(tPAout, axis=0)
    sumTPAout= tmpout[row_start:row_end, col_start:col_end]

    if before:
        file.write("Iteration : " + str(execTime) + "\n\n")
        file.write("######## Before Dissolution : ########\n")
        file.write(line)

        file.write("\ntPAin\n\n")

        for row in sumTPAin.transpose():
            file.write(" ".join(map(str, row.tolist())) + "\n")
                
        file.write("\ntPAout\n\n")
        
        for row in sumTPAout.transpose():
            file.write(" ".join(map(str, row.tolist())) + "\n")

        file.write("\nK\n\n")
                
        for row in K[0,row_start:row_end, col_start:col_end].transpose():
            file.write(" ".join(map(str, row.tolist())) + "\n")
    
    if not(before):
        file.write("\n######## After Dissolution ########\n")
        file.write("\ntPAin\n\n")

        for row in sumTPAin.transpose():
            file.write(" ".join(map(str, row.tolist())) + "\n")
        
        file.write("\ntPAout\n\n")
        
        for row in sumTPAout.transpose():
            file.write(" ".join(map(str, row.tolist())) + "\n")

        file.write("\nK\n\n")
                
        for row in K[0,row_start:row_end, col_start:col_end].transpose():
            file.write(" ".join(map(str, row.tolist())) + "\n")

        dotline = "\n"
        for i in range(121): dotline+="-"
        dotline += "\n\n"

        file.write(dotline)

# Determining how much the clot dissolve by 
# looking at the leftmost index where the average
# value of K is below 90% of initial
def leftMostK(K, Clot, clotMask):
    rows, cols = where(clotMask)
    row_start, row_end = rows.min(), rows.max() + 1
    col_start, col_end = cols.min(), cols.max() + 1

    Kmean = mean(K[0,row_start:row_end, col_start:col_end].transpose(), axis=0)
    threshold = 0.1 * Clot.K_initial[0]

    index = argmax(Kmean > threshold)

    return index

def plotKLeftMost(mainDir, kleftmost, lattice, clot):
    plt.clf()
    x = arange(0,len(kleftmost),1)
    plt.plot(x, kleftmost, "x")
    plt.title("Evolution of clot dissolution, it =" + str(lattice.maxIter))
    plt.xlabel("Time = dt*100 [s]")
    plt.ylabel("Left-most site index where K > 10%"+" of initial K")
    plt.ylim(0,clot.clotSize)

    name = mainDir + "/ClotEvolution_it=" + str(lattice.maxIter) + ".png"
    plt.savefig(name, bbox_inches="tight")


def plotTPAConcentration(directory, lattice, clot, tPAin, K, bounceback, execTime):
    nx = lattice.nx
    tubeSize = lattice.tubeSize
    clotCoord = clot.coord

    plt.clf()
    plotBiais = int(tubeSize*1.5)
    x_full = arange(1+plotBiais,clot.coord[0],1)

    
    tpaMean = mean(tPAin, axis=0).transpose()

    fig, (ax0, ax1, ax2,) = plt.subplots(3, 1, figsize=(8, 10))
    fig.suptitle("System state at iteration " + str(execTime), y=0.93)

    im0 = ax0.imshow(tpaMean, cmap=cm.Reds)
    ax0.set_title("Full system (Clot = [" + str(clotCoord[0]) + "," + str(clotCoord[1]) + "])")
    ax0.set_ylabel("Y coordinates")

    # Graphic display zone
    ax0.plot([1+plotBiais,clot.coord[0]+50],[1,1], color="red", linestyle="--", linewidth=1)
    ax0.plot([1+plotBiais,clot.coord[0]+50],[1+tubeSize,1+tubeSize], color="red", linestyle="--", linewidth=1)
    ax0.plot([1+plotBiais,1+plotBiais],[1,1+tubeSize], color="red", linestyle="--", linewidth=1)
    ax0.plot([clot.coord[0],clot.coord[0]],[1,1+tubeSize], color="blue", linestyle="--", linewidth=1)
    ax0.plot([clot.coord[0]+50,clot.coord[0]+50],[1,1+tubeSize], color="red", linestyle="--", linewidth=1)
    fig.colorbar(im0, ax=ax0, orientation='vertical', label="Velocity [m/s]", pad=0.02)

    # Clot
    cmap1 = ListedColormap([(1, 1, 1, 0),  # Fully transparent for zero values (RGBA)
                       (0, 0, 1, 1),  # Blue for non-zero values
                       (0, 1, 1, 1)]) # yellow for higher values

    ax0.imshow(K[0].transpose(), cmap=cmap1, alpha=1)

    # Bounceback borders
    cmap2 = ListedColormap([(1, 1, 1, 0),  # Fully transparent for 0 (RGBA)
                       (0, 0, 0, 1)]) # Black for 1 (RGBA)
    borders = zeros((lattice.nx,lattice.ny))
    borders[bounceback] = 1

    ax0.imshow(borders.transpose(), cmap=cmap2)

    tpaMean = mean(tPAin, axis=0)

    # tpa populations on x axis
    tPA_x = mean(tpaMean[(1+plotBiais):clot.coord[1]+50,1:(1+tubeSize+1)].transpose(), axis=0)
    # print(tpaMean.shape)
    # print(tpaMean[(1+plotBiais):clot.coord[0],1:(1+tubeSize+1)].shape)
    x_tpa1 = arange(1+plotBiais, 1+plotBiais+len(tPA_x), 1)

    ax1.plot(x_tpa1,tPA_x)
    ax1.set_ylabel("tPA mean density on axis = x")
    ax1.set_xlabel("x coordinates")
    ax1.set_ylim([0,.75])
    ax1.axvline(x=clotCoord[0], color='blue', linestyle='--', linewidth=1)
    ax1.axvline(x=clotCoord[1], color='blue', linestyle='--', linewidth=1)

    
    tPA_y = mean(tpaMean[(1+plotBiais):clot.coord[0],1:(1+tubeSize)].transpose(), axis=1)
    y_full = arange(0,len(tPA_y),1)
    ax2.plot(y_full, tPA_y)
    ax2.set_ylim([0,.75])
    ax2.set_ylabel("tPA mean density on axis = y")
    ax2.set_xlabel("y coordinates")

    # Saving
    name = directory + "/" + "tPA_population_" + str(execTime)
    plt.savefig(name, bbox_inches='tight')
    plt.close()

