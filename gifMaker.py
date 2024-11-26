import imageio
import numpy as np


directory = "Monitoring/FF_loop_D2Q4_tPA_dissolve_260x200_viscosity=0.01_Rho=2.5_rhoTPA_ini=1_d=1_F=[0, -0.0001]_K=[0.001, 0.001]_it=40000/Clot_force"
maxIter = 5000
plots = 10

index = np.arange(8250,36000,250)
print(index)

frames = []
# t = maxIter//plots

print("Making Gif ...")
for i in range(len(index)):
    # print(i)
    # num = "{0:0=5d}".format(i)
    image = imageio.v2.imread(f"./" + directory + "/clot_FF_it=_" + str(index[i]) + ".png")
    # image = imageio.v2.imread(f"/clot_FF_it=_" + str(index[i]) + ".png")

    frames.append(image)

imageio.mimsave("./" + directory + "/CLOT_DISSOLUTION.gif", frames, duration = 80)

# print("\nDone.")