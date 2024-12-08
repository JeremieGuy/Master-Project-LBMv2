import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

# Function to generate a random color
def random_color():
    return (random.random(), random.random(), random.random())

path  = "./Monitoring"
path += "/FF_branch=31_PBB__viscosity=0.01_Rho=2.5_rhoTPA=1_d=1_g=0.5_F=[0, -0.0001]_K=[0.001, 0.001]_it=110000"
path += "/clot_leftmost"

name =  "/KLeftMost_corrected.csv"

# Load the CSV file
data = pd.read_csv(path + name)

# Create equation output file
equations = open(path + "/equations.txt", 'w')

# Extract x and y
x = data['x'].values
y = data['y'].values

# Visualize the data
plt.scatter(x, y, label='K leftmost index', color='blue', alpha=0.3)

# Colors
colors = ["red", "green", "yellow", "purple"]

# Degrees
degrees = [1, 2, 3, 4]


for i in range(len(degrees)):
    d = degrees[i]
    # Fit a polynomial
    coefficients = np.polyfit(x, y, d)

    # Generate the polynomial function
    polynomial = np.poly1d(coefficients)

    # Generate fitted curve
    x_fit = np.linspace(min(x), max(x), 500)
    y_fit = polynomial(x_fit)



    # Plot the fitted curve
    plt.plot(x_fit, y_fit, label=f'Polynomial (degree {d})', color=colors[i])
    

    print(f"Degree {d} : ")
    print("Polynomial coefficients: ", coefficients)
    print("Polynomial equation: ", polynomial)

    equations.write(f"Degree {d} : \n")
    equations.write("Polynomial coefficients: " + str(coefficients) + "\n")
    equations.write("Polynomial equation: " + str(polynomial) + "\n\n")



plt.legend()
plt.title("K leftmost evolution fit")
plt.ylabel("Lefmost coordinate")
plt.xlabel("Iterations dt*100 [s]")
output = "/Polynomail_fit.png"
plt.savefig(path + output, bbox_inches='tight')
plt.show()

equations.close()
