import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import random
from math import sqrt, pi


def creat_dataframe(name, low, high):
    """extract relevant data from files and store it in DataFrame"""

    m = []
    for r in range(low, high, 5):
        for theta in range(70, 161):
            f = open(path_to_folder + "/" + name + ".r" + "%.2f" %
                     (r/100) + "theta%s" % theta + ".0.out", "r")
            for line in f:
                if "SCF Done:" in line:
                    a = line.split()
                    m.append([float("%.2f" % (r/100)), theta, float(a[4])])
            f.close()

    df = pd.DataFrame(m, columns=["r", "theta", "Energy"])

    return df


def plot_PES(dataframe, name):
    """3D plot of the Potential Energy Surface, and store the plot as """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(molecule_df["r"], molecule_df["theta"], molecule_df["Energy"], cmap='plasma')
    ax.set_title("Potential Eenergy Surface for " + name)
    ax.set_xlabel('r/angstroms')
    ax.set_ylabel("theta/degrees")
    ax.set_zlabel("Energy/hartrees")
    plt.savefig("PES_" + name + ".png", dpi=100)
    plt.show()
    plt.close()

    return None


def get_equilibrium_info(dataframe):
    """find bond length, bond angle and energy at equilibrium state"""

    energy = dataframe["Energy"].min()
    r = dataframe["r"][dataframe.Energy == energy]  
    theta = dataframe["theta"][dataframe.Energy == energy]
    r = r.to_string(index=False)
    theta = theta.to_string(index=False)

    print("Energy at equilibrium geometry: %s Hartrees" % energy)
    print("H-X bond length at equilibrium geometry: %s angstroms" % r)
    print("X-H-X bond angle at equilibrium geometry: %s degrees" % theta)
    print("--------------Section-2-----------------")

    return energy, r, theta



def quadratic_fitting(dataframe):
    """fit the 10 data closest to minimum energy to the quardratic equation and get the coefficients"""

    r = []
    theta = []
    Energy_fixed_theta = []
    Energy_fixed_r = []

    eq_theta_subsection = dataframe.loc[dataframe["theta"] == float(theta_eq)].sort_values(by=["Energy"]).reset_index(drop=True) #truncate the section of dataframe with theta_eq, rank in ascending energy
    eq_r_subsection = dataframe.loc[dataframe["r"] == float(r_eq)].sort_values(by=["Energy"]).reset_index(drop=True) #truncate the section of dataframe with r_eq, rank in ascending energy
    

    for i in range(0, 10):
        r.append(float(eq_theta_subsection["r"][i]) * angstrom_to_m) #find r value of 10 lowest energy at theta_eq 
        theta.append(float(eq_r_subsection["theta"][i]) * degree_to_rad) #find theta value of 10 lowest energy at r_eq 
        Energy_fixed_theta.append((float(eq_theta_subsection["Energy"][i])) * hartree_to_J) #get corresponding energy
        Energy_fixed_r.append(float(eq_r_subsection["Energy"][i]) * hartree_to_J)

    constant_of_r = np.polyfit(r, Energy_fixed_theta, 2)
    constant_of_theta = np.polyfit(theta, Energy_fixed_r, 2)

    k_r = constant_of_r[0] * 2 #as 0.5*k_r=constant_of_r[0]
    k_theta = constant_of_theta[0]*2 
    print("The fitting of quardratic equation gives k_r = %d and k_theta = %s " % (k_r, k_theta))
    print("---------------Section-3---------------")

    return r, theta, Energy_fixed_theta, Energy_fixed_r, constant_of_r, constant_of_theta, k_r, k_theta



def plot_fitted_graph(xval, yval, coeff, fixed, xaxis, name):
    """plot to visualise the quality of quardratic fitting"""

    xfit = np.arange(min(xval)-1e-11, max(xval) + 1e-11, (max(xval)-min(xval))/1e3)
    yfit = np.polyval(coeff, xfit)
    plt.scatter(xval, yval, color="red", marker="x")
    plt.plot(xfit, yfit)
    plt.xlabel(xaxis)
    plt.ylabel("Energy/J")
    plt.title(name +" quadratic fitting curve for " + fixed)
    plt.savefig(name +" quadratic fitting curve for " + fixed + ".png", dpi=100)
    plt.show()
    plt.close()

    return None


def frequency_calculation (k1,k2):
    # unit is cm-1 hence need to divide by c
    freq_1 = sqrt(k1/(2*mu_to_kg))/(2*pi*c)
    freq_2 = sqrt(k2/(0.5*mu_to_kg*(float(r_eq) * angstrom_to_m)**2))/(2*pi*c)
    print("symmetric streching frequency:%s cm-1" % freq_1)
    print("bending frequency: %s cm-1" % freq_2)
    
    return None



# some constant for unit conversion
mu_to_kg = 1.66053886e-27
hartree_to_J = 4.35974417e-18
angstrom_to_m = 1e-10
degree_to_rad = pi/180
c = 2.99792458e10  # in cm


# get input from user
flag = False
while not flag:
    molecule = str(input("Enter H2O or H2S:"))
    if molecule in ["H2O", "H2S"]:
        flag = True
    else:
        print("Sorry, we only have information for H2O and H2S. Please try again")

path_to_folder = str(input("Enter the path to the folder containing the out files:"))

print(molecule + " has the following information:")
print("-------------Section-1--------------")

if molecule == "H2O":
    r_min, r_max = [70, 191]
elif molecule == "H2S":
    r_min, r_max = [60, 181]

#plot Potential Energy Surface and find the global minimum
molecule_df = creat_dataframe(molecule, r_min, r_max)
plot_PES(molecule_df, molecule)
energy_eq, r_eq, theta_eq = get_equilibrium_info(molecule_df)

#fitting, plotting and frequency calculation
fitting_data_r, fitting_data_theta, fitting_data_Energy_constant_theta, fitting_data_Energy_constant_r, coeff_r, coeff_theta, k_r, k_theta = quadratic_fitting(molecule_df) 
plot_fitted_graph(fitting_data_r, fitting_data_Energy_constant_theta, coeff_r, "fixed theta", "Bond length/angstrom", molecule)
plot_fitted_graph(fitting_data_theta, fitting_data_Energy_constant_r, coeff_theta, "fixed r", "Bond angle/degree", molecule)
frequency_calculation (k_r, k_theta)
print("---------------END-------------------")

