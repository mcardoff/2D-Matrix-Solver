import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from infinitesquarewell import InfiniteSquareWell
from potentials import PotentialType
from generatehamiltonian import compute_hamiltonian
import time

def solve_problem(potential_choice=PotentialType.square, potential_amplitude=0.0, e_vals=5,
                  x_min=-1.0, x_max=1.0, y_min=-1.0, y_max=1.0):
    # Define ISW basis using predefined min/max
    # TODO: Change this so you can have rectangular space, where x and y not necessary same min/max
    ISW = InfiniteSquareWell(energy_eigenvals=e_vals, well_x_min=x_min, well_x_max=x_max, well_y_min=y_min, well_y_max=y_max)

    # Extract potential meshgrid
    V = potential_choice.get_potential(ISW, potential_amplitude)

    # compute the hamiltonian, it is a (e_vals*e_vals)x(e_vals*e_vals) array
    H = compute_hamiltonian(V,ISW)  
    # diagonalize it
    eigenvals, eigenvecs = la.eig(H)
    
    # New eigenfunctions are colummn of eigenvec array times original corresponding eigenfunc
    # note that they are in the same order
    newfuncs = []
    for col in np.transpose(eigenvecs):
        lin_combination = np.zeros(V.shape)
        for (func, val) in zip(ISW.basis_funcs.values(), col):
            lin_combination += func*val
        newfuncs.append(lin_combination)
    
    x,y = ISW.xvals,ISW.yvals

    zipped = zip(eigenvals,newfuncs)
    sorted_zip = sorted(zipped)
    sorted_newfuncs = []
    for (_, func) in sorted_zip:
        sorted_newfuncs.append(func)

    # returns:
    # x :: valid x values for plotting
    # y :: valid y values for plotting
    # V :: Potential
    # newfuncs :: new eigenfuncs from diagonalized H
    # eigenvals :: new eigenvals of diagonalized H
    return (x,y,V,sorted_newfuncs,sorted(eigenvals))
    
def main():
    # plot 2d wavefunction
    x, y, V, funcs, vals = solve_problem()
    for (energy,func) in zip(vals,funcs):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x, y, func, cmap='viridis')
        plt.title(f'2D Eigenfunc: {energy}')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    plt.ioff()
if __name__ == '__main__':
    main()
