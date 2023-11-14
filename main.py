import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from infinitesquarewell import InfiniteSquareWell
from potentials import PotentialType
from generatehamiltonian import compute_hamiltonian

def solve_problem(potential_choice=PotentialType.square, potential_amplitude=0.0, e_vals=5, l_bnd=-1.0, r_bnd=1.0):
    ISW = InfiniteSquareWell(energy_eigenvals=e_vals, well_min=l_bnd, well_max=r_bnd)
    potential = PotentialType.square
    V = potential.get_potential(ISW, potential_amplitude)
    return ISW, V
    
def main():
    # plot 2d wavefunction
    ISW, V = solve_problem()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    nx = 1
    ny = 2
    print(type(V))
    #ax.plot_surface(ISW.xvals, ISW.yvals, V,
    #                cmap='viridis')
    #plt.title(f'2D Potential')
    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.show()

if __name__ == '__main__':
    main()
