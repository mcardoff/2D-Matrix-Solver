import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from infinitesquarewell import InfiniteSquareWell
from potentials import PotentialType
from generatehamiltonian import compute_hamiltonian

def solve_problem(text_obj, potential_choice, potential_amplitude,
                  e_vals=10, l_bnd=0.0, r_bnd=1.0):
    ISW = InfiniteSquareWell(energy_eigenvals=5, well_min=0.0, well_max=1.0)
    potential = potential_choice
    V = potential.get_potential(ISW, potential_amplitude)
    
def main():
    # plot 2d wavefunction
    ISW = InfiniteSquareWell(energy_eigenvals=5, well_min=0.0, well_max=1.0)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(ISW.xvals, ISW.yvals, ISW.basis_funcs[1][1],
                    cmap='viridis')
    plt.title(f'2D Wavefunction ($n_x={2}$, $n_y={2}$)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

if __name__ == '__main__':
    main()
