# libraries
import numpy as np
import numpy.linalg as la
import tkinter

# matplotlib backends
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

# file imports
from infinitesquarewell import InfiniteSquareWell
from potentials import PotentialType
from incdecbutton import IncDecButton
from generatehamiltonian import compute_hamiltonian

def solve_problem(potential_choice=PotentialType.square, potential_amplitude=0.0, e_vals=5,
                  x_min=-1.0, x_max=1.0, y_min=-1.0, y_max=1.0):
    # Define ISW basis using supplied values
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

def _quit(root):
    root.quit()  # stops mainloop
    root.destroy()

    
def main():
    global canvas, root

    # set up tkinter window
    root = tkinter.Tk()
    root.geometry("1000x600")
    root.wm_title("2-D Schrodinger")

    # Default choice is square well upon start
    potential_choice = PotentialType.square
    potential_amp = 0.0

    # add matplotlib hook to tk
    fig = Figure(figsize=(5, 4), dpi=100)
    subfig = fig.add_subplot(111, projection='3d')
    fig.suptitle(potential_choice.name)

    # solve the problem with initial choice
    x, y, V, funcs, vals = solve_problem(potential_choice=potential_choice, potential_amplitude=potential_amp)

    # connect matplotlib hook to tk root
    canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
    canvas.draw()

    # Quit button to exit
    quit_button = tkinter.Button(master=root, text="Quit", command=lambda: _quit(root))

     # Helper class that has button functions
    inc_dec = IncDecButton(subfig, canvas, x, y, funcs, V)
    inc_dec.init_plot()

    # prev eigenfunction
    prev_button = tkinter.Button(
        master=root, text="Prev Plot", command=lambda: inc_dec.dec_selector())

    # next eigenfunction
    next_button = tkinter.Button(
        master=root, text="Next Plot", command=lambda: inc_dec.inc_selector())


    canvas.get_tk_widget().pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=1)
    prev_button.pack(side=tkinter.TOP)
    next_button.pack(side=tkinter.TOP)
    quit_button.pack(side=tkinter.BOTTOM)

    tkinter.mainloop()

#       for (energy,func) in zip(vals,funcs):
#        fig = plt.figure()
#        ax = fig.add_subplot(111, projection='3d')
#        ax.plot_surface(x, y, func, cmap='viridis')
#        plt.title(f'2D Eigenfunc: {energy}')
#        plt.xlabel('x')
#        plt.ylabel('y')
#        plt.show()
#
#    plt.ioff()
if __name__ == '__main__':
    main()
