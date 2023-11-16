# libraries
import numpy as np
import numpy.linalg as la
import tkinter

# matplotlib backends
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

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

def create_validated_entry(master=None, width=5, validate=None, def_val=""):
    return_entry = tkinter.Entry(
        master, validate="key", validatecommand=(validate, '%P'), width=width)
    return_entry.insert(tkinter.END, def_val)
    return return_entry

def _on_item_select(list_box, button_obj, e_text_obj, amp_text_obj,
                    fig, x_min_obj, x_max_obj, y_min_obj, y_max_obj,
                    e_val_obj, event):
    """When an item in list_box is selected, recalculate the problem."""
    # global canvas, root

    potential = [enum for enum in PotentialType][list_box.curselection()[0]]
    potential_amp = float(amp_text_obj.get())
    x_min, x_max = float(x_min_obj.get()), float(x_max_obj.get())
    y_min, y_max = float(y_min_obj.get()), float(y_max_obj.get())
    e_vals = int(e_val_obj.get())

    # re-solve the problem
    x, y, V, funcs, vals = solve_problem(
        potential_choice=potential, potential_amplitude=potential_amp,
        e_vals=e_vals, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)

    # update figure title
    fig.suptitle(potential.to_string())

    # update button class
    button_obj.update_vals(x, y, funcs, V)

    # update energy text
    _format_energy_text(e_text_obj, vals)

def _validate_float(test):
    """Validate whether or not test is a valid floating number input."""
    # Ensure only one minus sign and decimal point
    replace_chars = ['.', '-']
    for char in replace_chars:
        test = test.replace(char, "", 1)
    return (test.isdigit() or test == "")


def _validate_int(test):
    """Validate whether or not test is a valid integer input."""
    # isdigit is siffucient
    return (test.isdigit() or test == "")

def _format_energy_text(text_obj, energy_vals):
    """Clear text currently in text_obj, replace with energy_vals."""
    text_obj.delete("1.0", "end")

    energy_string = ""
    for (i, val) in enumerate(np.sort(energy_vals)):
        num = "0{}".format(i+1) if i+1 < 10 else str(i+1)
        energy_string += f"E_{num} = {val:.2f}\n"

    text_obj.insert(tkinter.END, energy_string)

def _quit(root):
    root.quit()  # stops mainloop
    root.destroy()

    
def main():
    global canvas, root

    # set up tkinter window and needed sub windows
    root = tkinter.Tk()
    root.geometry("1000x600")
    root.wm_title("2-D Schrodinger")
    x_min_max_frame = tkinter.Frame(root)
    y_min_max_frame = tkinter.Frame(root)
    next_prev_frame = tkinter.Frame(root)

    # validaters for float and int
    reg_f = root.register(_validate_float)
    reg_i = root.register(_validate_int)

    # Default choice is square well upon start
    potential_choice = PotentialType.square
    potential_amp = 0.0

    # Text field containing potential amplitude
    amp_label_text = tkinter.StringVar(value="Potential Amp:")
    amp_label = tkinter.Label(root, textvariable=amp_label_text, height=2)

    amp_text = tkinter.Entry(
        root, validate="key", validatecommand=(reg_f, '%P'))
    amp_text.insert(tkinter.END, "0.0")

    # add matplotlib hook to tk
    fig = Figure(figsize=(4, 4), dpi=100)
    subfig = fig.add_subplot(111, projection='3d')
    fig.suptitle(potential_choice.name)
    
    # connect matplotlib hook to tk root
    canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=True)
    
    # change dimension of Hamiltonian
    eig_entry = create_validated_entry(master=root, validate=reg_i, def_val="5")

    # Quit button to exit
    quit_button = tkinter.Button(master=root, text="Quit", command=lambda: _quit(root))

    # solve the problem with initial choice
    x, y, V, funcs, vals = solve_problem(potential_choice=potential_choice, potential_amplitude=potential_amp)

    # Text containing energy eigenvalues:
    e_text = tkinter.Text(root, height=3, width=20)
    _format_energy_text(e_text, vals)    

    # Helper class that has button functions
    inc_dec = IncDecButton(subfig, canvas, x, y, funcs, V)
    inc_dec.init_plot()

    # prev eigenfunction
    prev_button = tkinter.Button(
        master=next_prev_frame, text="Prev Plot", command=lambda: inc_dec.dec_selector())

    # next eigenfunction
    next_button = tkinter.Button(
        master=next_prev_frame, text="Next Plot", command=lambda: inc_dec.inc_selector())

    # plot potential on top of original
    potential_button = tkinter.Button(
        master=root,
        text="Plot Potential",
        command=lambda: inc_dec.plot_potential())

    # change well minimum and maximum
    x_min_entry = create_validated_entry(master=x_min_max_frame, validate=reg_f, def_val=str(min(x[0])))
    x_max_entry = create_validated_entry(master=x_min_max_frame, validate=reg_f, def_val=str(min(x[-1])))
    y_min_entry = create_validated_entry(master=y_min_max_frame, validate=reg_f, def_val=str(min(y[0])))
    y_max_entry = create_validated_entry(master=y_min_max_frame, validate=reg_f, def_val=str(min(y[-1])))

    # listbox to pick potential
    potential_options = [potential.name for potential in PotentialType]
    list_items = tkinter.Variable(value=potential_options)
    listbox = tkinter.Listbox(
        root,
        listvariable=list_items,
        height=3,
        exportselection=False,
        selectmode=tkinter.BROWSE)
    listbox.selection_set(first=0)  # set first selection by default
    
    listbox.bind('<<ListboxSelect>>',
                 lambda x: _on_item_select(
                     listbox, inc_dec, e_text, amp_text,
                     fig, x_min_entry, y_max_entry, x_min_entry, y_max_entry,
                     eig_entry, x))
    
    # labels
    eng_label_text = tkinter.StringVar(value="Energy Values:")
    pot_label_text = tkinter.StringVar(value="2-D Potential:")
    x_bound_label_text = tkinter.StringVar(value="Well x Bounds:")
    y_bound_label_text = tkinter.StringVar(value="Well y Bounds:")
    eig_label_text = tkinter.StringVar(value="Energy Eigenvals:")
    energy_label = tkinter.Label(root, textvariable=eng_label_text, height=2)
    pot_label = tkinter.Label(root, textvariable=pot_label_text, height=2)
    x_bound_label = tkinter.Label(x_min_max_frame, textvariable=x_bound_label_text, height=2)
    y_bound_label = tkinter.Label(y_min_max_frame, textvariable=y_bound_label_text, height=2)
    eig_label = tkinter.Label(root, textvariable=eig_label_text, height=2)
    
    # pack buttons
    next_prev_frame.pack(side=tkinter.TOP, fill=None, expand=False)
    prev_button.pack(in_=next_prev_frame, side=tkinter.LEFT)
    next_button.pack(in_=next_prev_frame, side=tkinter.RIGHT)
    potential_button.pack(side=tkinter.TOP)
    pot_label.pack(side=tkinter.TOP)
    listbox.pack(side=tkinter.TOP)
    eig_label.pack(side=tkinter.TOP)
    eig_entry.pack(side=tkinter.TOP)
    amp_label.pack(side=tkinter.TOP)
    amp_text.pack(side=tkinter.TOP)
    x_min_max_frame.pack(side=tkinter.TOP, expand=False)
    y_min_max_frame.pack(side=tkinter.TOP, expand=False)
    x_bound_label.pack(in_=x_min_max_frame)
    x_min_entry.pack(in_=x_min_max_frame, side=tkinter.LEFT)
    x_max_entry.pack(in_=x_min_max_frame, side=tkinter.RIGHT)
    y_bound_label.pack(in_=y_min_max_frame)
    y_min_entry.pack(in_=y_min_max_frame, side=tkinter.LEFT)
    y_max_entry.pack(in_=y_min_max_frame, side=tkinter.RIGHT)
    energy_label.pack(side=tkinter.TOP)
    e_text.pack(side=tkinter.TOP)
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
