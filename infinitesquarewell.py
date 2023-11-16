"""Define an infinite square well basis object for the eigenfunctions."""
import numpy as np
import matplotlib.pyplot as plt


def main():
    """Test if the correct basis functions are generated."""
    # tests if you can generate a ISW object, plots basis
    ISW = InfiniteSquareWell()
    for func in ISW.basis_funcs:
        plt.plot(ISW.xvals, func)
    plt.show()


class InfiniteSquareWell:
    """Class which contains all important information about the ISW."""

    def __init__(self, well_x_min=0.0, well_x_max=1.0, well_y_min=0.0, well_y_max=1.0, steps=200,
                 energy_eigenvals=5, hbar=1.0, mass=1.0):
        """Initialize given width, mass, number of evals and resolution."""
        # values set by user
        self.well_x_min, self.well_x_max = min(well_x_min, well_x_max), max(well_x_min, well_x_max) # x bounds of well
        self.well_y_min, self.well_y_max = min(well_y_min, well_y_max), max(well_y_min, well_y_max) # y bounds of well
        self.well_x_width, self.well_y_width = abs(well_x_max - well_x_min), abs(well_y_max - well_y_min)
        self.well_area = self.well_x_width * self.well_y_width # for convenience when calculating matrix elements
        self.energy_eigenvals = energy_eigenvals # number of energy eigenvals to use in basis
        self.steps = steps

        # used in generation
        self.basis_funcs = {}
        self.eigenvals = {}

        # Set to 1 because we can
        self.hbar = hbar
        self.mass = mass
        
        # set x and y values
        temp_xvals = np.linspace(self.well_x_min, self.well_x_max, self.steps+1)
        temp_yvals = np.linspace(self.well_y_min, self.well_y_max, self.steps+1)
        self.xvals, self.yvals = np.meshgrid(temp_xvals, temp_yvals)
        del temp_xvals
        del temp_yvals
        
        # for convenience, coordinate pairs of the above:
        self.coord_pairs = np.column_stack((self.xvals.ravel(),self.yvals.ravel())) 

        # do everything in initializer
        self.generate_basis_funcs()

    def basis_2D(self, x, y, n_x, n_y):
        PI = np.pi
        L1 = self.well_x_width
        L2 = self.well_y_width
        psi_x = np.sqrt(2/L1)*np.sin(n_x*PI*(x-self.well_x_min)/L1)
        psi_y = np.sqrt(2/L2)*np.sin(n_y*PI*(y-self.well_x_min)/L2)
        return psi_x * psi_y
    
    def generate_basis_funcs(self):
        """Generate eigenfunctions of zero potential well."""
        # Quick pneumonics
        PI = np.pi
        L1 = self.well_x_width
        L2 = self.well_y_width
        # ISW eigenvalues are natural numbers
        for n in range(1, self.energy_eigenvals+1):
            for m in range(1, self.energy_eigenvals+1):
                # analytic formulae
                energy = ((n/L1)**2 + (m/L2)**2) * (self.hbar * PI)**2 / (2.0*self.mass)
                eigenfunc = self.basis_2D(self.xvals, self.yvals, n, m)
                # Add to dicts
                self.eigenvals[(n,m)] = energy
                self.basis_funcs[(n,m)] = eigenfunc
            


if __name__ == "__main__":
    main()
