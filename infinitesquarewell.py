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

    def __init__(self, well_min=0.0, well_max=1.0, steps=200,
                 energy_eigenvals=5, hbar=1.0, mass=1.0):
        """Initialize given width, mass, number of evals and resolution."""
        # values set by user
        self.well_min = well_min
        self.well_max = well_max
        self.well_width = abs(well_max - well_min)
        self.steps = steps
        self.energy_eigenvals = energy_eigenvals
        self.step_size = self.well_width / steps

        # used in generation
        self.basis_funcs = {}
        self.eigenvals = {}
        self.xvals = []
        self.yvals = []

        # Set to 1 because we can
        self.hbar = hbar
        self.mass = mass
        
        # set x and y values
        self.xvals = np.linspace(self.well_min, self.well_max, self.steps+1)
        self.yvals = np.linspace(self.well_min, self.well_max, self.steps+1)
        self.xvals, self.yvals = np.meshgrid(self.xvals, self.yvals)

        # do everything in initializer
        self.generate_basis_funcs()

    def basis_2D(self, x, y, n_x, n_y):
        PI = np.pi
        L = self.well_width
        psi_x = np.sqrt(2/L)*np.sin(n_x*PI*(x-self.well_min)/L)
        psi_y = np.sqrt(2/L)*np.sin(n_y*PI*(y-self.well_min)/L)
        return psi_x * psi_y
    
    def generate_basis_funcs(self):
        """Generate eigenfunctions of zero potential well."""
        # know how to generate the infinite square well basis,
        # can base everything off that
        # Quick pneumonics
        PI = np.pi
        L = self.well_width
        # ISW eigenvalues are natural numbers
        for n in range(1, self.energy_eigenvals+1):
            for m in range(1, self.energy_eigenvals+1):
                # analytic formulae
                energy = (n*n + m*m) * (self.hbar * PI / L) ** 2 / (2.0*self.mass)
                eigenfunc = self.basis_2D(self.xvals, self.yvals, n, m)
                # Add to dicts
                self.eigenvals[(n,m)] = energy
                self.basis_funcs[(n,m)] = eigenfunc
            


if __name__ == "__main__":
    main()
