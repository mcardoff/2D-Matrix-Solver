"""Class to organize and work with various potentials easily."""
import matplotlib.pyplot as plt
import numpy as np
from enum import Enum, auto
from infinitesquarewell import InfiniteSquareWell


def main():
    """Test If potential plot is correct."""
    ISW = InfiniteSquareWell()
    potential = PotentialType.linear
    V = potential.get_potential(ISW,50.0)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(ISW.xvals, ISW.yvals, V, cmap='viridis')
    plt.title(f'2D Potential: {potential.to_string()}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()


def general_well(ISW, f):
    """Place an infinite barrier at bounds, evals func provided elsewhere."""
    MXVAL = 10000
    ret = np.array([f(x,y) for [x,y] in ISW.coord_pairs]).reshape(ISW.xvals.shape)
    ret[0,:] = ret[:,0] = ret[-1,:] = ret[:,-1] = MXVAL # change boundaries to be MXVAL
    return ret

def square(ISW, amplitude):
    return general_well(ISW, lambda x, y: x*0 + amplitude)

def linear(ISW, amplitude):
    """Particle in an electric field along the x axis."""
    assert(isinstance(ISW, InfiniteSquareWell))
    return general_well(ISW, lambda x, y: amplitude * x)


def quadratic(ISW, amplitude):
    """Half-harmonic oscillator potential."""
    assert(isinstance(ISW, InfiniteSquareWell))
    return general_well(ISW, lambda x, y: amplitude * ((x*x) + (y*y)))


def centered_quadratic(ISW, amplitude):
    """Quadratic potential barrier in the center of the well."""
    assert(isinstance(ISW, InfiniteSquareWell))

    def centered(x,y):
        x_width, y_width = ISW.well_x_width, ISW.well_y_width
        x_mid = (ISW.well_x_max - abs(ISW.well_x_min)) / 2.0
        y_mid = (ISW.well_y_max - abs(ISW.well_y_min)) / 2.0
        x_offset, y_offset = x - x_mid, y - y_mid
        if abs(x_offset) < 0.25 * x_width and abs(y_offset) < 0.25 * y_width:
            return amplitude * (x_offset**2 + y_offset**2)
        else:
            return 0
    return general_well(ISW, centered)


def square_barrier(ISW, amplitude):
    """Square-shaped potential barrier."""
    def sqb(x, y):
        x_width, y_width = ISW.well_x_width, ISW.well_y_width
        x_mid = (ISW.well_x_max - abs(ISW.well_x_min)) / 2.0
        y_mid = (ISW.well_y_max - abs(ISW.well_y_min)) / 2.0
        x_offset, y_offset = x - x_mid, y - y_mid
        if abs(x_offset) < 0.1 * x_width and abs(y_offset) < 0.1 * y_width:
            return amplitude
        else:
            return 0.0
    return general_well(ISW, sqb)


def square_plus_linear(ISW, amplitude):
    """Flat Potential that turns linear after a bit."""
    def spl(x, y):
        mid = (ISW.well_x_max + ISW.well_x_min) / 2.0
        offset = x - mid
        if offset < 0:
            return 0.0
        else:
            return amplitude * (offset)
    return general_well(ISW, spl)


def triangle_barrier(ISW, amplitude):
    """Triangle-Shaped Potential barrier."""
    def triangle(x):
        width = ISW.well_width
        mid = (ISW.well_max + ISW.well_min) / 2.0
        offset = x - mid
        if abs(offset) < 0.25*width:
            return -amplitude*(abs(offset) - 0.25*width)
        else:
            return 0.0
    return general_well(ISW, triangle)


def coupled_quadratic(ISW, amplitude):
    """Multiple quadratic potentials next to each other."""
    def cq(x):
        width = ISW.well_width
        mid = (ISW.well_max + ISW.well_min) / 2.0
        offset = x - mid
        if abs(offset) < 0.25 * width:
            return amplitude * (abs(offset) - 0.125*width) ** 2
        else:
            return 0.0
    return general_well(ISW, cq)


def kronig_penney(ISW, amplitude):
    """Kronig-Penney Potential to model solids."""
    def kp(x):
        # 3 barriers -> 0.25 0.5 0.75
        offset = x - ISW.well_min
        num_barriers = 5
        spacing = ISW.well_width / (num_barriers + 1)
        bar_wid = 2 / ((num_barriers + 1))
        # n bars of wid d equally spaced between min and max
        # positions of bars determined by width / (num_barriers+1)
        # if (x - left_lim) = w / (n+1), delta site
        for i in range(1, num_barriers+1):
            if i*spacing - bar_wid < offset and offset < i*spacing + bar_wid:
                return amplitude
        else:
            return 0
    return general_well(ISW, kp)

def hydrogen(ISW, amplitude):
    """Hydrogen Atom potential, ampltude equivalent to charge."""
    def h_atom(x, y):
        r = np.sqrt(x*x + y*y)
        return -amplitude / r if abs(r) > 1e-3 else h_atom(2e-3, 2e-3)
    return general_well(ISW, h_atom)

def lennard_jones(ISW, amplitude):
    """Lennard-Jones Potential, seen here: https://en.wikipedia.org/wiki/Lennard-Jones_potential."""
    def lj(x,y):
        r = np.sqrt(x*x+y*y)
        return amplitude * (r**-12 - r**-6) if abs(r) > 1e-3 else lj(2e-3, 2e-3)
    return general_well(ISW, lj)


class PotentialType(Enum):
    """Enumeration which contains all working potential types."""

    square = auto()              # WORKING
    linear = auto()              # WORKING
    quadratic = auto()           # WORKING
    centered_quadratic = auto()  # WORKING
    square_barrier = auto()      # WORKING
    square_plus_linear = auto()  # WORKING
    triangle_barrier = auto()    # BROKEN
    coupled_quadratic = auto()   # BROKEN
    kronig_penney = auto()       # BROKEN
    hydrogen = auto()            # WORKING
    lennard_jones = auto()       # WORKING

    def get_potential(self, ISW, amplitude):
        """From enum type, return the proper potential to compute."""
        assert(isinstance(ISW, InfiniteSquareWell))
        if self is PotentialType.square:
            return square(ISW, amplitude)
        elif self is PotentialType.linear:
            return linear(ISW, amplitude)
        elif self is PotentialType.quadratic:
            return quadratic(ISW, amplitude)
        elif self is PotentialType.centered_quadratic:
            return centered_quadratic(ISW, amplitude)
        elif self is PotentialType.square_barrier:
            return square_barrier(ISW, amplitude)
        elif self is PotentialType.square_plus_linear:
            return square_plus_linear(ISW, amplitude)
        elif self is PotentialType.triangle_barrier:
            return triangle_barrier(ISW, amplitude)
        elif self is PotentialType.coupled_quadratic:
            return coupled_quadratic(ISW, amplitude)
        elif self is PotentialType.kronig_penney:
            return kronig_penney(ISW, amplitude)
        elif self is PotentialType.hydrogen:
            return hydrogen(ISW, amplitude)
        elif self is PotentialType.lennard_jones:
            return lennard_jones(ISW, amplitude)

    def to_string(self):
        """Return a properly formatted Potential Name."""
        if self is PotentialType.square:
            return "Square Well"
        elif self is PotentialType.linear:
            return "Linear Well"
        elif self is PotentialType.quadratic:
            return "Quadratic Potential"
        elif self is PotentialType.centered_quadratic:
            return "Centered Quadratic Potential"
        elif self is PotentialType.square_barrier:
            return "Square Barrier"
        elif self is PotentialType.square_plus_linear:
            return "Square + Linear Well"
        elif self is PotentialType.triangle_barrier:
            return "Triangle Barrier"
        elif self is PotentialType.coupled_quadratic:
            return "Double Quadratic Potential"
        elif self is PotentialType.kronig_penney:
            return "Kronig-Penney Potential"
        elif self is PotentialType.hydrogen:
            return "Hydrogen Atom Potential"
        elif self is PotentialType.lennard_jones:
            return "Lennard-Jones Potential"
        else:
            return ""


if __name__ == "__main__":
    main()
