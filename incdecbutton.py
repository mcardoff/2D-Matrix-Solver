"""Create blank class to hold functions attached to buttons in main."""


class IncDecButton:
    """Class to track inc and dec of eigenfunctions."""

    def __init__(self, subfig, canvas, x, y, funcs, V):
        """Initialize necessary variables, attached to tkinter window."""
        # TODO Implement extrema calculation
        self.selector = 0                   # which one do we show
        self.subfig = subfig                # where to we put it
        self.canvas = canvas                # where do we draw
        self.x, self.y = x, y               # same x vals
        self.funcs = funcs                  # library of functions
        self.potential = V                  # potential
        self.energy_eigenvals = len(funcs)  # set in main
        self.calc_extrema()                 # x and y limits

    # plot the function func, clearing previous plot and resetting limits
    def replot(self, func, clear):
        """Clear the plot and plot func in its place."""
        if clear:
            self.subfig.clear()

        self.subfig.set_xlabel("x")
        self.subfig.set_ylabel("y")
        self.subfig.set_zlabel("Wavefunction")
        self.subfig.set_xlim(self.x_min, self.x_max)
        self.subfig.set_ylim(self.y_min, self.y_max)
        self.subfig.set_zlim(-self.max_val, self.max_val)
        self.subfig.plot_surface(self.x, self.y, func, cmap='viridis')
        self.canvas.draw()

    def inc_selector(self):
        """Show 'next' plot, loop to beginning at the end."""
        if self.selector < len(self.funcs)-1:
            self.selector += 1
        else:
            self.selector = 0

        self.replot(self.funcs[self.selector], True)

    def dec_selector(self):
        """Show 'previous' plot, loop to end if first."""
        if self.selector > 0:
            self.selector -= 1
        else:
            self.selector = len(self.funcs)-1

        self.replot(self.funcs[self.selector], True)

    def plot_potential(self):
        """Show potential on top of current plot."""
        self.replot(self.potential, False)

    def init_plot(self):
        """Show first plot in the sequence, ignore value of selector."""
        self.replot(self.funcs[0], True)

    def calc_extrema(self):
        """Calculate and recalculate maxes."""
        flat_funcs = [item for sublist in self.funcs for subsublist in sublist for item in subsublist]
        self.max_val = max(map(abs, flat_funcs))            # z limits
        self.x_max, self.x_min = self.x[-1,-1], self.x[0,0] # Set in ISW
        self.y_max, self.y_min = self.y[-1,-1], self.y[0,0] # Set in ISW

    def update_vals(self, x, y, funcs, V):
        """Set values of x, funcs, V and recalculate maxes."""
        self.x, self.y = x, y
        self.funcs = funcs
        self.potential = V
        self.calc_extrema()
        self.init_plot()
