"""Generate the hamiltonian given a basis and potential."""
from infinitesquarewell import InfiniteSquareWell


def mel(psil, V, psir, ISW):
    """Compute matrix element using average value theorem."""
    assert(isinstance(ISW, InfiniteSquareWell))
    # discrete inner product: < left | V | right >
    el = sum(l*v*r for (l, v, r) in zip(psil, V, psir))
    # readjust for avg val thm
    return float(ISW.well_width**2 * el / ISW.steps)

def compute_hamiltonian(V, ISW):
    """Compute discretized hamiltonian."""
    assert(isinstance(ISW, InfiniteSquareWell))
    hamiltonian = []
    # Note, the ket |i j> is a single eigenfunction
    for i in range(ISW.energy_eigenvals**2):
        row = []  # one row of the hamiltonian matrix
        for j in range(ISW.energy_eigenvals**2):
            psil, psir = ISW.basis_funcs[i], ISW.basis_funcs[j]
            el = mel(psil, V, psir, ISW)
            if i == j:  # diagonal elements get kinetic term
                el += ISW.eigenvals[i]
            row.append(el)
            el = 0.0
        hamiltonian.append(row)
    return hamiltonian
