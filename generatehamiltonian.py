"""Generate the hamiltonian given a basis and potential."""
import numpy as np
from infinitesquarewell import InfiniteSquareWell


def mel(psil, V, psir, ISW):
    """Compute matrix element using average value theorem."""
    assert(isinstance(ISW, InfiniteSquareWell))
    # discrete inner product: < left | V | right >
    el = np.mean(psil*V*psir)
    # readjust for avg val thm
    return float(el * ISW.well_width**2)

def compute_hamiltonian(V, ISW):
    """Compute discretized hamiltonian."""
    assert(isinstance(ISW, InfiniteSquareWell))
    hamiltonian = []
    for i in ISW.eigenvals.keys(): # note: i is a tuple 
        row = []  # one row of the hamiltonian matrix
        for j in ISW.eigenvals.keys():
            psil, psir = ISW.basis_funcs[i], ISW.basis_funcs[j]
            el = mel(psil, V, psir, ISW)
            if i == j:  # diagonal elements get kinetic term
                el += ISW.eigenvals[i]
            row.append(el)
            el = 0.0
        hamiltonian.append(row)
    return np.array(hamiltonian)
