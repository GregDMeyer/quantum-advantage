
import cirq

class ToffoliSimulator:

    allowed_gates = [
        cirq.TOFFOLI,
        cirq.CNOT,
        cirq.X
    ]

    def __init__(self, c, qubits=None):
        self.circuit = c

        if qubits is None:
            qubits = list(c.all_qubits())
        self.qubits = qubits
        
        self.phase = 1

    def simulate(self, state):
        """
        Simulate the circuit, returning the result as an integer

        Arguments
        ---------

        state : dict
            A dictionary mapping qubit objects to their initial state. Will
            be modified to contain the result.
        """
        for moment in self.circuit:
            for op in moment:
                if not any(op.gate is g for g in self.allowed_gates):
                    raise ValueError(f"unsupported gate type '{op.gate}'")

                if all(state[q] for q in op.qubits[:-1]):
                    state[op.qubits[-1]] ^= 1

        return state


def int_to_state(x, reg):
    if x.bit_length() > len(reg):
        raise ValueError("integer too large for register")

    rtn = {
        q : (x>>i)&1 for i,q in enumerate(reg)
    }

    return rtn

def state_to_int(state, reg):
    rtn = 0
    for q in reg[::-1]:
        rtn <<= 1
        rtn |= state[q]
    return rtn
