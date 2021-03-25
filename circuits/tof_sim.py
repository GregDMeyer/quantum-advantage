
import cirq

class ToffoliSimulator:

    allowed_gates = {
        cirq.TOFFOLI,
        cirq.CNOT,
        cirq.X,
    }

    phase_error_gates = {
        cirq.Y,
        cirq.Z
    }

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
        phase = 1
        for moment in self.circuit:
            for op in moment:
                if op.gate in self.allowed_gates:
                    if all(state[q.name] for q in op.qubits[:-1]):
                        state[op.qubits[-1].name] ^= 1

                elif op.gate in self.phase_error_gates:
                    if state[op.qubits[0].name] == 1:
                        phase *= -1

                    # Y also flips the bit
                    if op.gate is cirq.Y:
                        state[op.qubits[-1].name] ^= 1

                else:
                    raise ValueError(f"unsupported gate type '{op.gate}'")

        # don't need to return state; it is modified in-place
        return phase


def int_to_state(x, reg):
    if x.bit_length() > len(reg):
        raise ValueError("integer too large for register")

    rtn = {
        q.name : (x>>i)&1 for i,q in enumerate(reg)
    }

    return rtn

def state_to_int(state, reg):
    rtn = 0
    for q in reg[::-1]:
        rtn <<= 1
        rtn |= state[q.name]
    return rtn
