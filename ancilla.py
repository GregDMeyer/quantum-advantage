
import cirq
from tof_sim import int_to_state

class Ancilla(cirq.NamedQubit):

    def __init__(self, idx):
        cirq.NamedQubit.__init__(self, f'anc{idx}')
        self._isgarbage = False

    @property
    def isgarbage(self):
        return self._isgarbage

    def discard(self):
        self._isgarbage = True

class AncillaManager:

    def __init__(self):
        self._qubits = []
        
    def new(self):
        q = Ancilla(len(self._qubits))
        self._qubits.append(q)
        return q

    def new_register(self, n):
        rtn = [Ancilla(len(self._qubits)+i) for i in range(n)]
        self._qubits += rtn
        return rtn
    
    def discard(self, qubit):
        qubit.discard()

    def all(self):
        return self._qubits

    def init_state(self):
        return int_to_state(0, self._qubits)
        
    def __len__(self):
        return len(self._qubits)
