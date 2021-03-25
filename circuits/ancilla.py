
import cirq
from tof_sim import int_to_state

class Ancilla(cirq.NamedQubit):

    def __init__(self, idx):
        cirq.NamedQubit.__init__(self, f'anc{idx}')
        self._isgarbage = False

    @property
    def isgarbage(self):
        return self._isgarbage

    def _discard(self):
        self._isgarbage = True

class AncillaManager:

    def __init__(self):
        self._qubits = []
        self._n_active = 0
        self._max_ancillas = 0

    def new(self):
        q = Ancilla(len(self._qubits))
        self._qubits.append(q)
        self.n_active += 1
        return q

    def new_register(self, n):
        rtn = [Ancilla(len(self._qubits)+i) for i in range(n)]
        self._qubits += rtn
        self.n_active += n
        return rtn

    def discard(self, qubits):
        if isinstance(qubits, Ancilla):
            qubits = [qubits]

        for qubit in qubits:
            qubit._discard()

        self.n_active -= len(qubits)

    def all(self):
        return self._qubits

    def init_state(self):
        return int_to_state(0, self._qubits)

    def all_discarded(self):
        return self.n_active == 0

    def max_ancilla_usage(self):
        return self._max_ancillas

    @property
    def n_active(self):
        return self._n_active

    @n_active.setter
    def n_active(self, val):
        self._n_active = val
        if self.n_active > self._max_ancillas:
            self._max_ancillas = self.n_active

    def __len__(self):
        return len(self._qubits)
