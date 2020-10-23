
import cirq
from ancilla import Ancilla

def full_adder(circ, A, B, Cin, Cout):
    """
    appends a "full adder" circuit to circuit circ.

    inputs:   A  B         Cin  Cout
    outputs:  A  A+B+Cin%2 Cin  Cout+(A+B+C/2)

    the qubits A and Cin are measured in the H basis and zeroed 
    after the adder is performed.
    """
    circ.append(cirq.TOFFOLI(A, B, Cout))
    circ.append(cirq.CNOT(A, B))
    circ.append(cirq.TOFFOLI(B, Cin, Cout))
    circ.append(cirq.CNOT(Cin, B))

def add_int(circ, A, B, ancillas):
    """
    adds the n qubit register A to m qubit register B, with n <= m

    inputs:   A  B
    outputs:  A  B+A % 2^n
    """
    if len(A) > len(B):
        raise ValueError("register A too long to add to register B")

    # TODO: use a half adder for first one
    cin = ancillas.new()
    nB = len(B)
    for i, a in enumerate(A):
        cout = ancillas.new()
        full_adder(circ, a, B[i], cin, cout)
        ancillas.discard(cin)
        cin = cout

def square(circ, A, B, ancillas):
    """
    squares the integer in the A register, and adds the result to
    the B register.

    inputs:   A  B
    outputs:  A  B+A^2
    """
    if len(B) < 2*len(A):
        raise ValueError("register B not long enough to store square")

    for i in range(len(A)):
        cin = ancillas.new()
        for j in range(i, len(A)):
            if i == j:
                a = A[i]
            else:
                a = ancillas.new()
                circ.append(cirq.TOFFOLI(A[i], A[j], a))

            b_idx = i+j+(i!=j)
                
            cout = ancillas.new()
            full_adder(circ, a, B[b_idx], cin, cout)
            ancillas.discard(cin)
            cin = cout

            if i == j:
                # we need to carry an extra bit when i == j, because we are
                # skipping one
                cout = ancillas.new()
                circ.append(cirq.TOFFOLI(B[b_idx+1], cin, cout))
                circ.append(cirq.CNOT(cin, B[b_idx+1]))
                ancillas.discard(cin)
                cin = cout
            else:
                ancillas.discard(a)

        # record the last carry bit
        if i < len(A) - 1:
            circ.append(cirq.CNOT(cin, B[len(A)+i+1]))

        ancillas.discard(cin)
