
import cirq

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
    cout = ancillas.new()
    nB = len(B)
    for i, a in enumerate(A[::-1]):
        full_adder(circ, a, B[nB-i-1], cin, cout)
        ancillas.discard(cin)
        cin = cout
        cout = ancillas.new()
