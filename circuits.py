
import cirq
from ancilla import Ancilla

def half_adder(circ, A, B, Cout):
    """
    appends a half adder to circuit circ.

    inputs:  A  B    Cout
    outputs: A  A+B  Cout+(A+B)/2
    """
    circ.append(cirq.TOFFOLI(A, B, Cout))
    circ.append(cirq.CNOT(A, B))

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
    outputs:  A  B+A % 2^m
    """
    if len(A) > len(B):
        raise ValueError("register A too long to add to register B")

    for i, a in enumerate(A):
        cout = ancillas.new()
        if i == 0:
            half_adder(circ, a, B[i], cout)
        else:
            full_adder(circ, a, B[i], cin, cout)
            ancillas.discard(cin)
        cin = cout

    # need to carry to the end of B
    for b in B[len(A):]:
        cout = ancillas.new()
        half_adder(circ, cin, b, cout)
        ancillas.discard(cin)
        cin = cout

# TODO: add a flag to skip extra carries if we know B is empty
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

        # finish performing the carries
        b_idx += 1
        for b in B[b_idx:]:
            cout = ancillas.new()
            half_adder(circ, cin, b, cout)
            ancillas.discard(cin)
            cin = cout

        ancillas.discard(cin)

def schoolbook_mult(circ, A, B, C, ancillas):
    """
    applies schoolbook multiplication

    inputs:   A  B  C
    outputs:  A  B  C+A*B
    """
    if len(C) < len(A)+len(B):
        raise ValueError("register C not large enough to store result")

    for i,a in enumerate(A):
        cin = ancillas.new()
        for j,b in enumerate(B):
            d = ancillas.new()
            circ.append(cirq.TOFFOLI(a, b, d))
                
            cout = ancillas.new()
            full_adder(circ, d, C[i+j], cin, cout)
            ancillas.discard(cin)
            cin = cout

            ancillas.discard(d)

        # finish performing the carries
        for c in C[i+len(B):]:
            cout = ancillas.new()
            half_adder(circ, cin, c, cout)
            ancillas.discard(cin)
            cin = cout

        ancillas.discard(cin)
