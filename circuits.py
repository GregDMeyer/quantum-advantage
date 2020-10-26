
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

def copy_register(circ, A, B):
    """
    XOR the contents of A into B (thus copying them if B is empty)

    inputs:   A  B
    outputs:  A  A xor B
    """
    if len(A) > len(B):
        raise ValueError("register B must be at least as long as A")

    for a, b in zip(A, B):
        circ.append(cirq.CNOT(a, b))
    
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

def add_classical_int(circ, x, A, ancillas):
    """
    add the classical integer x to register A

    inputs:  A
    outputs: A+x
    """
    if x.bit_length() > len(A):
        raise ValueError("x too large to add to A")

    cin = ancillas.new()
    for a in A:
        cout = ancillas.new()
        if x&1 == 0:
            half_adder(circ, cin, a, cout)
        else:
            # it's a half adder + 1
            circ.append(cirq.X(cin))
            circ.append(cirq.X(a))
            half_adder(circ, cin, a, cout)
            circ.append(cirq.X(a))
            circ.append(cirq.X(cout))

        ancillas.discard(cin)
        cin = cout
        x >>= 1
        
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

# TODO: add flag if we know C is zero?
def karatsuba_mult(circ, A, B, C, ancillas, cutoff=4):
    """
    applies karatsuba multiplication

    inputs:   A  B  C
    outputs:  A  B  C+A*B
    """
    if len(C) < len(A)+len(B):
        raise ValueError("register C not large enough to store result")

    # base case
    if any(l < cutoff for l in (len(A), len(B))):
        schoolbook_mult(circ, A, B, C, ancillas)
        return
    
    AB_break = min(len(A), len(B))//2
    C_break = 2*AB_break

    A_low = A[:AB_break]
    B_low = B[:AB_break]
    C_low = C[:C_break]
    
    A_high = A[AB_break:]
    B_high = B[AB_break:]
    C_high = C[C_break:]

    # can't use C_low on this first one, because we need to carry
    # all the way to the end of C
    karatsuba_mult(circ, A_low,  B_low,  C, ancillas, cutoff)
    karatsuba_mult(circ, A_high, B_high, C_high, ancillas, cutoff)

    A_sum = ancillas.new_register(len(A_high)+1)
    B_sum = ancillas.new_register(len(B_high)+1)
    C_mid = ancillas.new_register(len(A_sum)+len(B_sum))

    copy_register(circ, A_low, A_sum)
    add_int(circ, A_high, A_sum, ancillas)
    
    copy_register(circ, B_low, B_sum)
    add_int(circ, B_high, B_sum, ancillas)

    # sum C_low and C_high, then negate them
    # in 2s complement, this is a bit flip + 1
    copy_register(circ, C_low, C_mid)
    add_int(circ, C_high, C_mid, ancillas)
    for c in C_mid:
        circ.append(cirq.X(c))
    add_classical_int(circ, 1, C_mid, ancillas)

    # finally add in the product
    karatsuba_mult(circ, A_sum, B_sum, C_mid, ancillas, cutoff)

    add_int(circ, C_mid, C[AB_break:], ancillas)
    
    ancillas.discard(A_sum)
    ancillas.discard(B_sum)
    ancillas.discard(C_mid)

# TODO: add a flag to skip extra carries if we know B is empty
def schoolbook_square(circ, A, B, ancillas):
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
                b_idx += 1
                cout = ancillas.new()
                half_adder(circ, cin, B[b_idx], cout)
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

# TODO: add flag if we know C is zero?
def karatsuba_square(circ, A, C, ancillas, cutoff=4):
    """
    applies karatsuba multiplication to square A

    inputs:   A  C
    outputs:  A  C+A*B
    """
    if len(C) < 2*len(A):
        raise ValueError("register C not large enough to store result")

    # base case
    if len(A) < cutoff:
        schoolbook_square(circ, A, C, ancillas)
        return
    
    A_break = len(A)//2
    A_low = A[:A_break]
    A_high = A[A_break:]

    C_mid = C[A_break+1:]   # +1 to handle factor of 2 in mult
    C_high = C[2*A_break:]

    karatsuba_square(circ, A_low,         C,      ancillas, cutoff)
    karatsuba_mult(circ,   A_low, A_high, C_mid,  ancillas, cutoff)
    karatsuba_square(circ, A_high,        C_high, ancillas, cutoff)
