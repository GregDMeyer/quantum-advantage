
import cirq

def half_adder(circ, A, B, Cout):
    """
    appends a half adder to circuit circ.

    inputs:  A  B    Cout
    outputs: A  A+B  Cout+(A+B)/2
    """
    circ.my_append(cirq.TOFFOLI(A, B, Cout))
    circ.my_append(cirq.CNOT(A, B))

def full_adder(circ, A, B, Cin, Cout):
    """
    appends a "full adder" circuit to circuit circ.

    inputs:   A  B         Cin  Cout
    outputs:  A  A+B+Cin%2 Cin  Cout+(A+B+C/2)

    the qubits A and Cin are measured in the H basis and zeroed
    after the adder is performed.
    """
    circ.my_append(cirq.TOFFOLI(A, B, Cout))
    circ.my_append(cirq.CNOT(A, B))
    circ.my_append(cirq.TOFFOLI(B, Cin, Cout))
    circ.my_append(cirq.CNOT(Cin, B))

def copy_register(circ, A, B):
    """
    XOR the contents of A into B (thus copying them if B is empty)

    inputs:   A  B
    outputs:  A  A xor B
    """
    if len(A) > len(B):
        raise ValueError("register B must be at least as long as A")

    for a, b in zip(A, B):
        circ.my_append(cirq.CNOT(a, b))

def add_int(circ, A, B, ancillas, allow_overflow=False):
    """
    adds the n qubit register A to m qubit register B, with n <= m

    inputs:   A  B
    outputs:  A  B+A % 2^m
    """
    if not allow_overflow and len(A) > len(B):
        raise ValueError("register A too long to add to register B")

    for i, a in enumerate(A):
        if i >= len(B):
            break

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

def add_classical_int(circ, x, A, ancillas, control=None):
    """
    add the classical integer x to register A

    inputs:  A
    outputs: A+x
    """
    if x.bit_length() > len(A):
        raise ValueError("x too large to add to A")

    if control is None:
        Xgate = cirq.X
    else:
        Xgate = lambda q: cirq.CNOT(control, q)

    cin = ancillas.new()
    for a in A:
        cout = ancillas.new()
        if x&1 == 0:
            half_adder(circ, cin, a, cout)
        else:
            # it's a half adder + 1
            circ.my_append(Xgate(cin))
            circ.my_append(Xgate(a))
            half_adder(circ, cin, a, cout)
            circ.my_append(Xgate(a))
            circ.my_append(Xgate(cout))

        ancillas.discard(cin)
        cin = cout
        x >>= 1

def lessthan_classical(circ, A, x, b, ancillas):
    """
    inputs:  A  b
    outputs: A  b+(A<x)
    """
    if x.bit_length() > len(A):
        # x is certainly larger than A
        circ.my_append(cirq.X(b))
        return

    # whether they are equal, starts as 1
    eq = ancillas.new()
    circ.my_append(cirq.X(eq))

    for i,a in reversed(list(enumerate(A))):
        eq_out = ancillas.new()
        circ.my_append(cirq.CNOT(eq, eq_out))
        if (x&(1<<i)) == 0:
            # if a is 1, they are not equal
            circ.my_append(cirq.TOFFOLI(eq, a, eq_out))
        else:
            # if a is 0 and they are equal so far, a is less and they are not equal
            circ.my_append(cirq.X(a))
            circ.my_append(cirq.TOFFOLI(eq, a, b))
            circ.my_append(cirq.TOFFOLI(eq, a, eq_out))
            circ.my_append(cirq.X(a))

        ancillas.discard(eq)
        eq = eq_out

    ancillas.discard(eq)

def schoolbook_mult(circ, A, B, C, ancillas, C_zero=False):
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
            circ.my_append(cirq.TOFFOLI(a, b, d))

            cout = ancillas.new()
            full_adder(circ, d, C[i+j], cin, cout)
            ancillas.discard(cin)
            cin = cout

            ancillas.discard(d)

        # finish performing the carries
        if not C_zero:
            for c in C[i+len(B):]:
                cout = ancillas.new()
                half_adder(circ, cin, c, cout)
                ancillas.discard(cin)
                cin = cout
        else:
            circ.my_append(cirq.CNOT(cin, C[i+len(B)]))

        ancillas.discard(cin)

# TODO: add flag if we know C is zero?
def karatsuba_mult(circ, A, B, C, ancillas, cutoff=None, C_zero=False):
    """
    applies karatsuba multiplication

    inputs:   A  B  C
    outputs:  A  B  C+A*B
    """
    if len(C) < len(A)+len(B):
        raise ValueError("register C not large enough to store result")

    if cutoff is None:
        _cutoff = 15  # found by optimization
    else:
        _cutoff = cutoff

    if _cutoff < 3:
        # otherwise we end up infinitely recursing
        raise ValueError("cutoff must be >= 3")

    # base case
    if any(l <= _cutoff for l in (len(A), len(B))):
        schoolbook_mult(circ, A, B, C, ancillas, C_zero=C_zero) # TODO
        return

    AB_break = min(len(A), len(B))//2

    A_low = A[:AB_break]
    B_low = B[:AB_break]

    A_high = A[AB_break:]
    B_high = B[AB_break:]

    C_mid  = ancillas.new_register(len(A_high)+len(B_high)+2)

    # just using C_mid here to avoid extra allocation of C_low
    karatsuba_mult(circ, A_low, B_low, C_mid, ancillas, cutoff, C_zero=True)
    if C_zero:
        copy_register(circ, C_mid[:2*AB_break], C[:2*AB_break])
    else:
        add_int(circ, C_mid, C, ancillas)  # C += C_low

    C_high_size = len(A_high)+len(B_high)
    if C_zero:
        C_high = C[2*AB_break:2*AB_break+C_high_size]
    else:
        C_high = ancillas.new_register(C_high_size)

    karatsuba_mult(circ, A_high, B_high, C_high, ancillas, cutoff, C_zero=True)
    add_int(circ, C_high, C_mid, ancillas)

    if not C_zero:
        add_int(circ, C_high, C[2*AB_break:], ancillas) # C += C_high
        ancillas.discard(C_high)

    A_sum = ancillas.new_register(len(A_high)+1)
    B_sum = ancillas.new_register(len(B_high)+1)

    copy_register(circ, A_low, A_sum)
    add_int(circ, A_high, A_sum, ancillas)

    copy_register(circ, B_low, B_sum)
    add_int(circ, B_high, B_sum, ancillas)

    # negate the sum of C_low and C_high (stored in C_mid)
    # in 2s complement, this is a bit flip + 1
    for c in C_mid:
        circ.my_append(cirq.X(c))
    add_classical_int(circ, 1, C_mid, ancillas)

    # finally add in the product of A_sum and B_sum
    karatsuba_mult(circ, A_sum, B_sum, C_mid, ancillas, cutoff)
    add_int(circ, C_mid, C[AB_break:], ancillas)

    ancillas.discard(A_sum)
    ancillas.discard(B_sum)
    ancillas.discard(C_mid)

def schoolbook_classical_mult(circ, a, B, C, ancillas, allow_overflow=False, C_zero=False):
    """
    applies schoolbook multiplication, with classical input a

    inputs:   B  C
    outputs:  B  C+a*B
    """
    len_a = a.bit_length() if a>=0 else len(C)

    if not allow_overflow and len(C) < len_a+len(B):
        raise ValueError("register C not large enough to store result")

    for i in range(len_a):
        if a&1:
            cin = ancillas.new()
            for j,b in enumerate(B):
                if i+j >= len(C):
                    break
                cout = ancillas.new()
                full_adder(circ, b, C[i+j], cin, cout)
                ancillas.discard(cin)
                cin = cout

            # finish performing the carries
            if not C_zero:
                for c in C[i+len(B):]:
                    cout = ancillas.new()
                    half_adder(circ, cin, c, cout)
                    ancillas.discard(cin)
                    cin = cout
            elif len(C) > i+len(B):
                circ.my_append(cirq.CNOT(cin, C[i+len(B)]))

            ancillas.discard(cin)

        a >>= 1

def karatsuba_classical_mult(circ, a, B, C, ancillas, cutoff=None, allow_overflow=False, C_zero=False):
    """
    applies karatsuba multiplication with one classical input

    inputs:   B  C
    outputs:  B  C+a*B
    """
    len_a = a.bit_length() if a>=0 else len(C)

    if not allow_overflow and len(C) < len_a+len(B):
        raise ValueError("register C not large enough to store result")

    if cutoff is None:
        _cutoff = 24
    else:
        _cutoff = cutoff

    if _cutoff < 4:
        # otherwise we end up infinitely recursing
        raise ValueError("cutoff must be >= 4")

    # if C has the same length as a and B, it doesn't make sense to do the full karatsuba
    if allow_overflow and len(C)>4:
        if len_a > len(B):
            B_break = len(C) // 2
            A_break = len(C) - B_break
        else:
            A_break = len(C) // 2
            B_break = len(C) - A_break

        a_low = a & ((1<<A_break)-1)
        B_low = B[:B_break]

        a_high = a>>A_break
        B_high = B[B_break:]

        karatsuba_classical_mult(circ, a_low, B_low, C, ancillas, cutoff, C_zero=C_zero)
        karatsuba_classical_mult(circ, a_low, B_high, C[B_break:], ancillas, cutoff, allow_overflow=True)
        karatsuba_classical_mult(circ, a_high, B_low, C[A_break:], ancillas, cutoff, allow_overflow=True)
        return

    # base case
    if any(l <= _cutoff for l in (len_a, len(B), len(C))):
        schoolbook_classical_mult(circ, a, B, C, ancillas, allow_overflow=allow_overflow, C_zero=C_zero)
        return

    AB_break = min(len_a, len(B))//2

    a_low = a & ((1<<AB_break)-1)
    B_low = B[:AB_break]

    a_high = a>>AB_break
    B_high = B[AB_break:]

    C_mid  = ancillas.new_register(a_high.bit_length()+len(B_high)+2)

    # just using C_mid here to avoid extra allocation of C_low
    karatsuba_classical_mult(circ, a_low, B_low, C_mid, ancillas, cutoff, C_zero=True)
    if C_zero:
        copy_register(circ, C_mid[:2*AB_break], C[:2*AB_break])
    else:
        add_int(circ, C_mid, C, ancillas)  # C += C_low

    C_high_size = a_high.bit_length()+len(B_high)
    if C_zero:
        C_high = C[2*AB_break:2*AB_break+C_high_size]
    else:
        C_high = ancillas.new_register(C_high_size)

    karatsuba_classical_mult(circ, a_high, B_high, C_high, ancillas, cutoff, C_zero=True)
    add_int(circ, C_high, C_mid, ancillas)

    if not C_zero:
        add_int(circ, C_high, C[2*AB_break:], ancillas) # C += C_high
        ancillas.discard(C_high)

    a_sum = a_low + a_high

    B_sum = ancillas.new_register(len(B_high)+1)
    copy_register(circ, B_low, B_sum)
    add_int(circ, B_high, B_sum, ancillas)

    # negate the sum of C_low and C_high (stored in C_mid)
    # in 2s complement, this is a bit flip + 1
    for c in C_mid:
        circ.my_append(cirq.X(c))
    add_classical_int(circ, 1, C_mid, ancillas)

    # finally add in the product of A_sum and B_sum
    karatsuba_classical_mult(circ, a_sum, B_sum, C_mid, ancillas, cutoff)

    add_int(circ, C_mid, C[AB_break:], ancillas)

    ancillas.discard(B_sum)
    ancillas.discard(C_mid)

def schoolbook_square(circ, A, B, ancillas, C_zero=False):
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
                circ.my_append(cirq.TOFFOLI(A[i], A[j], a))

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
        if not C_zero:
            for b in B[b_idx:]:
                cout = ancillas.new()
                half_adder(circ, cin, b, cout)
                ancillas.discard(cin)
                cin = cout
        elif len(B) > b_idx:
            circ.my_append(cirq.CNOT(cin, B[b_idx]))

        ancillas.discard(cin)

# TODO: add flag if we know C is zero
def karatsuba_square(circ, A, C, ancillas, cutoff=None, C_zero=False):
    """
    applies karatsuba multiplication to square A

    inputs:   A  C
    outputs:  A  C+A*B
    """
    if len(C) < 2*len(A):
        raise ValueError("register C not large enough to store result")

    if cutoff is None:
        _cutoff = 15
    else:
        _cutoff = cutoff

    # base case
    if len(A) <= _cutoff:
        schoolbook_square(circ, A, C, ancillas, C_zero=C_zero)
        return

    A_break = len(A)//2
    A_low = A[:A_break]
    A_high = A[A_break:]

    if C_zero:
        C_low = C[:2*A_break]
    else:
        C_low = ancillas.new_register(2*A_break)

    karatsuba_square(circ, A_low, C_low, ancillas, cutoff, C_zero=True)

    if not C_zero:
        add_int(circ, C_low, C, ancillas)
        ancillas.discard(C_low)
        C_high = ancillas.new_register(2*(len(A)-A_break))
    else:
        C_high = C[2*A_break:]

    karatsuba_square(circ, A_high, C_high, ancillas, cutoff, C_zero=True)

    if not C_zero:
        add_int(circ, C_high, C[2*A_break:], ancillas)
        ancillas.discard(C_high)

    C_mid = ancillas.new_register(len(A))
    karatsuba_mult(circ, A_low, A_high, C_mid, ancillas, cutoff, C_zero=True)
    add_int(circ, C_mid, C[A_break+1:], ancillas)
    ancillas.discard(C_mid)

def montgomery_reduce(circ, T, ancillas, N, mult=None):

    if len(T) < 2*N.bit_length() + 1:
        raise ValueError("T too small for montgomery reduction")

    if not N&1:
        raise ValueError("N must be odd")

    if mult is None:
        mult = karatsuba_classical_mult

    r = N.bit_length()  # R is closest power of 2 greater than N
    R = 1 << r
    _, Np = extended_gcd(R, N)

    m = ancillas.new_register(r)
    #print('m=Np*T')
    mult(circ, Np, T[:r], m, ancillas, allow_overflow=True, C_zero=True)
    #print('T+=Nm')
    mult(circ, N, m, T, ancillas)
    ancillas.discard(m)

    b = ancillas.new()
    lessthan_classical(circ, T[r:], N, b, ancillas)
    circ.my_append(cirq.X(b))
    #print('T-=N')
    add_classical_int(circ, -N, T[r:], ancillas, b)

    return R

# adapted from https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
def extended_gcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1

    while r != 0:
        quotient = old_r // r
        old_r, r = r, (old_r - quotient * r)
        old_s, s = s, (old_s - quotient * s)
        old_t, t = t, (old_t - quotient * t)

    # we want a*x - b*y = 1 and 0 < x < b and 0 < y < a
    if old_s < 0:
        old_s += b

    if old_t > 0:
        old_t -= a

    return old_s, -old_t

def x2_mod_N(circ, N, x, y, ancillas, mult_type):
    if mult_type == "karatsuba":
        karatsuba_square(circ, x, y, ancillas, C_zero=True)
        R = montgomery_reduce(circ, y, ancillas, N, mult=karatsuba_classical_mult)

    elif mult_type == "schoolbook":
        schoolbook_square(circ, x, y, ancillas, C_zero=True)
        R = montgomery_reduce(circ, y, ancillas, N, mult=schoolbook_classical_mult)

    else:
        raise ValueError(f"unknown mult type '{mult_type}'")

    return R
