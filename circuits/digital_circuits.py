
from itertools import repeat
import cirq


def x2_mod_N(N, x, y, ancillas, mult_type, threes=0):

    factor = 3**threes

    # make sure we have enough room
    xmax = factor * N
    x_len = xmax.bit_length()  # could be -1?
    y_len = (2*xmax + 1).bit_length()
    if len(x) < x_len or len(y) < y_len:
        raise ValueError('registers not large enough')

    three_gen = (times_three(x, ancillas) for _ in range(threes))

    if mult_type == "karatsuba":
        mult_gen = karatsuba_square(x, y, ancillas, C_zero=True)
        R, reduce_gen = montgomery_reduce(y, ancillas, factor**2 * N,
                                          mult=karatsuba_classical_mult)

    elif mult_type == "schoolbook":
        mult_gen = schoolbook_square(x, y, ancillas, C_zero=True)
        R, reduce_gen = montgomery_reduce(y, ancillas, factor**2 * N,
                                          mult=schoolbook_classical_mult)

    else:
        raise ValueError(f"unknown mult type '{mult_type}'")

    return R, (three_gen, mult_gen, reduce_gen)


# TODO: can we trim these a bit? I think so.
def get_registers(n, factor):
    extra_bits = factor.bit_length()
    x_reg = cirq.NamedQubit.range(n+extra_bits, prefix="x")
    y_reg = cirq.NamedQubit.range(2*(n+2*extra_bits)+1, prefix="y")
    return x_reg, y_reg


def add_int(A, B, ancillas, allow_overflow=False):
    """
    adds the n qubit register A to m qubit register B, with n <= m

    inputs:   A  B
    outputs:  A  B+A % 2^m
    """
    if not allow_overflow and len(A) > len(B):
        raise ValueError("register A too long to add to register B")

    cin = None
    for i, a in enumerate(A):
        if i >= len(B):
            break

        cout = ancillas.new()
        if cin is None:
            yield half_adder(a, B[i], cout)
        else:
            yield full_adder(a, B[i], cin, cout)
            ancillas.discard(cin)
        cin = cout

    # need to carry to the end of B
    for b in B[len(A):]:
        cout = ancillas.new()
        yield half_adder(cin, b, cout)
        ancillas.discard(cin)
        cin = cout

    ancillas.discard(cin)


def times_three(A, ancillas):
    """
    multiply the value in register A by 3, in-place
    """
    # can't really check for appropriate register size
    # user just needs to make sure that A is big enough
    # (should have room for two extra bits)

    cin = None
    prev = None
    for a in A[:-1]:
        cout = ancillas.new()
        new_prev = ancillas.new()
        yield cirq.CNOT(a, new_prev)
        if cin is not None:
            yield full_adder(prev, a, cin, cout)
            ancillas.discard(cin)
            ancillas.discard(prev)
        prev = new_prev
        cin = cout

    # carry to the end---either prev or cin may be 1, but not both
    yield cirq.CNOT(cin, A[-1])
    yield cirq.CNOT(prev, A[-1])
    ancillas.discard(cin)
    ancillas.discard(prev)


def add_classical_int(x, A, ancillas, control=None):
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
            yield half_adder(cin, a, cout)
        else:
            # it's a half adder + 1
            yield Xgate(cin)
            yield Xgate(a)
            yield half_adder(cin, a, cout)
            yield Xgate(a)
            yield Xgate(cout)

        ancillas.discard(cin)
        cin = cout
        x >>= 1

    ancillas.discard(cin)


def lessthan_classical(A, x, b, ancillas):
    """
    inputs:  A  b
    outputs: A  b+(A<x)
    """
    if x.bit_length() > len(A):
        # x is certainly larger than A
        yield cirq.X(b)
        return

    # whether they are equal, starts as 1
    eq = ancillas.new()
    yield cirq.X(eq)

    for i, a in reversed(list(enumerate(A))):
        eq_out = ancillas.new()
        yield cirq.CNOT(eq, eq_out)
        if x & (1<<i) == 0:
            # if a is 1, they are not equal
            yield cirq.TOFFOLI(eq, a, eq_out)
        else:
            # if a is 0 and they are equal so far, a is less and
            # they are not equal
            yield cirq.X(a)
            yield cirq.TOFFOLI(eq, a, b)
            yield cirq.TOFFOLI(eq, a, eq_out)
            yield cirq.X(a)

        ancillas.discard(eq)
        eq = eq_out

    ancillas.discard(eq)


def schoolbook_mult(A, B, C, ancillas, C_zero=False):
    """
    applies schoolbook multiplication

    inputs:   A  B  C
    outputs:  A  B  C+A*B
    """
    check_mult_sizes(A, B, C)

    for i, a in enumerate(A):
        cin = ancillas.new()
        for j, b in enumerate(B):
            d = ancillas.new()
            yield cirq.TOFFOLI(a, b, d)

            cout = ancillas.new()
            yield full_adder(d, C[i+j], cin, cout)
            ancillas.discard(cin)
            cin = cout

            ancillas.discard(d)

        # finish performing the carries
        if not C_zero:
            for c in C[i+len(B):]:
                cout = ancillas.new()
                yield half_adder(cin, c, cout)
                ancillas.discard(cin)
                cin = cout
        else:
            yield cirq.CNOT(cin, C[i+len(B)])

        ancillas.discard(cin)


def schoolbook_classical_mult(a, B, C, ancillas,
                              allow_overflow=False, C_zero=False):
    """
    applies schoolbook multiplication, with classical input a

    inputs:   B  C
    outputs:  B  C+a*B
    """
    if not allow_overflow:
        check_mult_sizes(a, B, C)

    # TODO: should just use add here, instead of reimplementing it?
    for i in range(a.bit_length() if a >= 0 else len(C)):
        if a&1:
            cin = ancillas.new()
            for j, b in enumerate(B):
                if i+j >= len(C):
                    break
                cout = ancillas.new()
                yield full_adder(b, C[i+j], cin, cout)
                ancillas.discard(cin)
                cin = cout

            # finish performing the carries
            if not C_zero:
                for c in C[i+len(B):]:
                    cout = ancillas.new()
                    yield half_adder(cin, c, cout)
                    ancillas.discard(cin)
                    cin = cout
            elif len(C) > i+len(B):
                yield cirq.CNOT(cin, C[i+len(B)])

            ancillas.discard(cin)

        a >>= 1


def schoolbook_square(A, B, ancillas, C_zero=False):
    """
    squares the integer in the A register, and adds the result to
    the B register.

    inputs:   A  B
    outputs:  A  B+A^2
    """
    check_mult_sizes(A, A, B)

    for i in range(len(A)):
        cin = ancillas.new()
        for j in range(i, len(A)):
            if i == j:
                a = A[i]
            else:
                a = ancillas.new()
                yield cirq.TOFFOLI(A[i], A[j], a)

            b_idx = i + j + (i!=j)

            cout = ancillas.new()
            yield full_adder(a, B[b_idx], cin, cout)
            ancillas.discard(cin)
            cin = cout

            if i == j:
                # we need to carry an extra bit when i == j, because we are
                # skipping one
                b_idx += 1
                cout = ancillas.new()
                yield half_adder(cin, B[b_idx], cout)
                ancillas.discard(cin)
                cin = cout
            else:
                ancillas.discard(a)

        # finish performing the carries
        b_idx += 1
        if not C_zero:
            for b in B[b_idx:]:
                cout = ancillas.new()
                yield half_adder(cin, b, cout)
                ancillas.discard(cin)
                cin = cout
        elif len(B) > b_idx:
            yield cirq.CNOT(cin, B[b_idx])

        ancillas.discard(cin)


def karatsuba_mult(A, B, C, ancillas, cutoff=None, C_zero=False):
    """
    applies karatsuba multiplication

    inputs:   A  B  C
    outputs:  A  B  C+A*B
    """
    check_mult_sizes(A, B, C)

    # default found via optimize_cutoff.py
    _cutoff = get_cutoff(cutoff, default=15)

    # base case
    if any(l <= _cutoff for l in (len(A), len(B))):
        yield schoolbook_mult(A, B, C, ancillas, C_zero=C_zero)
        return

    AB_break = min(len(A), len(B))//2

    A_low = A[:AB_break]
    B_low = B[:AB_break]

    A_high = A[AB_break:]
    B_high = B[AB_break:]

    C_mid = ancillas.new_register(len(A_high)+len(B_high)+2)

    # just using C_mid here to avoid extra allocation of C_low
    yield karatsuba_mult(A_low, B_low, C_mid, ancillas, cutoff, C_zero=True)
    if C_zero:
        yield copy_register(C_mid[:2*AB_break], C[:2*AB_break])
    else:
        yield add_int(C_mid, C, ancillas)  # C += C_low

    C_high_size = len(A_high)+len(B_high)
    if C_zero:
        C_high = C[2*AB_break:2*AB_break+C_high_size]
    else:
        C_high = ancillas.new_register(C_high_size)

    yield karatsuba_mult(A_high, B_high, C_high, ancillas, cutoff, C_zero=True)
    yield add_int(C_high, C_mid, ancillas)

    if not C_zero:
        # C += C_high
        yield add_int(C_high, C[2*AB_break:], ancillas)
        ancillas.discard(C_high)

    A_sum = ancillas.new_register(len(A_high)+1)
    B_sum = ancillas.new_register(len(B_high)+1)

    yield copy_register(A_low, A_sum)
    yield add_int(A_high, A_sum, ancillas)

    yield copy_register(B_low, B_sum)
    yield add_int(B_high, B_sum, ancillas)

    # negate the sum of C_low and C_high (stored in C_mid)
    # in 2s complement, this is a bit flip + 1
    for c in C_mid:
        yield cirq.X(c)
    yield add_classical_int(1, C_mid, ancillas)

    # finally add in the product of A_sum and B_sum
    yield karatsuba_mult(A_sum, B_sum, C_mid, ancillas, cutoff)
    yield add_int(C_mid, C[AB_break:], ancillas)

    ancillas.discard(A_sum)
    ancillas.discard(B_sum)
    ancillas.discard(C_mid)

def karatsuba_classical_mult(a, B, C, ancillas, cutoff=None, allow_overflow=False, C_zero=False):
    """
    applies karatsuba multiplication with one classical input

    inputs:   B  C
    outputs:  B  C+a*B
    """
    if not allow_overflow:
        check_mult_sizes(a, B, C)

    _cutoff = get_cutoff(cutoff, default=24)

    # if C has the same length as a and B, it doesn't
    # make sense to do the full karatsuba
    len_a = a.bit_length() if a >= 0 else len(C)
    if allow_overflow and len(C) > 4:
        if len_a > len(B):
            B_break = len(C) // 2
            A_break = len(C) - B_break
        else:
            A_break = len(C) // 2
            B_break = len(C) - A_break

        a_low = a & ((1<<A_break)-1)
        B_low = B[:B_break]

        a_high = a >> A_break
        B_high = B[B_break:]

        yield karatsuba_classical_mult(a_low, B_low, C,
                                       ancillas, cutoff, C_zero=C_zero)

        yield karatsuba_classical_mult(a_low, B_high, C[B_break:],
                                       ancillas, cutoff, allow_overflow=True)

        yield karatsuba_classical_mult(a_high, B_low, C[A_break:],
                                       ancillas, cutoff, allow_overflow=True)
        return

    # base case
    if any(l <= _cutoff for l in (len_a, len(B), len(C))):
        yield schoolbook_classical_mult(a, B, C, ancillas,
                                        allow_overflow=allow_overflow,
                                        C_zero=C_zero)
        return

    AB_break = min(len_a, len(B))//2

    a_low = a & ((1<<AB_break)-1)
    B_low = B[:AB_break]

    a_high = a>>AB_break
    B_high = B[AB_break:]

    C_mid = ancillas.new_register(a_high.bit_length()+len(B_high)+2)

    # just using C_mid here to avoid extra allocation of C_low
    yield karatsuba_classical_mult(a_low, B_low, C_mid,
                                   ancillas, cutoff, C_zero=True)
    if C_zero:
        yield copy_register(C_mid[:2*AB_break], C[:2*AB_break])
    else:
        yield add_int(C_mid, C, ancillas)  # C += C_low

    C_high_size = a_high.bit_length()+len(B_high)
    if C_zero:
        C_high = C[2*AB_break:2*AB_break+C_high_size]
    else:
        C_high = ancillas.new_register(C_high_size)

    yield karatsuba_classical_mult(a_high, B_high, C_high,
                                   ancillas, cutoff, C_zero=True)
    yield add_int(C_high, C_mid, ancillas)

    if not C_zero:
        # C += C_high
        yield add_int(C_high, C[2*AB_break:], ancillas)
        ancillas.discard(C_high)

    a_sum = a_low + a_high

    B_sum = ancillas.new_register(len(B_high)+1)
    yield copy_register(B_low, B_sum)
    yield add_int(B_high, B_sum, ancillas)

    # negate the sum of C_low and C_high (stored in C_mid)
    # in 2s complement, this is a bit flip + 1
    for c in C_mid:
        yield cirq.X(c)
    yield add_classical_int(1, C_mid, ancillas)

    # finally add in the product of A_sum and B_sum
    yield karatsuba_classical_mult(a_sum, B_sum, C_mid, ancillas, cutoff)

    yield add_int(C_mid, C[AB_break:], ancillas)

    ancillas.discard(B_sum)
    ancillas.discard(C_mid)


def karatsuba_square(A, C, ancillas, cutoff=None, C_zero=False):
    """
    applies karatsuba multiplication to square A

    inputs:   A  C
    outputs:  A  C+A*A
    """
    check_mult_sizes(A, A, C)
    _cutoff = get_cutoff(cutoff, default=15)

    # base case
    if len(A) <= _cutoff:
        yield schoolbook_square(A, C, ancillas, C_zero=C_zero)
        return

    A_break = len(A)//2
    A_low = A[:A_break]
    A_high = A[A_break:]

    if C_zero:
        C_low = C[:2*A_break]
    else:
        C_low = ancillas.new_register(2*A_break)

    yield karatsuba_square(A_low, C_low, ancillas, cutoff, C_zero=True)

    if not C_zero:
        yield add_int(C_low, C, ancillas)
        ancillas.discard(C_low)
        C_high = ancillas.new_register(2*(len(A)-A_break))
    else:
        C_high = C[2*A_break:]

    yield karatsuba_square(A_high, C_high, ancillas, cutoff, C_zero=True)

    if not C_zero:
        yield add_int(C_high, C[2*A_break:], ancillas)
        ancillas.discard(C_high)

    C_mid = ancillas.new_register(len(A))
    yield karatsuba_mult(A_low, A_high, C_mid, ancillas, cutoff, C_zero=True)
    yield add_int(C_mid, C[A_break+1:], ancillas)
    ancillas.discard(C_mid)


def montgomery_reduce(T, ancillas, N, mult=None):
    """
    Perform montgomery reduction on the register T, for some modulus N.

    Takes the 2n digit value T and replaces it with the n-digit
    value (T/R) mod N, followed by n zeros.

    Returns the (classical) value R, which should be multiplied by
    T/R mod N after measurement to ultimately yield T mod N classically.

    The T register should also have one extra leading qubit that is
    initialized to zero (and is returned as zero).
    """
    if len(T) < 2*N.bit_length() + 1:
        raise ValueError("T too small for montgomery reduction")

    # required for extended gcd to work
    if not N & 1:
        raise ValueError("N must be odd")

    if mult is None:
        mult = karatsuba_classical_mult

    r = N.bit_length()  # R is closest power of 2 greater than N
    R = 1 << r
    _, Np = extended_gcd(R, N)

    m = ancillas.new_register(r)

    # need to put all the generators in a list so we can return R also
    ops = []

    # m = Np*T
    ops.append(mult(Np, T[:r], m, ancillas, allow_overflow=True, C_zero=True))

    # T += Nm
    ops.append(mult(N, m, T, ancillas))
    ancillas.discard(m)

    # if T >= N
    b = ancillas.new()
    ops.append(lessthan_classical(T[r:], N, b, ancillas))
    ops.append(cirq.X(b))

    # then T -= N
    ops.append(add_classical_int(-N, T[r:], ancillas, b))
    ancillas.discard(b)

    return R, ops


def extended_gcd(a, b):
    """
    Compute x and y such that ax - by = 1
    (This function is completely classical)
    """
    # adapted from https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm

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


def check_mult_sizes(A, B, C):
    """
    bounds check for arguments to multiplication functions
    """

    # if it has a len, we consider it a register, otherwise an int
    try:
        len_a = len(A)
    except TypeError:
        # also, bit_length doesn't make sense for negative integers
        if A < 0:
            raise ValueError("allow_overflow required for negative integers") \
                from None
        len_a = A.bit_length()

    if len(C) < len_a+len(B):
        raise ValueError("register C not large enough to store result")


def get_cutoff(cutoff_in, default):
    """
    find and validate the cutoff for recursing in Karatsuba algorithms

    default values should be found using optimize_cutoffs.py
    """
    if cutoff_in is None:
        rtn = default
    else:
        rtn = cutoff_in

    if rtn < 4:
        # otherwise we end up infinitely recursing
        raise ValueError("cutoff must be >= 4")

    return rtn


def half_adder(A, B, Cout):
    """
    inputs:  A  B    Cout
    outputs: A  A+B  Cout+(A+B)/2
    """
    yield cirq.TOFFOLI(A, B, Cout)
    yield cirq.CNOT(A, B)


def full_adder(A, B, Cin, Cout):
    """
    inputs:   A  B         Cin  Cout
    outputs:  A  A+B+Cin%2 Cin  Cout+(A+B+C/2)

    the qubits A and Cin are measured in the H basis and zeroed
    after the adder is performed.
    """
    yield cirq.TOFFOLI(A, B, Cout)
    yield cirq.CNOT(A, B)
    yield cirq.TOFFOLI(B, Cin, Cout)
    yield cirq.CNOT(Cin, B)


def copy_register(A, B):
    """
    XOR the contents of A into B (thus copying them if B is empty)

    inputs:   A  B
    outputs:  A  A xor B
    """
    if len(A) > len(B):
        raise ValueError("register B must be at least as long as A")

    for a, b in zip(A, B):
        yield cirq.CNOT(a, b)
