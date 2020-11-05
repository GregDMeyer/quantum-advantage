
import random
from sys import argv

import cirq
from main import describe
from ancilla import AncillaManager
from circuits import karatsuba_mult, schoolbook_mult
from circuits import karatsuba_square, schoolbook_square
from circuits import karatsuba_classical_mult, schoolbook_classical_mult

def compare_counts(cutoff, karatsuba, schoolbook):
    ck = cirq.Circuit()
    cs = cirq.Circuit()

    karatsuba(ck, cutoff)
    schoolbook(cs, cutoff)

    k_gates, _ = describe(ck)
    s_gates, _ = describe(cs)

    return k_gates, s_gates

def main():
    cutoff = 3

    if argv[1] == 'mult':
        karatsuba, schoolbook = get_mult_fns()
    elif argv[1] == 'square':
        karatsuba, schoolbook = get_square_fns()
    elif argv[1] == 'classical':
        karatsuba, schoolbook = get_classical_mult_fns()
    else:
        raise ValueError(f"unknown mult type '{argv[1]}'")

    k_count = 1
    s_count = 0

    print("cutoff   k_counts   s_counts")

    while k_count >= s_count:
        cutoff += 1
        k_count, s_count = compare_counts(cutoff, karatsuba, schoolbook)
        print(f"{str(cutoff).rjust(6)} {str(k_count).rjust(10)} {str(s_count).rjust(10)}")

def get_mult_fns():
    def get_registers(cutoff):
        n = cutoff+1
        a_reg = cirq.NamedQubit.range(n, prefix="a")
        b_reg = cirq.NamedQubit.range(n, prefix="b")
        c_reg = cirq.NamedQubit.range(2*n, prefix="c")
        ancillas = AncillaManager()
        return a_reg, b_reg, c_reg, ancillas

    def karatsuba(circ, cutoff):
        karatsuba_mult(circ, *get_registers(cutoff), cutoff=cutoff)

    def schoolbook(circ, cutoff):
        schoolbook_mult(circ, *get_registers(cutoff))

    return karatsuba, schoolbook

def get_square_fns():
    def get_registers(cutoff):
        n = cutoff+1
        a_reg = cirq.NamedQubit.range(n, prefix="a")
        c_reg = cirq.NamedQubit.range(2*n, prefix="c")
        ancillas = AncillaManager()
        return a_reg, c_reg, ancillas

    def karatsuba(circ, cutoff):
        karatsuba_square(circ, *get_registers(cutoff), cutoff=cutoff)

    def schoolbook(circ, cutoff):
        schoolbook_square(circ, *get_registers(cutoff))

    return karatsuba, schoolbook

def get_classical_mult_fns():

    # do a few iterations for averaging
    iters = 4

    def get_args(cutoff, it):
        n = cutoff+1

        random.seed(0xC0DEBEE5+it)
        a = random.randint(2**(n-1), 2**n-1)

        b_reg = cirq.NamedQubit.range(n, prefix="b")
        c_reg = cirq.NamedQubit.range(2*n, prefix="c")
        ancillas = AncillaManager()
        return a, b_reg, c_reg, ancillas

    def karatsuba(circ, cutoff):
        for it in range(iters):
            karatsuba_classical_mult(circ, *get_args(cutoff, it), cutoff=cutoff)

    def schoolbook(circ, cutoff):
        for it in range(iters):
            schoolbook_classical_mult(circ, *get_args(cutoff, it))

    return karatsuba, schoolbook

if __name__ == '__main__':
    main()
