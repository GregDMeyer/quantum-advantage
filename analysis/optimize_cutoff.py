"""
Karatsuba multiplication is a recursive algorithm, and it is
generally optimal to switch to a different multiplication method
at some level of recursion. This script determines the optimal cutoff
for recursion by generating circuits using various cutoff values
and comparing circuit sizes to the schoolbook multiplication
algorithm.

Command-line options and usage can be viewed by passing the `-h` flag
(or inspecting the parse_args function below).

(c) Gregory D. Kahanamoku-Meyer 2021
"""

import random
from argparse import ArgumentParser
from itertools import chain

import cirq

# add circuit files to path
import sys
from os.path import abspath, dirname, join
sys.path.append(join(dirname(abspath(__file__)), '../circuits'))

from count_gates import describe, fast_flatten
from ancilla import AncillaManager
from digital_circuits import karatsuba_mult, schoolbook_mult
from digital_circuits import karatsuba_square, schoolbook_square
from digital_circuits import karatsuba_classical_mult, schoolbook_classical_mult


def parse_args():
    p = ArgumentParser(
        description='Optimize the cutoff for Karatsuba recursion.'
    )
    p.add_argument('op', choices=['mult', 'square', 'classical'],
                   default='mult', help='which operation to optimize')
    return p.parse_args()


def main():
    cutoff = 3

    args = parse_args()

    if args.op == 'mult':
        karatsuba, schoolbook = get_mult_fns()
    elif args.op == 'square':
        karatsuba, schoolbook = get_square_fns()
    elif args.op == 'classical':
        karatsuba, schoolbook = get_classical_mult_fns()

    k_count = 1
    s_count = 0

    print("cutoff   k_counts   s_counts")

    while k_count >= s_count:
        cutoff += 1
        k_count, s_count = compare_counts(cutoff, karatsuba, schoolbook)
        print(f"{str(cutoff).rjust(6)} {str(k_count).rjust(10)} "
              f"{str(s_count).rjust(10)}")


def compare_counts(cutoff, karatsuba, schoolbook):
    ck = karatsuba(cutoff)
    cs = schoolbook(cutoff)

    k_gates, *_ = describe(fast_flatten(ck))
    s_gates, *_ = describe(fast_flatten(cs))

    return k_gates, s_gates


def get_mult_fns():
    def get_registers(cutoff):
        n = cutoff+1
        a_reg = cirq.NamedQubit.range(n, prefix="a")
        b_reg = cirq.NamedQubit.range(n, prefix="b")
        c_reg = cirq.NamedQubit.range(2*n, prefix="c")
        ancillas = AncillaManager()
        return a_reg, b_reg, c_reg, ancillas

    def karatsuba(cutoff):
        return karatsuba_mult(*get_registers(cutoff), cutoff=cutoff)

    def schoolbook(cutoff):
        return schoolbook_mult(*get_registers(cutoff))

    return karatsuba, schoolbook


def get_square_fns():
    def get_registers(cutoff):
        n = cutoff+1
        a_reg = cirq.NamedQubit.range(n, prefix="a")
        c_reg = cirq.NamedQubit.range(2*n, prefix="c")
        ancillas = AncillaManager()
        return a_reg, c_reg, ancillas

    def karatsuba(cutoff):
        return karatsuba_square(*get_registers(cutoff), cutoff=cutoff)

    def schoolbook(cutoff):
        return schoolbook_square(*get_registers(cutoff))

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

    def karatsuba(cutoff):
        rtn = []
        for it in range(iters):
            rtn.append(
                karatsuba_classical_mult(
                    *get_args(cutoff, it),
                    cutoff=cutoff
                )
            )
        return chain(*rtn)

    def schoolbook(cutoff):
        rtn = []
        for it in range(iters):
            rtn.append(
                karatsuba_classical_mult(
                    *get_args(cutoff, it),
                )
            )
        return chain(*rtn)

    return karatsuba, schoolbook


if __name__ == '__main__':
    main()
