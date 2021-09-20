"""
This script generates circuits for computing x^2 mod N, and counts the number
of gates in the circuit, as well as gate depth and other useful quantities.

Command-line options and usage can be viewed by passing the `-h` flag
(or inspecting the parse_args function below).

(c) Gregory D. Kahanamoku-Meyer 2021
"""

from argparse import ArgumentParser
from random import randint
from sys import stderr

import cirq
from matplotlib import pyplot as plt

# add circuit files to path
import sys
from os.path import abspath, dirname, join
sys.path.append(join(dirname(abspath(__file__)), '../circuits'))

from digital_circuits import x2_mod_N
from phase_circuits import x2_mod_N_phase
from ancilla import AncillaManager


def parse_args():
    p = ArgumentParser(description='Count gates for x^2 mod N circuits.')
    p.add_argument('ns', type=lambda lst: [int(x) for x in lst.split(',')],
                   help='a comma separated list of problem sizes n=log2(N)')
    p.add_argument('-p', choices=['gates', 'depth', 'off'], default='gates',
                   help='which parameter to plot, or "off"')
    p.add_argument('-d', action='store_true', help='compute depth')
    p.add_argument('--no-decompose', action='store_true',
                   help='do not decompose Toffoli into Clifford+T')

    all_methods = ['karatsuba', 'schoolbook', 'narrow', 'fast']
    p.add_argument('--methods', type=lambda x: x.split(','), default=all_methods)

    args = p.parse_args()

    for m in args.methods:
        if m not in all_methods:
            raise ValueError(f"unrecognized method '{m}'")

    if args.p == 'depth' and not args.d:
        raise ValueError("Pass the -d flag to plot circuit depth.")

    return args


def describe(c, decompose=True):

    if not isinstance(c, cirq.Circuit):
        # add one layer of iteration---we are pretending
        # it is one giant moment
        c = (c,)

    depth = 0
    tot_gates = 0
    t_gates = 0
    for m in c:
        has_tof = False
        for op in m:
            if isinstance(op, int):
                tot_gates += op
                continue

            if decompose and op.gate is cirq.TOFFOLI:
                tot_gates += 15
                t_gates += 7
                has_tof = True
            else:
                tot_gates += 1

        if decompose and has_tof:
            depth += 12
        else:
            depth += 1

    return (tot_gates, t_gates, depth)


def plot(ns, ydata, yname):
    for impl, vals in ydata.items():
        plt.plot(ns, vals, label=impl)

    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel('n')
    plt.ylabel(yname)

    plt.legend()

    plt.tight_layout()
    plt.show()


def fast_flatten(root):
    '''
    flatten an op-tree
    '''
    if isinstance(root, (int, cirq.Operation)):
        yield root
    else:
        for subtree in root:
            yield from fast_flatten(subtree)


def shim_gates():
    '''
    This function causes circuit building to directly yield gate counts
    instead of actually constructing gate objects. This yields much better
    performance when only the gate counts are desired (the circuit does not
    actually need to be constructed).
    '''
    gate_gens = [
        'ZPowGate',
        'CZPowGate',
        'CCZPowGate'
    ]

    def fake_gate_gen(*args, **kwargs):
        return lambda *qubits: 1

    for gate_name in gate_gens:
        setattr(cirq, gate_name, fake_gate_gen)


def main():
    args = parse_args()

    digital_impls = [
        'schoolbook',
        'karatsuba'
    ]

    phase_impls = [
        'narrow',
        'fast'
    ]

    header = ['impl', 'n', 'gates', 't_gates', 'qubits', 'measurements']
    if args.d:
        header.append('depth')
    print(','.join(header))

    ydata = {impl: [] for impl in digital_impls+phase_impls}
    for n in args.ns:
        # for a fair comparison across implementations
        N = randint(0, 2**n-1) | 1 | (1<<(n-1))

        for impl in args.methods:

            x_reg = cirq.NamedQubit.range(n, prefix="x")
            ancillas = AncillaManager()

            if impl in digital_impls:
                y_reg = cirq.NamedQubit.range(2*n+1, prefix="y")
                R, circ_gen = x2_mod_N(N, x_reg, y_reg, ancillas, impl)

            elif impl in phase_impls:
                y_reg = cirq.NamedQubit.range(n+3, prefix="y")
                circ_gen = x2_mod_N_phase(N, x_reg, y_reg, ancillas, impl)

            else:
                raise ValueError("unknown implementation")

            # only actually build the circuit if we need the machinery
            # to parallelize gates. otherwise we can just iterate through the
            # generator.
            if args.d:
                c = cirq.Circuit(circ_gen, strategy=cirq.InsertStrategy.EARLIEST)
            else:
                shim_gates()
                c = fast_flatten(circ_gen)

            if not ancillas.all_discarded():
                print(f"qubit memory leak! {ancillas.n_active} "
                      "ancillas not discarded", file=stderr)

            gates, t_gates, depth = describe(c, not args.no_decompose)
            meas = ancillas.count_measurements()
            gates += meas

            qubits = len(x_reg) + len(y_reg) + ancillas.max_ancilla_usage()
            print(f"{impl}, {n}, {gates}, {t_gates}, {qubits}, {meas}"
                  + (f", {depth}" if args.d else ""))

            if args.p != 'off':
                ydata[impl].append({'gates': gates,
                                    'depth': depth}[args.p])

    if args.p != 'off':
        plot(args.ns, ydata, args.p)


if __name__ == '__main__':
    main()
