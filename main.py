
from argparse import ArgumentParser
from random import randint
import cirq
from digital_circuits import x2_mod_N
from phase_circuits import x2_mod_N_phase
from ancilla import AncillaManager

from matplotlib import pyplot as plt

def parse_args():
    p = ArgumentParser(description='Count gates for x^2 mod N circuits.')
    p.add_argument('ns', type=lambda lst: [int(x) for x in lst.split(',')],
                   help='a comma separated list of problem sizes to compute')
    p.add_argument('-p', choices=['gates', 'depth', 'off'], default='gates',
                   help='which parameter to plot, or "off"')
    p.add_argument('--method', choices=['schoolbook', 'karatsuba', 'narrow', 'fast'])

    args = p.parse_args()
    return args

def describe(c):

    depth = 0
    tot_gates = 0
    for m in c:
        has_tof = False
        for op in m:
            if op.gate is cirq.TOFFOLI:
                tot_gates += 15
                has_tof = True
            else:
                tot_gates += 1

        if has_tof:
            depth += 12
        else:
            depth += 1

    return (tot_gates, depth)

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

    header = ['impl', 'n', 'gates', 'depth']
    print(','.join(header))

    ydata = {impl:[] for impl in digital_impls+phase_impls}
    for n in args.ns:
        # for a fair comparison across implementations
        N = randint(0, 2**n-1) | 1 | (1<<(n-1))
        for impl in digital_impls:
            if args.method is not None and impl != args.method:
                continue

            c = cirq.Circuit()
            x_reg = cirq.NamedQubit.range(n, prefix="x")
            y_reg = cirq.NamedQubit.range(2*n+1, prefix="y")
            ancillas = AncillaManager()

            x2_mod_N(c, N, x_reg, y_reg, ancillas, impl)

            gates, depth = describe(c)
            print(f"{impl}, {n}, {gates}, {depth}")

            if args.p != 'off':
                ydata[impl].append({'gates':gates,
                                    'depth':depth}[args.p])

        for impl in phase_impls:
            if args.method is not None and impl != args.method:
                continue

            x_reg = cirq.NamedQubit.range(n, prefix="x")
            y_reg = cirq.NamedQubit.range(n+3, prefix="y")
            ancillas = AncillaManager()

            c = x2_mod_N_phase(N, x_reg, y_reg, ancillas, impl)

            gates, depth = describe(c)
            print(f"{impl}, {n}, {gates}, {depth}")

            if args.p != 'off':
                ydata[impl].append({'gates':gates,
                                    'depth':depth}[args.p])

    if args.p != 'off':
        plot(args.ns, ydata, args.p)

if __name__ == '__main__':
    main()
