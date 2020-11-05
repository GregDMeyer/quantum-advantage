
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
    p.add_argument('-d', action='store_true', help='compute depth')

    all_methods = ['karatsuba', 'schoolbook', 'narrow', 'fast']
    p.add_argument('--methods', type=lambda x: x.split(','), default=all_methods)

    args = p.parse_args()

    for m in args.methods:
        if m not in all_methods:
            raise ValueError(f"unrecognized method '{m}'")

    return args

def describe(c):

    depth = 0
    tot_gates = 0
    t_gates = 0
    for m in c:
        has_tof = False
        for op in m:
            if op.gate is cirq.TOFFOLI:
                tot_gates += 15
                t_gates += 7
                has_tof = True
            else:
                tot_gates += 1

        if has_tof:
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

    header = ['impl', 'n', 'gates', 't_gates']
    if args.d:
        header.append('depth')
    print(','.join(header))

    ydata = {impl:[] for impl in digital_impls+phase_impls}
    for n in args.ns:
        # for a fair comparison across implementations
        N = randint(0, 2**n-1) | 1 | (1<<(n-1))

        for impl in args.methods:

            c = cirq.Circuit()

            # adjust the default insertstrategy
            # this is kind of dumb, but alas
            if args.d:
                c.append = lambda arg: cirq.Circuit.append(c, arg, cirq.InsertStrategy.EARLIEST)
            else:
                c.append = lambda arg: cirq.Circuit.append(c, arg, cirq.InsertStrategy.NEW)

            x_reg = cirq.NamedQubit.range(n, prefix="x")
            ancillas = AncillaManager()

            if impl in digital_impls:
                y_reg = cirq.NamedQubit.range(2*n+1, prefix="y")
                x2_mod_N(c, N, x_reg, y_reg, ancillas, impl)

            elif impl in phase_impls:
                y_reg = cirq.NamedQubit.range(n+3, prefix="y")
                x2_mod_N_phase(c, N, x_reg, y_reg, ancillas, impl)

            else:
                raise ValueError("unknown implementation")

            gates, t_gates, depth = describe(c)
            print(f"{impl}, {n}, {gates}, {t_gates}" + (f", {depth}" if args.d else ""))

            if args.p != 'off':
                ydata[impl].append({'gates':gates,
                                    'depth':depth}[args.p])

    if args.p != 'off':
        plot(args.ns, ydata, args.p)

if __name__ == '__main__':
    main()
