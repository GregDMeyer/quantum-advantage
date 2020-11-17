
from argparse import ArgumentParser
import cirq
from phase_circuits import x2_mod_N_phase
from ancilla import AncillaManager

def parse_args():
    p = ArgumentParser(description='Check fidelity of phase circuits.')

    p.add_argument('N', type=int, help='integer modulus')
    p.add_argument('--extra-bits', type=int, default=2, help='extra bits to improve QFT accuracy')
    p.add_argument('-t', choices=['narrow', 'fast'], default='narrow')
    args = p.parse_args()

    return args

def parse_results(results, N, registers):
    rtn = {}
    m = results.state_vector()

    for i,v in enumerate(m):
        x, y, anc = idx_to_vals(i, registers)
        true_y = round(y/(2**len(registers[1])) * N)
        prob = (v*v.conj()).real

        # if prob != 0:
        #     print(x, y, anc, prob)

        # make sure we uncomputed everything we needed to
        if anc != 0:
            assert prob == 0

        k = (true_y, x)
        if k not in rtn:
            rtn[k] = 0
        rtn[k] += prob

    return rtn

def idx_to_vals(idx, registers):
    rtn = []
    for reg in registers:
        v = idx & ((1<<len(reg))-1)
        rtn.append(v)
        idx >>= len(reg)
    return rtn

def check_a(results, N):
    good = 0
    for (y, x), prob in results.items():
        if abs(y - (x**2 % N)) < 1:
            good += prob
    return good

def main():
    args = parse_args()
    N = args.N

    domain_max = (N+1)//2 - 1
    x = cirq.NamedQubit.range(domain_max.bit_length(), prefix='x')
    y = cirq.NamedQubit.range(N.bit_length()+args.extra_bits, prefix='y')
    ancillas = AncillaManager()

    simulator = cirq.Simulator()
    circuit = cirq.Circuit()
    x2_mod_N_phase(circuit, N, x, y, ancillas, args.t)

    registers = x, y, ancillas.all()
    result = simulator.simulate(circuit, qubit_order=[q for r in registers for q in r][::-1])
    parsed = parse_results(result, N, registers)

    good_frac = check_a(parsed, N)

    print(f'success prob.: {good_frac:0.3f}')

if __name__ == '__main__':
    main()
