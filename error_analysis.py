
from random import random, randrange, choice
import argparse
import cirq
from sympy.ntheory.generate import randprime
from sympy.ntheory.residue_ntheory import legendre_symbol
import progressbar

from digital_circuits import x2_mod_N, get_registers
from main import fast_flatten
from tof_sim import ToffoliSimulator, int_to_state, state_to_int
from ancilla import AncillaManager

# TODO: add phase errors

def error_circuit(circuit, error_rate):
    '''
    The circuit, but with bit flip errors added randomly after each gate
    with probability error_rate
    '''
    for op in circuit:
        # randomly apply X gates to qubits involved with the gate
        if random() < error_rate:
            qubit = choice(op.qubits)
            yield cirq.X(qubit)

        yield op


def test_circuit(circuit, regs, error_rate, iters, p, q, R, factor):
    gates = list(fast_flatten(circuit))
    x_reg, y_reg = regs
    N = p*q

    wrong_fail = 0
    wrong_pass = 0

    widgets = [
        'Error simulation: ', progressbar.Percentage(),
        ' ', progressbar.Bar(),#marker=progressbar.RotatingMarker()),
        ' ', progressbar.ETA(),
    ]
    bar = progressbar.ProgressBar(widgets=widgets, max_value=iters).start()

    all_qubits = circuit.all_qubits()

    try:
        for i in range(iters):
            c = error_circuit(gates, error_rate)
            sim = ToffoliSimulator([c], qubits=all_qubits)

            x = 0
            # gcd(x, N) should be 1
            while x%p == 0 or x%q == 0:
                x = randrange(N//2)

            state = int_to_state(0, all_qubits)
            state.update(int_to_state(x, x_reg))
            sim.simulate(state, check=False)
            result = state_to_int(state, y_reg[R.bit_length()-1:])

            factor_N = factor**2 * N
            y = result*R % factor_N
            correct = y == (factor*x)**2 % factor_N
            fail = not valid_y(y, p, q, factor)
            if correct:
                if fail:
                    print('warning: correct answer failed check!')
            else:
                if fail:
                    wrong_fail += 1
                else:
                    wrong_pass += 1

            bar.update(i+1)

    except KeyboardInterrupt as e:
        iters = i-1

    bar.finish()

    right = iters - wrong_fail - wrong_pass
    return iters, wrong_fail, wrong_pass, right


def valid_y(y, p, q, factor):
    # it is a quadratic residue mod N if it is one mod both p and q
    if not all(legendre_symbol(y//factor**2, f) == 1 for f in (p, q)):
        return False

    # it must be divisible by factor**2
    if not y % factor**2 == 0:
        return False

    return True


def choose_primes(n):
    p_len = n//2
    q_len = n-p_len

    p_max = 2**p_len
    q_max = 2**q_len

    p = randprime(3*p_max//4, p_max)
    q = randprime(3*q_max//4, q_max)

    return (p, q)


def get_error_rate(circuit, circuit_fid):
    '''
    Compute a per-gate error rate required to obtain the desired circuit fidelity
    '''
    ngates = sum(1 for _ in circuit.all_operations())
    return 1 - circuit_fid**(1/ngates)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Estimate success rates with postselection'
    )

    parser.add_argument('-n', type=int, required=True, help='number of bits of N')
    parser.add_argument('-f', type=float, required=True, help='circuit fidelity')
    parser.add_argument('--iters', type=int, default=1000, help='number of iterations')
    parser.add_argument('-m', type=int, default=0, help='number of factors of 3 to include')

    return parser.parse_args()


def main():
    args = parse_args()

    p, q = choose_primes(args.n)
    N = p*q

    x_reg, y_reg = get_registers(args.n, 1)
    ancillas = AncillaManager()
    R, circ_gen = x2_mod_N(N, x_reg, y_reg, ancillas, 'karatsuba')
    circuit = cirq.Circuit(circ_gen, strategy=cirq.InsertStrategy.NEW)
    error_rate = get_error_rate(circuit, args.f)

    orig_circuit = circuit

    # we want the error rate of the original circuit, but to run on the multiplied one
    factor = 3**args.m

    if factor > 1:
        x_reg, y_reg = get_registers(args.n, factor)
        ancillas = AncillaManager()
        R, circ_gen = x2_mod_N(N, x_reg, y_reg, ancillas, 'karatsuba', threes=args.m)
        circuit = cirq.Circuit(circ_gen, strategy=cirq.InsertStrategy.NEW)

    results = test_circuit(circuit, (x_reg, y_reg), error_rate, args.iters, p, q, R, factor)
    iters, wrong_fail, wrong_pass, right = results

    count_width = len(str(iters))
    print(f'iterations:        {iters}')
    print(f'correct:           {str(right).rjust(count_width)} ({right/iters:0.04f})')
    print(f'wrong, rejected:   {str(wrong_fail).rjust(count_width)} ({wrong_fail/iters:0.04f})')
    print(f'wrong, accepted:   {str(wrong_pass).rjust(count_width)} ({wrong_pass/iters:0.04f})')
    print()
    print(f'Succ. without postselection: {right/iters:0.02f}')
    print(f'Succ. with postselection:    {right/(right+wrong_pass):0.02f}')
    print()

    fact_iters = iters/(right+wrong_pass)
    fact_circ = len(circuit)/len(orig_circuit)
    print(f'Runtime factor (iters):      {fact_iters:.02f}')
    print(f'Runtime factor (circuit):    {fact_circ:.02f}')
    print(f'Runtime factor (total):      {fact_iters*fact_circ:.02f}')

if __name__ == '__main__':
    main()
