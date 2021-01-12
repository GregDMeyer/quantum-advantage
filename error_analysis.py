
from random import random, randrange, choice
import argparse
import cirq
from sympy.ntheory.generate import randprime
from sympy.ntheory.residue_ntheory import legendre_symbol

from digital_circuits import x2_mod_N
from main import fast_flatten
from tof_sim import ToffoliSimulator, int_to_state, state_to_int
from ancilla import AncillaManager


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


def test_circuit(circuit, regs, error_rate, iters, p, q, R):
    gates = list(fast_flatten(circuit))
    x_reg, y_reg = regs
    N = p*q

    wrong_fail = 0
    wrong_pass = 0
    for _ in range(iters):
        c = error_circuit(gates, error_rate)
        sim = ToffoliSimulator([c], qubits=circuit.all_qubits())

        x = 0
        # gcd(x, N) should be 1
        while x%p == 0 or x%q == 0:
            x = randrange(N//2)

        state = int_to_state(0, circuit.all_qubits())
        state.update(int_to_state(x, x_reg))
        sim.simulate(state, check=False)
        result = state_to_int(state, y_reg[len(x_reg):])

        y = result*R % N
        correct = y == x**2 % N
        fail = not valid_y(y, p, q)
        if correct:
            if fail:
                print('warning: correct answer failed check!')
        else:
            if fail:
                wrong_fail += 1
            else:
                wrong_pass += 1

    right = iters - wrong_fail - wrong_pass
    return wrong_fail, wrong_pass, right


def valid_y(y, p, q):
    # it is a quadratic residue mod N if it is one mod both p and q
    if not legendre_symbol(y, p) == 1:
        return False
    if not legendre_symbol(y, q) == 1:
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

    return parser.parse_args()


def main():
    args = parse_args()

    p, q = choose_primes(args.n)
    N = p*q

    x_reg = cirq.NamedQubit.range(args.n, prefix="x")
    y_reg = cirq.NamedQubit.range(2*args.n+1, prefix="y")

    ancillas = AncillaManager()
    R, circ_gen = x2_mod_N(N, x_reg, y_reg, ancillas, 'karatsuba')

    circuit = cirq.Circuit(circ_gen, strategy=cirq.InsertStrategy.NEW)

    error_rate = get_error_rate(circuit, args.f)

    results = test_circuit(circuit, (x_reg, y_reg), error_rate, args.iters, p, q, R)
    wrong_fail, wrong_pass, right = results

    count_width = len(str(args.iters))
    print(f'iterations:        {args.iters}')
    print(f'correct:           {str(right).rjust(count_width)} ({right/args.iters:0.04f})')
    print(f'wrong, rejected:   {str(wrong_fail).rjust(count_width)} ({wrong_fail/args.iters:0.04f})')
    print(f'wrong, accepted:   {str(wrong_pass).rjust(count_width)} ({wrong_pass/args.iters:0.04f})')
    print()
    print(f'Succ. without postselection: {right/args.iters:0.02f}')
    print(f'Succ. with postselection:    {right/(right+wrong_pass):0.02f}')
    print(f'Runtime factor:              {args.iters/(right+wrong_pass):.02f}')

if __name__ == '__main__':
    main()
