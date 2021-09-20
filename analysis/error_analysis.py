"""
This script numerically analyzes the effectiveness of the post-selection
scheme described in the paper associated with this software.

Command-line options and usage can be viewed by passing the `-h` flag
(or inspecting the parse_args function below).

(c) Gregory D. Kahanamoku-Meyer 2021
"""

from random import random, randrange, choice, seed
import argparse
from copy import deepcopy
from math import cos, pi, atan
from sys import stdout, stderr

import cirq
from sympy.ntheory.generate import randprime
from sympy.ntheory.residue_ntheory import legendre_symbol
import progressbar

# add circuit files to path
import sys
from os.path import abspath, dirname, join
sys.path.append(join(dirname(abspath(__file__)), '../circuits'))

from count_gates import fast_flatten
from digital_circuits import x2_mod_N, get_registers, extended_gcd
from tof_sim import ToffoliSimulator, int_to_state, state_to_int
from ancilla import AncillaManager


def error_circuit(circuit, error_rate, rand_seed=None):
    '''
    The circuit, but with bit flip errors added randomly after each gate
    with probability error_rate
    '''
    seed(a=rand_seed)
    err_gates = [cirq.X, cirq.Y, cirq.Z]
    for op in circuit:
        # randomly apply X gates to qubits involved with the gate
        if random() < error_rate:
            qubit = choice(op.qubits)
            gate = choice(err_gates)
            yield gate(qubit)

        yield op


def test_circuit(circuit, regs, error_rate, iters, p, q, R, factor, fidelity):
    gates = list(fast_flatten(circuit))
    x_reg, y_reg = regs
    N = p*q

    widgets = [
        'Error simulation: ', progressbar.Percentage(),
        ' ', progressbar.Bar(),#marker=progressbar.RotatingMarker()),
        ' ', progressbar.ETA(),
    ]
    bar = progressbar.ProgressBar(widgets=widgets, max_value=iters, fd=stdout).start()

    all_qubits = circuit.all_qubits()

    # so that we get a different run each time
    # start off with the system time seed
    seed_delta = randrange(2**32)

    post_data = {
        'total' : 0,
        'good_x': 0,
        'good_mz': 0,
        'good_mx': 0,
    }

    all_data = deepcopy(post_data)

    try:
        for i in range(iters):
            x0 = 0
            # gcd(x, N) should be 1
            while x0%p == 0 or x0%q == 0:
                x0 = randrange(N//2)

            x1 = compute_collision(x0, p, q)
            assert x0 != x1 and x0**2 % N == x1**2 % N

            combined_x = 0
            combined_y = 0
            combined_phase = 1
            n_pass = 0
            n_correct = 0

            for x in (x0, x1):
                c = error_circuit(gates, error_rate, rand_seed=i^seed_delta)
                sim = ToffoliSimulator([c], qubits=all_qubits)

                state = int_to_state(0, all_qubits)
                state.update(int_to_state(x, x_reg))
                phase = sim.simulate(state)
                x_result = state_to_int(state, x_reg)
                y_result = state_to_int(state, y_reg[R.bit_length()-1:])

                combined_x ^= x_result
                combined_y ^= y_result
                combined_phase *= phase

                factor_N = factor**2 * N
                y = y_result*R % factor_N
                fail = not valid_y(y, p, q, factor)

                correct = y == (factor*x)**2 % factor_N
                correct = correct and x_result == factor*x

                data = [all_data]
                if not fail:
                    data += [post_data]
                    n_pass += 1

                if correct:
                    n_correct += 1
                    if fail:
                        print('Warning: rejected correct value!')

                for d in data:
                    d['total'] += 1
                    if correct:
                        d['good_x'] += 1     # every x case
                        d['good_mz'] += 1   # Z polarized
                    else:
                        d['good_mz'] += 0.5  # Z polarized and lucky

            # finally, handle the case of m branch, X polarized

            # first case: both y and x were correct!
            # only possible error is a phase error
            if n_correct == 2:
                # either neither flipped, or they both did!
                if combined_phase == 1:
                    all_data['good_mx'] += 2
                    if n_pass == 2:
                        post_data['good_mx'] += 2

            # otherwise, at least one y or x was wrong. in this case we will
            # never get a correct superposition, and the single-qubit state
            # will be basically unrelated to what it is supposed to be.
            # therefore we give it 50% chance that you just happen to get it
            # right.
            # also putting a dent into the probability is the fact that if you
            # don't provide a valid y you auto-fail.
            else:
                all_data['good_mx'] += n_pass*0.5
                post_data['good_mx'] += n_pass*0.5

            bar.update(i+1)

    except KeyboardInterrupt:
        pass

    bar.finish()

    # add measurement factor to the m ones
    for d in (all_data, post_data):
        pz = d['good_mz']/d['total']
        px = d['good_mx']/d['total']

        if d is post_data:
            print()
            print(f'Good mz: {pz:0.5f}')
            print(f'Good mx: {px:0.5f}')

            theta = compute_optimal_theta(fidelity, factor)

            print(f'Optimal theta: {theta/pi:0.5f}*pi')

        else:
            theta = pi/8

        cz = cos(theta)**2
        cx = cos(theta-pi/4)**2

        d['good_m'] = cz*pz + (1-cz)*(1-pz) + cx*px + (1-cx)*(1-px)
        d['good_m'] *= d['total']/2

    return all_data, post_data


def compute_optimal_theta(f, factor):
    """
    compute the optimal measurement angle, given fidelity f and factor
    applied for postselection
    """
    ngood = f
    f1 = f**(1/3)
    nphase = (1-f1)*f1**2
    nerr = 0.25 * (1/(factor**2)) * (1-f1**2)

    ngood, nphase, nerr = [x/(ngood+nphase+nerr) for x in (ngood, nphase, nerr)]

    px = ngood + nphase
    pmz = ngood + nphase + 0.5*nerr
    pmx = ngood + 0.5*(nerr + nphase)

    return 0.5*atan((2*pmx-1)/(2*pmz-1))


def valid_y(y, p, q, factor):
    # it is a quadratic residue mod N if it is one mod both p and q
    if not all(legendre_symbol(y//factor**2, f) == 1 for f in (p, q)):
        return False

    # it must be divisible by factor**2
    if not y % factor**2 == 0:
        return False

    return True


def compute_collision(x0, p, q):
    mp = pow(x0, (p+1)//2, p)
    mq = pow(x0, (q+1)//2, q)
    yp, yq = extended_gcd(p, q)
    yq = p-yq  # the implementation in digital_circuits inverts it

    N = p*q
    xp0 = (yp*p*mq + yq*q*mp) % N
    if xp0 > N//2:
        xp0 = N-xp0

    if xp0 != x0:
        return xp0

    # otherwise, it was the other one
    xp1 = (yp*p*mq - yq*q*mp) % N
    if xp1 > N//2:
        xp1 = N-xp1

    return xp1


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

    parser.add_argument('-p', type=int, required=True, help='prime p')
    parser.add_argument('-q', type=int, required=True, help='prime q')
    parser.add_argument('-f', type=float, required=True, help='circuit fidelity')
    parser.add_argument('--iters', type=int, default=1000, help='number of iterations')
    parser.add_argument('-m', type=int, default=0, help='number of factors of 3 to include')

    return parser.parse_args()


def trunc_int(x):
    '''
    Return a string representing an integer, with ellipsis in the middle if
    the integer is huge.
    '''
    s = str(x)
    if len(s) < 12:
        return s

    return s[:4] + '...' + s[-4:]

def main():
    args = parse_args()

    N = args.p*args.q
    n = N.bit_length()

    print("Parameters:")
    print(f"N = {trunc_int(N)} = {trunc_int(args.p)}*{trunc_int(args.q)}")
    print(f"n = {n}", file=stderr)
    print(f"f = {args.f}", file=stderr)
    print(f"m = {args.m}", file=stderr)

    x_reg, y_reg = get_registers(n, 1)
    ancillas = AncillaManager()
    R, circ_gen = x2_mod_N(N, x_reg, y_reg, ancillas, 'karatsuba')
    circuit = cirq.Circuit(circ_gen, strategy=cirq.InsertStrategy.NEW)
    error_rate = get_error_rate(circuit, args.f)

    orig_circuit = circuit

    # we want the error rate of the original circuit, but to run on the multiplied one
    if args.m > 0:
        x_reg, y_reg = get_registers(n, 3**args.m)
        ancillas = AncillaManager()
        R, circ_gen = x2_mod_N(N, x_reg, y_reg, ancillas, 'karatsuba', threes=args.m)
        circuit = cirq.Circuit(circ_gen, strategy=cirq.InsertStrategy.NEW)

    results = test_circuit(circuit, (x_reg, y_reg), error_rate, args.iters, args.p, args.q, R, 3**args.m, args.f)
    all_data, post_data = results

    print()
    print('Results, no post-selection:')
    print_results(all_data)
    print()

    print('Results, with post-selection: ')
    print(f'samples: {post_data["total"]}', file=stderr)
    print_results(post_data, use_stderr=True)
    print()

    fact_iters = all_data['total'] / post_data['total']
    fact_circ = len(circuit) / len(orig_circuit)
    print(f'Runtime factor (iters):   {fact_iters:.02f}')
    print(f'Runtime factor (circuit): {fact_circ:.02f}')
    print('Runtime factor (total):   ', end='', flush=True)
    print(f'{fact_iters*fact_circ:.02f}', file=stderr)


def print_results(data, use_stderr=False):
    px = data['good_x'] / data['total']
    pm = data['good_m'] / data['total']
    result = px + 4*pm - 4

    print(f' px = {px:0.4f}')
    print(f' pm = {pm:0.4f}')
    print(' px + 4*pm - 4 = ', end='', flush=True)

    result_str = f'{result:0.4f}'
    if use_stderr:
        print(result_str, file=stderr)
    else:
        print(result_str)

    if result > 0:
        print('✓ exceeded the classical bound!')
    else:
        print('✗ did not exceed classical bound.')


if __name__ == '__main__':
    main()
