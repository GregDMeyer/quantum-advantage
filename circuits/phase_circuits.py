"""
This file contains functions for generating quantum circuits
for x^2 mod N that compute the function in the phase, and then
perform an IQFT to transfer the result to a quantum register.
These circuits are diagrammed in the manuscript associated with
this software.

The functions here return generators of the circuit's gates (an "op tree"
in the terminology used by Cirq). It can be converted to a Cirq circuit by
being passed to the cirq.Circuit(op_tree) constructor.

(c) Gregory D. Kahanamoku-Meyer 2021
"""

import itertools
import cirq


def x2_mod_N_phase(N, x, y, ancillas, circuit_type):
    '''
    Generate the full circuit
    '''
    for q in itertools.chain(x, y):
        yield cirq.H(q)

    if circuit_type == 'fast':
        yield x2modN_fast(x, y, ancillas, N)
        yield iqft(y)

    elif circuit_type == 'narrow':
        for k, yk in enumerate(y):
            yield x2modN_narrow_single(x, [yk], ancillas, N, len(y)-k-1)
            yield cirq.H(yk)
            for m, ym in list(enumerate(y))[k+1:]:
                phase = -(1/2**(m-k))
                yield CZpow(phase, ym, yk)

    else:
        raise ValueError('unknown circuit type')


def x2modN_narrow_single(x, y, ancillas, N, bit_idx):
    for i, xi in enumerate(x):
        for j, xj in list(enumerate(x))[i:]:
            if i == j:
                controls = [xi]
                fctr = 1
            else:
                controls = [xi, xj]
                fctr = 2

            phase = fctr*2*2**(i+j+bit_idx)/N
            yield phase_add(phase, y, controls)


def x2modN_fast(x, y, ancillas, N):
    n = len(x)
    counter = ancillas.new_register(n.bit_length())
    for q in counter:
        yield cirq.H(q)

    for l in range(2*n-1):
        qubit_pairs = [(i, l-i) for i in range(max(0, l-n+1), l//2 + 1)]

        yield count(x, counter, qubit_pairs, sign=+1)
        yield iqft(counter)

        for k, ck in enumerate(counter):
            phase = 2*2**(k+l) / N
            yield phase_add(phase, y, [ck])

        # uncompute
        yield qft(counter)
        yield count(x, counter, qubit_pairs, sign=-1)

    for q in counter:
        yield cirq.H(q)

    ancillas.discard(counter)


def count(x, counter, pairs, sign):
    for i,j in pairs:
        if i == j:
            controls = [x[i]]
            fctr = 1
        else:
            controls = [x[i], x[j]]
            fctr = 2

        phase = sign*fctr*(2**(1-len(counter)))
        yield phase_add(phase, counter, controls)


def phase_add(in_phase, target, controls):
    for k, yk in enumerate(target):
        phase = 2**(len(target) - k - 1) * in_phase % 2
        yield CZpow(phase, yk, *controls)


def iqft(reg):
    for k, yk in enumerate(reg):
        yield cirq.H(yk)
        for m, ym in list(enumerate(reg))[k+1:]:
            phase = -(1/2**(m-k))
            yield CZpow(phase, ym, yk)


def qft(reg):
    for k, yk in list(enumerate(reg))[::-1]:
        for m, ym in list(enumerate(reg))[k+1:]:
            phase = 1/2**(m-k)
            yield CZpow(phase, ym, yk)
        yield cirq.H(yk)


def CZpow(phase, *qubits):
    nqubits = len(qubits)

    # we put most common case first
    if nqubits == 3:
        base_gate = cirq.CCZPowGate(exponent=phase)
    elif nqubits == 2:
        base_gate = cirq.CZPowGate(exponent=phase)
    elif nqubits == 1:
        base_gate = cirq.ZPowGate(exponent=phase)
    else:
        raise ValueError('too many control qubits')

    return base_gate(*qubits)
