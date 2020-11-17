
import itertools
import cirq

def x2modN_narrow_single(circ, x, y, ancillas, N, bit_idx):
    for i, xi in enumerate(x):
        for j, xj in list(enumerate(x))[i:]:
            if i == j:
                controls = [xi]
                fctr = 1
            else:
                controls = [xi, xj]
                fctr = 2

            phase = fctr*2*2**(i+j+bit_idx)/N
            phase_add(circ, phase, y, controls)

def x2modN_fast(circ, x, y, ancillas, N):
    n = len(x)
    counter = ancillas.new_register(n.bit_length())
    for q in counter:
        circ.append(cirq.H(q))

    for m in range(2*n-1):
        qubit_pairs = [(i,m-i) for i in range(max(0, m-n+1), m//2 + 1)]

        count(circ, x, counter, qubit_pairs, sign=+1)
        iqft(circ, counter)

        for k, ck in enumerate(counter):
            phase = 2*2**(k+m) / N
            phase_add(circ, phase, y, [ck])

        # uncompute
        qft(circ, counter)
        count(circ, x, counter, qubit_pairs, sign=-1)

    for q in counter:
        circ.append(cirq.H(q))

    ancillas.discard(counter)

def count(circ, x, counter, pairs, sign):
    for i,j in pairs:
        if i == j:
            controls = [x[i]]
            fctr = 1
        else:
            controls = [x[i], x[j]]
            fctr = 2

        phase = sign*fctr*(2**(1-len(counter)))
        phase_add(circ, phase, counter, controls)

def phase_add(circ, in_phase, target, controls):
    for k, yk in enumerate(target):
        phase = 2**(len(target) - k - 1) * in_phase % 2
        circ.append((cirq.Z(yk)**phase).controlled_by(*controls))

# def phase_increment(circ, reg, controls):
#     for k, yk in enumerate(reg):
#         phase = 2**(-k)
#         circ.append((cirq.Z(yk)**phase).controlled_by(controls))

def iqft(circ, reg):
    for k, yk in enumerate(reg):
        circ.append(cirq.H(yk))
        for m, ym in list(enumerate(reg))[k+1:]:
            phase = -(1/2**(m-k))
            circ.append(cirq.CZ(ym, yk)**phase)

def qft(circ, reg):
    for k, yk in list(enumerate(reg))[::-1]:
        for m, ym in list(enumerate(reg))[k+1:]:
            phase = 1/2**(m-k)
            circ.append(cirq.CZ(ym, yk)**phase)
        circ.append(cirq.H(yk))

def x2_mod_N_phase(circ, N, x, y, ancillas, circuit_type):
    '''
    Generate the full circuit
    '''
    for q in itertools.chain(x, y):
        circ.append(cirq.H(q))

    if circuit_type == 'fast':
        x2modN_fast(circ, x, y, ancillas, N)
        iqft(circ, y)

    elif circuit_type == 'narrow':
        for k, yk in enumerate(y):
            x2modN_narrow_single(circ, x, [yk], ancillas, N, len(y)-k-1)
            circ.append(cirq.H(yk))
            for m, ym in list(enumerate(y))[k+1:]:
                phase = -(1/2**(m-k))
                circ.append(cirq.CZ(ym, yk)**phase)

    else:
        raise ValueError('unknown circuit type')

    return circ
