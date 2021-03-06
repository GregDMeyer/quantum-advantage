"""
This script uses the Python "unittest" module to explicitly
check correctness of the "digital" quantum circuits for Karatsuba
and schoolbook multiplication.

By using the polynomial-time quantum simulation defined in
../circuits/tof_sim.py, correctness is checked for circuits
consisting of hundreds of qubits.

(c) Gregory D. Kahanamoku-Meyer 2021
"""

import unittest

from itertools import product, combinations_with_replacement
import random
import cirq

# add circuit files to path
import sys
from os.path import abspath, dirname, join
sys.path.append(join(dirname(abspath(__file__)), '../circuits'))

from digital_circuits import full_adder, half_adder, add_int
from digital_circuits import times_three, add_classical_int, lessthan_classical
from digital_circuits import schoolbook_square, karatsuba_square
from digital_circuits import schoolbook_mult, karatsuba_mult
from digital_circuits import schoolbook_classical_mult, karatsuba_classical_mult
from digital_circuits import extended_gcd, montgomery_reduce
from digital_circuits import x2_mod_N, get_registers
from tof_sim import ToffoliSimulator, int_to_state, state_to_int
from ancilla import AncillaManager

NEW = cirq.InsertStrategy.NEW

class TestX2modN(unittest.TestCase):

    ns = [4, 8, 16, 50]
    threes = [0, 1, 2, 5]
    iters = 8
    methods = ['karatsuba', 'schoolbook']

    def setUp(self):
        random.seed(0xF00DCAFE)

    def test_x2modN(self):
        for n in self.ns:
            for three in self.threes:
                factor = 3**three
                x_reg, y_reg = get_registers(n, factor)
                for _, method in product(range(self.iters), self.methods):

                    N = random.randint(0, 2**n-1)  # we will have N = pq, but this is fine for testing
                    N |= 1           # N must be odd
                    N |= 1 << (n-1)  # N must be of length n

                    with self.subTest(N=N, method=method, three=three):

                        ancillas = AncillaManager()
                        R, circ_gen = x2_mod_N(N, x_reg, y_reg, ancillas, method, threes=three)
                        circ = cirq.Circuit(circ_gen, strategy=NEW)
                        sim = ToffoliSimulator(circ)

                        for _ in range(self.iters):
                            x = random.randint(0, N-1)
                            with self.subTest(x=x):
                                state = int_to_state(x, x_reg)
                                state.update(int_to_state(0, y_reg))
                                state.update(ancillas.init_state())
                                sim.simulate(state)
                                rx = state_to_int(state, x_reg)
                                ry = state_to_int(state, y_reg[R.bit_length()-1:]) # result is top bits of y register

                                self.assertEqual(rx, factor*x)

                                self.assertLess(ry, factor**2 * N)
                                self.assertEqual((R*ry) % (factor**2), 0)
                                true_y = (R*ry // factor**2) % N
                                self.assertEqual(true_y, x**2 % N)


class TestSubCircuits(unittest.TestCase):

    def test_half_adder(self):
        qubits = cirq.NamedQubit.range(3, prefix="q")
        c = cirq.Circuit(half_adder(*qubits), strategy=NEW)
        sim = ToffoliSimulator(c)

        for a, b in product((0, 1), repeat=2):
            with self.subTest(a=a, b=b):
                state = {x.name: v for x, v in zip(qubits, (a, b, 0))}
                sim.simulate(state)
                ra, rsum, rcout = [state[q.name] for q in qubits]

                self.assertEqual(ra, a)
                self.assertEqual(rsum, (a+b)%2)
                self.assertEqual(rcout, (a+b)//2)

    def test_full_adder(self):
        qubits = cirq.NamedQubit.range(4, prefix="q")
        c = cirq.Circuit(full_adder(*qubits), strategy=NEW)
        sim = ToffoliSimulator(c)

        for a, b, cin in product((0, 1), repeat=3):
            with self.subTest(a=a, b=b, cin=cin):
                state = {x.name: v for x, v in zip(qubits, (a, b, cin, 0))}
                sim.simulate(state)
                ra, rsum, rcin, rcout = [state[q.name] for q in qubits]

                self.assertEqual(ra, a)
                self.assertEqual(rsum, (a+b+cin)%2)
                self.assertEqual(rcin, cin)
                self.assertEqual(rcout, (a+b+cin)//2)


class TestArithmetic(unittest.TestCase):

    addition_test_cases = [
        (5, 5, 0, 0),
        (5, 5, 0, 1),
        (5, 5, 1, 0),
        (5, 5, 1, 1),
        (5, 5, 3, 0),
        (5, 5, 18, 30),
        (5, 5, 27, 27),
        (5, 5, 31, 20),
        (5, 5, 1, 31),
        (5, 5, 30, 12),
        (5, 8, 0, 0),
        (5, 8, 0, 255),
        (5, 8, 1, 255),
        (5, 8, 27, 130),
        (5, 8, 27, 250),
    ]

    def test_int_addition(self):
        for n, m, a, b in self.addition_test_cases:
            with self.subTest(a=a, b=b):

                a_reg = cirq.NamedQubit.range(n, prefix="a")
                b_reg = cirq.NamedQubit.range(m, prefix="b")
                ancillas = AncillaManager()

                c = cirq.Circuit(add_int(a_reg, b_reg, ancillas), strategy=NEW)
                sim = ToffoliSimulator(c)

                state = int_to_state(a, a_reg)
                state.update(int_to_state(b, b_reg))
                state.update(ancillas.init_state())
                sim.simulate(state)
                ra = state_to_int(state, a_reg)
                rb = state_to_int(state, b_reg)

                self.assertEqual(ra, a)
                self.assertEqual(rb, (a+b)%(2**m))

    def test_classical_int_addition(self):
        for _, m, a, b in self.addition_test_cases:
            with self.subTest(a=a, b=b):

                b_reg = cirq.NamedQubit.range(m, prefix="b")
                ancillas = AncillaManager()

                c = cirq.Circuit(add_classical_int(a, b_reg, ancillas), strategy=NEW)
                sim = ToffoliSimulator(c)

                state = int_to_state(b, b_reg)
                state.update(ancillas.init_state())
                sim.simulate(state)
                rb = state_to_int(state, b_reg)

                self.assertEqual(rb, (a+b)%(2**m))

    def test_int_addition_exceptions(self):
        a_reg = cirq.NamedQubit.range(5, prefix="a")
        b_reg = cirq.NamedQubit.range(6, prefix="b")
        ancillas = AncillaManager()

        c = cirq.Circuit(add_int(a_reg, b_reg, ancillas), strategy=NEW)

        with self.assertRaises(ValueError):
            next(add_int(b_reg, a_reg, ancillas))

        with self.assertRaises(ValueError):
            next(add_classical_int(32, a_reg, ancillas))

    def test_lessthan_classical(self):
        for _, m, a, b in self.addition_test_cases:
            with self.subTest(a=a, b=b):

                b_reg = cirq.NamedQubit.range(m, prefix="b")
                ancillas = AncillaManager()
                result = ancillas.new()

                c = cirq.Circuit(lessthan_classical(b_reg, a, result, ancillas), strategy=NEW)
                sim = ToffoliSimulator(c)

                state = int_to_state(b, b_reg)
                state.update(ancillas.init_state())
                sim.simulate(state)

                rb = state_to_int(state, b_reg)
                rresult = state_to_int(state, [result])

                self.assertEqual(rb, b)
                self.assertEqual(rresult, b<a)

    def test_square(self):
        n = 5
        a_test_cases = range(2**n)
        b_test_cases = (0, 1, 31, 2**(2*n)-1)

        square_methods = [
            ("schoolbook", schoolbook_square),
            ("karatsuba", lambda *args: karatsuba_square(*args, cutoff=4))
        ]

        for name, square in square_methods:
            with self.subTest(square_method=name):

                a_reg = cirq.NamedQubit.range(n, prefix="a")
                b_reg = cirq.NamedQubit.range(2*n, prefix="b")
                ancillas = AncillaManager()

                with self.assertRaises(ValueError):
                    next(square(a_reg, a_reg, ancillas))

                c = cirq.Circuit(square(a_reg, b_reg, ancillas), strategy=NEW)
                sim = ToffoliSimulator(c)

                for a in a_test_cases:
                    for b in b_test_cases:
                        with self.subTest(a=a, b=b):
                            state = int_to_state(a, a_reg)
                            state.update(int_to_state(b, b_reg))
                            state.update(ancillas.init_state())
                            sim.simulate(state)
                            ra = state_to_int(state, a_reg)
                            rb = state_to_int(state, b_reg)

                            self.assertEqual(ra, a)
                            self.assertEqual(rb, (b+a**2)%(2**(2*n)))

    def test_mult(self):
        n = 5
        m = 5
        test_cases = product(range(2**n), range(2**m))

        mult_methods = [
            ("schoolbook", schoolbook_mult),
            ("karatsuba", lambda *args: karatsuba_mult(*args, cutoff=4))
        ]

        for name, mult in mult_methods:
            with self.subTest(mult_method=name):

                a_reg = cirq.NamedQubit.range(n, prefix="a")
                b_reg = cirq.NamedQubit.range(m, prefix="b")
                c_reg = cirq.NamedQubit.range(m+n, prefix="c")
                ancillas = AncillaManager()

                circ1 = cirq.Circuit(mult(a_reg, b_reg, c_reg, ancillas), strategy=NEW)
                sim1 = ToffoliSimulator(circ1)

                circ2 = cirq.Circuit(mult(b_reg, a_reg, c_reg, ancillas), strategy=NEW)
                sim2 = ToffoliSimulator(circ2)

                for a, b in test_cases:
                    for c in (0, int(0.8122297*2**(n+m)), 2**(n+m)-1):
                        for case, sim in ("n<m", sim1), ("n>m", sim2):
                            with self.subTest(a=a, b=b, c=c, case=case):
                                state = int_to_state(a, a_reg)
                                state.update(int_to_state(b, b_reg))
                                state.update(int_to_state(c, c_reg))
                                state.update(ancillas.init_state())

                                sim.simulate(state)

                                ra = state_to_int(state, a_reg)
                                rb = state_to_int(state, b_reg)
                                rc = state_to_int(state, c_reg)

                                self.assertEqual(ra, a)
                                self.assertEqual(rb, b)
                                self.assertEqual(rc, (c+a*b)%(2**(n+m)))

    def test_classical_mult(self):
        n = 5
        a_test_cases = b_test_cases = range(2**n)

        mult_methods = [
            ("schoolbook", schoolbook_classical_mult),
            ("karatsuba", lambda *args: karatsuba_classical_mult(*args, cutoff=4))
        ]

        for name, mult in mult_methods:
            with self.subTest(mult_method=name):
                for a in a_test_cases:
                    b_reg = cirq.NamedQubit.range(n, prefix="b")
                    c_reg = cirq.NamedQubit.range(2*n, prefix="c")
                    ancillas = AncillaManager()

                    circ = cirq.Circuit(mult(a, b_reg, c_reg, ancillas), strategy=NEW)
                    sim = ToffoliSimulator(circ)

                    for b,c in product(b_test_cases, (0, int(0.8122297*2**(2*n)), 2**(2*n)-1)):
                        with self.subTest(a=a, b=b, c=c):
                            state = int_to_state(b, b_reg)
                            state.update(int_to_state(c, c_reg))
                            state.update(ancillas.init_state())

                            sim.simulate(state)

                            rb = state_to_int(state, b_reg)
                            rc = state_to_int(state, c_reg)

                            self.assertEqual(rb, b)
                            self.assertEqual(rc, (c+a*b)%(2**(2*n)))

    def test_times_three(self):
        n = 5
        test_cases = range(2**n)

        a_reg = cirq.NamedQubit.range(n+2, prefix="a")
        ancillas = AncillaManager()

        c = cirq.Circuit(times_three(a_reg, ancillas), strategy=NEW)
        sim = ToffoliSimulator(c)

        for a in test_cases:
            with self.subTest(a=a):
                state = int_to_state(a, a_reg)
                state.update(ancillas.init_state())
                sim.simulate(state)
                ra = state_to_int(state, a_reg)

                self.assertEqual(ra, 3*a)


class TestArithmeticLarge(unittest.TestCase):

    ns = [16, 30]
    iterations = 16

    def setUp(self):
        random.seed(0xBEEFCAFE)

    def test_int_addition(self):
        for n, m in combinations_with_replacement(self.ns, 2):
            a_reg = cirq.NamedQubit.range(n, prefix="a")
            b_reg = cirq.NamedQubit.range(m, prefix="b")
            ancillas = AncillaManager()

            c = cirq.Circuit(add_int(a_reg, b_reg, ancillas), strategy=NEW)
            sim = ToffoliSimulator(c)

            for _ in range(self.iterations):
                a = random.randint(0, 2**n-1)
                b = random.randint(0, 2**m-1)
                with self.subTest(a=a, b=b):
                    state = int_to_state(a, a_reg)
                    state.update(int_to_state(b, b_reg))
                    state.update(ancillas.init_state())
                    sim.simulate(state)
                    ra = state_to_int(state, a_reg)
                    rb = state_to_int(state, b_reg)

                    self.assertEqual(ra, a)
                    self.assertEqual(rb, (a+b)%(2**m))

    def test_classical_int_addition(self):
        for n, m in combinations_with_replacement(self.ns, 2):
            for _ in range(self.iterations):
                a = random.randint(0, 2**n-1)
                b = random.randint(0, 2**m-1)
                with self.subTest(a=a, b=b):

                    b_reg = cirq.NamedQubit.range(m, prefix="b")
                    ancillas = AncillaManager()

                    c = cirq.Circuit(add_classical_int(a, b_reg, ancillas), strategy=NEW)
                    sim = ToffoliSimulator(c)

                    state = int_to_state(b, b_reg)
                    state.update(ancillas.init_state())
                    sim.simulate(state)
                    rb = state_to_int(state, b_reg)

                    self.assertEqual(rb, (a+b)%(2**m))

    def test_square(self):

        square_methods = [
            ("schoolbook", schoolbook_square),
            ("karatsuba", lambda *args: karatsuba_square(*args, cutoff=4))
        ]

        for n in self.ns:
            for name, square in square_methods:
                with self.subTest(square_method=name):

                    a_reg = cirq.NamedQubit.range(n, prefix="a")
                    b_reg = cirq.NamedQubit.range(2*n, prefix="b")
                    ancillas = AncillaManager()

                    c = cirq.Circuit(square(a_reg, b_reg, ancillas), strategy=NEW)
                    sim = ToffoliSimulator(c)

                    for _ in range(self.iterations):
                        a = random.randint(0, 2**n-1)
                        b = random.randint(0, 2**(2*n)-1)
                        with self.subTest(a=a, b=b):
                            state = int_to_state(a, a_reg)
                            state.update(int_to_state(b, b_reg))
                            state.update(ancillas.init_state())
                            sim.simulate(state)
                            ra = state_to_int(state, a_reg)
                            rb = state_to_int(state, b_reg)

                            self.assertEqual(ra, a)
                            self.assertEqual(rb, (b+a**2)%(2**(2*n)))

    def test_mult(self):

        mult_methods = [
            ("schoolbook", schoolbook_mult),
            ("karatsuba", lambda *args: karatsuba_mult(*args, cutoff=4))
        ]

        for n,m in product(self.ns, repeat=2):
            for name, mult in mult_methods:
                with self.subTest(mult_method=name):

                    a_reg = cirq.NamedQubit.range(n, prefix="a")
                    b_reg = cirq.NamedQubit.range(m, prefix="b")
                    c_reg = cirq.NamedQubit.range(m+n, prefix="c")
                    ancillas = AncillaManager()

                    circ = cirq.Circuit(mult(a_reg, b_reg, c_reg, ancillas), strategy=NEW)
                    sim = ToffoliSimulator(circ)

                    for _ in range(self.iterations):
                        a = random.randint(0, 2**n-1)
                        b = random.randint(0, 2**m-1)
                        c = random.randint(0, 2**(n+m)-1)
                        with self.subTest(a=a, b=b, c=c):
                            state = int_to_state(a, a_reg)
                            state.update(int_to_state(b, b_reg))
                            state.update(int_to_state(c, c_reg))
                            state.update(ancillas.init_state())

                            sim.simulate(state)

                            ra = state_to_int(state, a_reg)
                            rb = state_to_int(state, b_reg)
                            rc = state_to_int(state, c_reg)

                            self.assertEqual(ra, a)
                            self.assertEqual(rb, b)
                            self.assertEqual(rc, (c+a*b)%(2**(n+m)))


class TestAncillas(unittest.TestCase):

    def test_new(self):
        a = AncillaManager()
        self.assertEqual(len(a), 0)

        q1 = a.new()
        self.assertTrue(isinstance(q1, cirq.NamedQubit))
        self.assertEqual(len(a), 1)

        q2 = a.new()
        self.assertTrue(isinstance(q2, cirq.NamedQubit))
        self.assertEqual(len(a), 2)

        self.assertNotEqual(q1, q2)

    def test_new_register(self):
        n = 5
        a = AncillaManager()

        qubits = a.new_register(n)
        self.assertEqual(len(qubits), n)
        self.assertTrue(all(isinstance(q, cirq.NamedQubit) for q in qubits))
        self.assertEqual(len(a), n)

        qubits2 = a.new_register(n)
        self.assertEqual(len(qubits), n)
        self.assertTrue(all(isinstance(q, cirq.NamedQubit) for q in qubits2))
        self.assertEqual(len(a), 2*n)

        self.assertTrue(not any(q in qubits for q in qubits2))


class TestUtils(unittest.TestCase):

    def test_int_state_conversion(self):
        n = 8
        test_cases = [84, 117, 225]
        r = cirq.NamedQubit.range(n, prefix="r")
        for x in test_cases:
            with self.subTest(x=x):
                self.assertEqual(x, state_to_int(int_to_state(x, r), r))


class TestMontgomery(unittest.TestCase):

    from math import gcd
    iters = 32

    def setUp(self):
        random.seed(0xFEEDBEEF)

    def test_extended_gcd(self):
        n = 16
        for _ in range(self.iters):

            # pick two random numbers with gcd 1
            a = random.randint(0, 2**n-1)
            b = random.randint(0, 2**n-1)
            g = self.gcd(a, b)
            a //= g
            b //= g

            x, y = extended_gcd(a, b)

            self.assertEqual(a*x % b, 1)
            self.assertEqual(b*y % a, a-1)

            self.assertGreaterEqual(x, 0)
            self.assertGreaterEqual(y, 0)
            self.assertLess(x, b)
            self.assertLess(y, a)

    def test_montgomery_reduce(self):
        n = 8

        mult_methods = [
            ("schoolbook", schoolbook_classical_mult),
            ("karatsuba", lambda *args, **kwargs: karatsuba_classical_mult(*args, cutoff=4, **kwargs))
        ]

        for _ in range(self.iters):
            N = random.randint(0, 2**n-1) | 1 | (1 << (n-1))# make sure it's odd and big
            T = random.randint(0, N*(2**n)-1)
            for name, mult in mult_methods:
                with self.subTest(T=T, N=N, mult_method=name):

                    T_reg = cirq.NamedQubit.range(2*n+1, prefix="T")
                    ancillas = AncillaManager()

                    R, circ_gen = montgomery_reduce(T_reg, ancillas, N, mult=mult)
                    circ = cirq.Circuit(circ_gen, strategy=NEW)
                    sim = ToffoliSimulator(circ)

                    state = int_to_state(T, T_reg)
                    state.update(ancillas.init_state())

                    sim.simulate(state)

                    rT = state_to_int(state, T_reg[n:])
                    rT_low = state_to_int(state, T_reg[:n])

                    self.assertEqual(rT_low, 0)
                    self.assertEqual((rT*R)%N, T%N)
                    self.assertLess(rT, N)

if __name__ == '__main__':
    unittest.main()
