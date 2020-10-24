
import unittest

from itertools import product
import cirq
from circuits import full_adder, half_adder, add_int, square
from tof_sim import ToffoliSimulator, int_to_state, state_to_int
from ancilla import AncillaManager

class TestSubCircuits(unittest.TestCase):

    def test_half_adder(self):
        qubits = cirq.LineQubit.range(3)
        c = cirq.Circuit()
        half_adder(c, *qubits)
        sim = ToffoliSimulator(c)

        for a, b in product((0, 1), repeat=2):
            with self.subTest(a=a, b=b):
                init_state = dict(zip(qubits, (a, b, 0)))
                result = sim.simulate(init_state)
                ra, rsum, rcout = [result[q] for q in qubits]

                self.assertEqual(ra, a)
                self.assertEqual(rsum, (a+b)%2)
                self.assertEqual(rcout, (a+b)//2)

    def test_full_adder(self):
        qubits = cirq.LineQubit.range(4)
        c = cirq.Circuit()
        full_adder(c, *qubits)
        sim = ToffoliSimulator(c)

        for a, b, cin in product((0, 1), repeat=3):
            with self.subTest(a=a, b=b, cin=cin):
                init_state = dict(zip(qubits, (a, b, cin, 0)))
                result = sim.simulate(init_state)
                ra, rsum, rcin, rcout = [result[q] for q in qubits]

                self.assertEqual(ra, a)
                self.assertEqual(rsum, (a+b+cin)%2)
                self.assertEqual(rcin, cin)
                self.assertEqual(rcout, (a+b+cin)//2)

class TestArithmetic(unittest.TestCase):
                
    def test_int_addition_same(self):
        test_cases = [
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

        for n, m, a, b in test_cases:
            with self.subTest(a=a, b=b):

                c = cirq.Circuit()
                a_reg = cirq.NamedQubit.range(n, prefix="a")
                b_reg = cirq.NamedQubit.range(m, prefix="b")
                ancillas = AncillaManager()
        
                add_int(c, a_reg, b_reg, ancillas)

                sim = ToffoliSimulator(c)
        
                state = int_to_state(a, a_reg)
                state.update(int_to_state(b, b_reg))
                state.update(ancillas.init_state())
                sim.simulate(state)
                ra = state_to_int(state, a_reg)
                rb = state_to_int(state, b_reg)

                self.assertEqual(ra, a)
                self.assertEqual(rb, (a+b)%(2**m))

    def test_int_addition_exceptions(self):
        c = cirq.Circuit()
        a_reg = cirq.NamedQubit.range(5, prefix="a")
        b_reg = cirq.NamedQubit.range(6, prefix="b")
        ancillas = AncillaManager()

        # this should not raise an exception
        add_int(c, a_reg, b_reg, ancillas)
                    
        with self.assertRaises(ValueError):
            add_int(c, b_reg, a_reg, ancillas)

    def test_square(self):
        n = 5
        test_cases = range(2**n)

        c = cirq.Circuit()
        a_reg = cirq.NamedQubit.range(n, prefix="a")
        b_reg = cirq.NamedQubit.range(2*n, prefix="b")
        ancillas = AncillaManager()

        with self.assertRaises(ValueError):
            square(c, a_reg, a_reg, ancillas)
        
        square(c, a_reg, b_reg, ancillas)

        sim = ToffoliSimulator(c)
        
        for a in test_cases:
            for b in (0, 2**n-1):
                with self.subTest(a=a, b=b):
                    state = int_to_state(a, a_reg)
                    state.update(int_to_state(b, b_reg))
                    state.update(ancillas.init_state())
                    sim.simulate(state)
                    ra = state_to_int(state, a_reg)
                    rb = state_to_int(state, b_reg)

                    self.assertEqual(ra, a)
                    self.assertEqual(rb, b+a**2)

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


class TestUtils(unittest.TestCase):

    def test_int_state_conversion(self):
        n = 8
        test_cases = [84, 117, 225]
        r = cirq.LineQubit.range(n)
        for x in test_cases:
            with self.subTest(x=x):
                self.assertEqual(x, state_to_int(int_to_state(x, r), r))
        
if __name__ == '__main__':
    unittest.main()
