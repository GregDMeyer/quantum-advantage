
import unittest

from itertools import product
import cirq
from circuits import full_adder, add_int
from tof_sim import ToffoliSimulator, int_to_state, state_to_int
from ancilla import AncillaManager

class TestSubCircuits(unittest.TestCase):

    def test_adder(self):
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
                
    def test_int_addition(self):
        n = 5
        test_cases = [
            (0, 0),
            (0, 1),
            (1, 0),
            (1, 1),
            (3, 0),
            (18, 30),
            (27, 27),
            (31, 20),
            (1, 31),
            (30, 12)
        ]

        c = cirq.Circuit()
        a_reg = cirq.NamedQubit.range(n, prefix="a")
        b_reg = cirq.NamedQubit.range(n, prefix="b")
        ancillas = AncillaManager()
        
        add_int(c, a_reg, b_reg, ancillas)

        sim = ToffoliSimulator(c)
        
        for a, b in test_cases:
            with self.subTest(a=a, b=b):
                state = int_to_state(a, a_reg)
                state.update(int_to_state(b, b_reg))
                state.update(ancillas.init_state())
                sim.simulate(state)
                ra = state_to_int(state, a_reg)
                rb = state_to_int(state, b_reg)

                self.assertEqual(ra, a)
                self.assertEqual(rb, (a+b)%(2**n))

    def test_int_addition_exceptions(self):
        c = cirq.Circuit()
        a_reg = cirq.NamedQubit.range(5, prefix="a")
        b_reg = cirq.NamedQubit.range(6, prefix="b")
        ancillas = AncillaManager()

        # this should not raise an exception
        add_int(c, a_reg, b_reg, ancillas)
                    
        with self.assertRaises(ValueError):
            add_int(c, b_reg, a_reg, ancillas)

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
