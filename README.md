## Quantum circuits for *x*<sup>2</sup> mod *N*

This repository contains implementations for four different quantum circuits that implement the unitary 𝒰 defined by

𝒰\|*x*⟩\|0⟩ = \|*x*⟩\|*x*<sup>2</sup> mod *N*⟩.

This unitary is the key step in the "proof of quantumness" protocol described in this paper: [arXiv:2104.00687](https://arxiv.org/abs/2104.00687).

The repository also contains code for analyzing the circuits.

#### Circuit constructions

The code implements four circuits. The first two are quantum implementations of classical multiplication circuits:

 - Quantum implementation of [Karatsuba multiplication](https://en.wikipedia.org/wiki/Karatsuba_algorithm)
 - Quantum implementation of ["schoolbook" multiplication](https://en.wikipedia.org/wiki/Multiplication_algorithm#Long_multiplication)

These circuits take advantage of the advantage protocol's ability to discard garbage bits without uncomputing (see manuscript).

The other two circuits compute the value in the phase, and then transfer it to the register using an inverse quantum Fourier transform. See manuscript for details.

### Requirements

 The code has been tested with Python 3.9, though should be compatible with other recent Python 3 versions. The dependencies can be installed with

```bash
pip3 install -r requirements.txt
```

### Usage

#### Implementations

The `circuits` directory contains files containing the functions to generate the circuits, as well as ancillary functions for circuit construction and simulation.

In `circuits/digital_circuits.py`, we define functions for generating the quantum circuits based on Karatsuba and schoolbook classical multiplication algorithms.

In `circuits/phase_circuits.py`, we define functions for generating quantum circuits that use controlled-controlled-phase gates to compute the value in phase space, and then perform an inverse quantum Fourier transform to transfer this phase into the register. See manuscript for details.

Examples of usage of the functions defined in these files can be seen in `analysis/count_gates.py` (see next section).

#### Analysis

In the `analysis` directory, there are three scripts:

 - `count_gates.py`: compute gate count, gate depth, and other relevant quantities for each of the implementations
 - `error_analysis.py`: analyze the performance of the quantum advantage protocol when subject to noise, including the performance of the postselection scheme described in the manuscript
 - `optimize_cutoff.py`: optimize the cutoff for recursion in the Karatsuba algorithm

#### Tests

 The `tests` directory contains tests.

 The "digital" circuits can be tested (using the `unittest` framework) by running `python3 test_digital_circuits.py`.

 The "phase" circuits can be checked using the `check_phase_circuits.py` script, which does not explicitly test anything itself, but outputs the accuracy of the circuit in reproducing the desired state (those circuits are approximate).