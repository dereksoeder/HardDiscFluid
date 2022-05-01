# Hard-disc fluid

Basic event-driven molecular dynamics (EDMD) simulation of a system of rigid, (mostly) smooth discs, to explore the hydrodynamic limit of a simple, classical, many-body system.
The test configurations generally follow Ishiwata _et al._, [Int. J. Mod. Phys. C __15__ (2004) 1413-1424](https://doi.org/10.1142/S0129183104006820).

## Build

Download the source files, or clone this repository:

    git clone https://github.com/dereksoeder/HardDiscFluid
    cd HardDiscFluid

On Linux, or more generally an environment with a GNU C++ Compiler (supporting C++11 or newer) installed, go to the directory containing the source files and run the following command:

    g++ *.cxx

An executable file (typically named `a.out`) will be created in the same directory if the build is successful.

## Run

After building the executable (`a.out` in these instructions), run it with the desired test number (_`testnumber`_) and any required arguments (_`args`_):

<pre><code>./a.out <i>testnumber</i> <i>args...</i></code></pre>

## Tests

Currently two tests are implemented:

### Test 0: Plane Poiseuille flow

Simulates a flow of discs between two "hot" walls (i.e., diffuse-reflection boundary conditions, to realize no-slip) on top and bottom, with periodic boundaries at left and right, and a constant, uniform gravitational acceleration to the right imitating a constant pressure gradient.  The steady-state flow under these conditions is [plane Poiseuille flow](https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_equation#Plane_Poiseuille_flow), a  profile of horizontal flow that is parabolic in _y_.

* **Usage:**  `./a.out 0 `_`density`_, where _`density`_ is a number from 0.2 to 0.5 expressing the number density of discs in the arena.
* **Output:**  A sequence of 30 rows of two fields (two columns) each: a _y_ coordinate followed by the mean horizontal flow of all discs observed in that horizontal strip or "lane" of the arena, averaged over all simulation time excluding an initial burn-in.
* **Function:**  `TestPlanePoiseuille`, in [`main.cxx`](main.cxx)
* **Notes:**  Expect this test to take a few hours depending on the density specified.  Currently there is no status or progress output.

### Test 1: Flow past a cylinder

Simulates a flow of discs around a single large, immovable disc referred to as the cylinder, which has a special "erratic" boundary condition that randomizes the reflection angles of impinging discs to effect no-slip.  The arena consists of two "hot" walls on top and bottom, and "hot" periodic boundaries at left and right that assign random thermally-distributed velocities to rightward-crossing discs; in all cases, thermal velocities are boosted horizontally by the inflow velocity.

* **Usage:**  `./a.out 1 `_`Re`_, where _`Re`_ is one of the following Reynolds numbers: 6, 11, 33, 66, or 106.  These Reynolds numbers correspond to configurations (arena size, cylinder placement, inflow velocity, disc density) adapted from Ishiwata _et al._
* **Output:**  A "header"-style line of the form `# `_`dragforce`_, where _`dragforce`_ is the average horizontal momentum exchanged with the cylinder per unit time, followed by _m_ lines of 2_n_ fields each, constituting an _m_ &times; _n_ &times; 2 grid expressing the mean flow velocity (&langle;_v<sub>x</sub>_&rangle;, &langle;_v<sub>y</sub>_&rangle;) in each 20 a.u. &times; 20 a.u. cell of the arena.  Both the drag force and the mean flow are averaged over all simulation time excluding an initial burn-in.
* **Function:**  `TestCylinder`, in [`main.cxx`](main.cxx)
* **Notes:**  Expect this test to take more than a day to complete.  Currently there is no status or progress output.
