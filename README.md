# Many-Body-Problems
Pre-req modules: `pygame.py`

A pygame-based visual simulation of bodies interacting under an inverse-square gravitational force law.

The python source code included in the `ManyBody.py` commit produces a simulation of the forward-time evolution of a number of bodies (with given initial positions and momenta) interacting under something like a Newtonian graviational potential.

Bodies are represented as window-blitted circles whose radius is proportional to the mass used to calculated gravitational forces.

KNOWN BUGS: An attracting steady-state 'glued' onto the edge of the screen can occur for some initial conditions.
