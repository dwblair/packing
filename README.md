<<<<<<< HEAD
## About

Statistical physics-y approaches to finding optimal packings (highest density arrangements of shapes within shapes)

Initially simple python code; the hope is to extend this with cython.
 
## Algorithm

### Verison 0.1 (current)
- Two types of moves: a) translation (of an individual circle); b) expansion / contraction (of all circles).  At each step, the algorithm chooses one of these moves (type "a" being much more likely); in the case of move type a), a random particle is chosen for an attempted random translation -- if it overlaps with another particle, then the move is rejected.  In the case of move type b), an attempt is made to expand or contract all particles (biased towards expansion); if no overlap, move is accepted (which is always the case if contraction)

### Version 0.2:
- Nearest neighbor list, or cell, overlap checking
- Contraction should be accepted with probability P = exp(-Beta dV) -- dV is change in volume of system
- Expansion and translation amounts should be periodically adjusted depending on relative size of particles in system -- do a 'thermalization' of these paramters, then do 'equilibration' ... repeat ...
- Parallel tempering -- do more than one system at a time
- Population annealing -- implement this algorithm
=======
irkit
=====

Raspberry Pi-based infrared/visible image compositing and webcam control (Public Lab)

* fswebcam must be installed
* try running "python snap.py -n 5 -d 5" to take 5 snapshots, separated by several seconds ...
>>>>>>> bf566dc8dea64c10bbea65ecae600ab63cf3b1bf
