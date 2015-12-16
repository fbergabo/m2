Phys 2200 midterm 2 project

In this program, we use the ising model to determine the temperature of the phase transition for a trangular two-dimensional lattice. This is done using a metropolis simulation which calculates the interaction energy between the nearest neightbors for each spin within the lattice. For each temperature, a ratio is created from Boltzmann factors. This ratio is used as the probability that a spin change that increases energy will be accepted. All spin changes that decrease energy will always be accepted. We use a "cold start", meaning that all the spins are intially aligned. Using the GSL random number generator taus2, we are able to simulate the apparent randomness of a real material. 

This program is an adaptation of a code used for the a similar situation with a square lattice. In the square lattice, each spin has four neighbors with which it interacts. In the triangular lattice, this has been increased to six neighbors. This change is also reflected in the calculation of the average energy per spin. In order to determine the temperature of the phase transition, the program runs the metropolis simulation for several different temperatures. The specific heat is calculated with each step, and the specific heat is then automatically plotted against the temperature once the program is completed. The phase transition temperature corresponds to the peak in the specific heat. The graph produced by this program is shown below:

![](graph.png?)

This program was in large part a collaboration with Connor Occhialini.
