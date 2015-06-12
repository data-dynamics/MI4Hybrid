# MI4Hybrid2 (An Informal Version)

1. StateSpace.m is a class for state space models with noise specification.
   * UnStateSpace.m is a class for state space models with noise specification and parameter uncertainty.
   * UnStateSpace.m is a subclass of StateSpace.m

2. ARXmodel is a class for ARX models with noise specification.
   * UnARXmodel.m is a class for ARX models with noise specification and parameter uncertainty.
   * UnARXmodel.m is a subclass of ARXmodel.m

3. polymodel.m is a class for (not switched) polynomial models with noise specification.
   * Unpolymodel.m is a class for (not switched) polynomial models with noise specification and parameter uncertainty.
   * Unpolymodel.m is a subclass of polymodel.m

4. bounded_noise.m is a function generates l_p norm bounded noise (a matrix) whose number of rows is the noise dimension and number of columns is the time horizon.

5. simulates.m is a simulation function which will simulate the systems using the models we made.

6. ARX_Invalidation_Single.m is a file running for testing the algorithm of system invalidation. The system model used in this code is an ARX model with a single mode (i.e. only one submodel).
