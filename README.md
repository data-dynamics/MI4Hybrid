#### MI4Hybrid

MI4Hybrid is a model invalidation toolbox for hybrid systems.

####Installation Instructions:

This toolbox can be used in MATLAB with the following necessary packages/softwares installed:

1. [YALMIP](http://users.isy.liu.se/johanl/yalmip/)
2. [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/)
3. [SparsePOP] (http://www.is.titech.ac.jp/~kojima/SparsePOP/)

 
####Contents:

#####In the folder "lib":

1. StateSpace.m is a class for state-space models with input, state, and noise specifications.
   * UnStateSpace.m is a class for state-space models with parameter uncertainty in addition.
   * UnStateSpace.m is a subclass of StateSpace.m

2. ARXmodel.m is a class for ARX models with input and noise specifications.
   * UnARXmodel.m is a class for ARX models with parameter uncertainty in addition.
   * UnARXmodel.m is a subclass of ARXmodel.m

3. PolyModel.m is a class for (non-switched) polynomial models with noise specifications.
   * UnPolyModel.m is a class for (non-switched) polynomial models parameter uncertainty in addition.
   * UnPolyModel.m is a subclass of PolyModel.m

4. bounded_noise.m is a function generating l_p norm bounded noise (a matrix) whose number of rows is the noise dimension and number of columns is the time horizon.

5. swarx_sim.m is a function that generates simulated I/O data for ARX models defined on ARXmodel.m or UnARXmodel.m
 
6. swss_sim.m is a function that generates simulated I/O data for state-space models defined on StateSpace.m or UnStateSpace.m

7. poly_sim.m is a function that generates simulated I/O data for polynomial models defined on PolyModel.m or UnPolyModel.m

8. invalidation_arx.m is a function that applies an invalidation algorithm to non-switched ARX models.
 
9. invalidation_ss.m is a function that applies an invalidation algorithm to non-switched state-space models.

10. invalidation_sarx_milp.m is a function that applies an invalidation algorithm to any switched or non-switched ARX models.

11. invalidation_swa_milp.m is a function that applies an invalidation algorithm to any switched or non-switched state-space models.

12. invalidation_poly.m is a function that appplies an invalidation algorithm to any certain or uncertain polynomial state-space models.

13. tdet_poly.m is a function that checks whether a fault model sysf is T-detectable for a system model sys for a given T.

#####In the folder "examples":

Examples for switched, non-switched ARX/state-space and polynomial model invalidation using different functions.

#####In the folder "extras":

Extra files are inside this folder.


####Related publications:
1. F. Harirchi and N. Ozay, "Model Invalidation for Switched Affine Systems with Applications to Fault and Anomaly Detection", IFAC ADHS, 2015.

2. F. Harirchi, Z. Luo and N. Ozay, "Model (In)validation and Fault Detection for Systems with Polynomial State-Space Models", ACC, 2016.

####Acknowledgments:
This research is supported in part by DARPA grant N66001-14-1-4045.
