#### MI4Hybrid

MI4Hybrid is a model invalidation toolbox for hybrid systems.

####Installation Instructions:

This toolbox can be used in MATLAB with the following necessary packages/softwares installed:
* [YALMIP](http://users.isy.liu.se/johanl/yalmip/)

For Polynomial State-Space Model Invalidation and T-Detectability: 
* [SparsePOP] (http://www.is.titech.ac.jp/~kojima/SparsePOP/)

For Switched Affine and Reggressive Model Invalidation and T-Detectability: 
* [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) or [GUROBI](http://www.gurobi.com/)
 
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

4. PWAModel.m is a class for (non-switched) piece-wise affine models without noise

5. bounded_noise.m is a function generating l_p norm bounded noise (a matrix) whose number of rows is the noise dimension and number of columns is the time horizon.

6. swarx_sim.m is a function that generates simulated I/O data for ARX models defined on ARXmodel.m or UnARXmodel.m
 
7. swss_sim.m is a function that generates simulated I/O data for state-space models defined on StateSpace.m or UnStateSpace.m

8. poly_sim.m is a function that generates simulated I/O data for polynomial models defined on PolyModel.m or UnPolyModel.m

9. pwa_sim.m is a function that generates simulated I/O data for piece-wise affine models defined on PWAModel.m

10. invalidation_arx.m is a function that applies an invalidation algorithm to non-switched ARX models.
 
11. invalidation_ss.m is a function that applies an invalidation algorithm to non-switched state-space models.

12. invalidation_sarx_milp.m is a function that applies an invalidation algorithm to any switched or non-switched ARX models.

13. invalidation_swa_milp_old.m is a function that applies an invalidation algorithm to any switched or non-switched state-space models.

14. invalidation_swa_milp.m is a function that applies an invalidation algorithm to any switched or non-switched state-space models.

15. invalidation_uswa_milp_old.m is a function that applies an invalidation algorithm to any switched or non-switched state-space models subject to parameter uncertainty. (old version, which works with the old examples)

16. invalidation_poly.m is a function that applies an invalidation algorithm to any certain or uncertain polynomial state-space models.

17. invalidation_pwa.m is a function that applies an invalidation algorithm to any certain piece-wise affine models.

18. Tdetect_swa_milp.m is a function that checks whether an SWA fault model sysf is T-detectable for an SWA system model sys for a given T.

19. Tdetect_uswa_milp.m is a function that checks whether an uncertain SWA fault model sysf is T-detectable for an uncertain SWA system model sys for a given T.

20. tdet_poly.m is a function that checks whether a polynomial fault model sysf is T-detectable for a polynomial system model sys for a given T.

#####In the folder "examples":

* Examples for switched, non-switched ARX/state-space, polynomial and PWA model invalidation using different functions.

* Examples for fault detection in SWA, uncertain SWA models as well as weak detectability are included.

* Examples of the publications cited below are also provided.

#####In the folder "extras":

Extra files are inside this folder.
1. SMT-based T-Detectability functions are located here.


####Related publications:
1. F. Harirchi and N. Ozay, "Model Invalidation for Switched Affine Systems with Applications to Fault and Anomaly Detection", IFAC ADHS, 2015.

2. F. Harirchi, Z. Luo and N. Ozay, "Model (In)validation and Fault Detection for Systems with Polynomial State-Space Models", ACC, 2016.

3. F. Harirchi and N. Ozay, "[Guaranteed Model-Based Fault Detection in Cyber-Physical Systems: A Model Invalidation Approach](https://arxiv.org/abs/1609.05921)", arXiv:1609.05921, 2016.

####Acknowledgments:
This research is supported in part by DARPA grant N66001-14-1-4045.
