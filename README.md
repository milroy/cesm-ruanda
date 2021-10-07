# cesm-ruanda

### CESM Root caUse Analysis of Numerical DiscrepAncy (CESM-RUANDA)
Welcome to the repository for code developed to locate 
the root causes of numerical discrepancy in the Community 
Earth System Model ([CESM&reg;](https://www.cesm.ucar.edu/)). 
The methods are described in detail in 
[Making Root Cause Analysis Feasible for Large Code Bases: A Solution Approach for a Climate Model](https://doi.org/10.1145/3307681.3325399),
and appear in [Keeping science on keel when software moves](https://doi.org/10.1145/3382037).

### Building the RUANDA Docker images
We are updating the code to make it more robust and 
reproducible with Docker. To build the RUANDA AST analysis image, run: 

```
git/cesm-ruanda$ docker build -t <image_name> -f ruanda.docker .
```

To build the ML analysis container, run:
```
git/cesm-ruanda$ docker build -t <image_name> -f ml-tools.docker . 
```

### Running the RUANDA Analysis

Following the workflow described in Making Root Cause Analysis Feasible for Large Code Bases: A Solution Approach for a Climate Model, 
you identify the CAM variables most affected by the cause(s) of the discrepancy.  
To do so using the L1-penalized logistic regression for the RANDMT experiment, execute:
```
docker run -it <image_name> python3 /root/cesm-ruanda/ml-tools/lassovars.py --ensSumUF /root/cesm-ruanda/data/sz750-ufect_T.nc --ensSumYR /root/cesm-ruanda/data/sz300.intel-gnu-pgi-rand1_V6.nc --ensemble /root/cesm-ruanda/data/cheyenne_noavx_dict_10ts_119.pickle --test /root/cesm-ruanda/data/uf_chey_randmt_10ts_119.pickle --timestep 0 --nRuns 119 --uf --regCoef 0.011 --standardize
```
which in the case of the RANDMT experiment, will produce output similar to the following: 
```
[...]
Ensemble betas:  [(b'snowhlnd', 0.19449419367867737), (b'flns', 0.0965858789787028), (b'qrl', 0.024443707174923203)]
Experiment betas:  [(b'flds', -0.09658587866942345), (b'taux', -0.06903271431400505)]
```
These are the CAM output variables identified by the technique as being "most-affected", together with their coefficient values. CESM does not use these names internally, so we need to find their internal names before proceeding to the AST directed-graph analysis.  You can find the mapping in `CAM-var-mapping.csv`.  
Note that you should use the name after the Fortran derived type selector, e.g.: `flwds` for `cam_out%flwds`.  
See Making Root Cause Analysis Feasible for Large Code Bases: A Solution Approach for a Climate Model, Section 4.2 for details.

With the identified variables and their mapped names, we can now run the AST analysis.  
To do so, execute (relative to the git/cesm-ruanda directory), e.g., for the RANDMT experiment:

```
$ docker run -it <image_name>:latest python -i cesm-ruanda/ast-processing/ast_graph.py --src cesm-ruanda/ast-processing/randmt/randmtsrc
```
You will be dropped into a Python REPL and further analysis can be performed 
on the CESM metagraph. To get the Girvan-Newman communities, variable eigencentralities, CAM modules, etc, you can run:
```
>>> ordered_mods, ordered_mods_procedures, ordered_mods_procedures_eigen, subgraph_nodes, eigens, communities = CESM.eigenCentrality(['flwds', 'wsx', 'snowhland', 'flns', 'qrl'], models=['cam'], girvan=True)
```
Where `ordered_mods` are the modules ordered by maximum contained node eigencentrality, `ordered_mods_procedures` 
also includes the subroutine or functions (ordered by call depth), `ordered_mods_procedures_eigen` includes the 
eigenvalue of the node with the max eigenvalue, `subgraph_nodes` are all the nodes in the induced subgraph, 
`eigens` are the subgraph nodes sorted in descending order by their eigencentralities, and `communities` 
are the communities detected by the Girvan-Newman algorithm.

There are many more functions available within the CESM metagraph, which you can print via ```>>> help(CESM)```.


### Details on the ML techniques and distribution comparison
#### `lassovars.py`
This technique uses L1-penalized logistic regression to identify CAM output variables "most affected" by the experimental change.  The Python script uses scikit learn's [LogisticRegression](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html?highlight=logisticregression#sklearn.linear_model.LogisticRegression) linear model with L1 penalty and the saga solver.  The following arguments are required:
```
--test the file with the experimental output data
--ensemble the file with the ensemble output data
--ensSumUF the Ultra-Fast ensemble summary file
--ensSumYR the yearly ensemble summary file
--timestep the CAM time step at which to compare the outputs
--nRuns the number of runs to use
--regCoef the L1 penalty to use (smaller values identify fewer variables)
--uf (if specified) these are Ultra-Fast runs
--yr (if specified) these are yearly runs
```

#### `distributions.py`
This technique compares CAM output variable distributions between a set of experimental runs and a set of ensemble runs.  After standardizing the variables, it outputs variables with "non-overlapping" distributions, i.e., those whose interquartile ranges do not overlap.  The following arguments are required:
```
--sumFile the ensemble summary file
--testFile the file with the test output data
--tstep the CAM time step at which to compare the variables' distributions
--ensFile the file with the ensemble output data
```
For example, here is the output for the RANDMT experiment at time step 0:
```
$ docker run -it <image name> python3 /root/cesm-ruanda/ml-tools/distributions.py --sumFile /root/cesm-ruanda/data/sz750-ufect_T.nc --testFile /root/cesm-ruanda/data/uf_chey_randmt_10ts_119.pickle --tstep 0 --ensFile /root/cesm-ruanda/data/cheyenne_noavx_dict_10ts_119.pickle
[(b'TMQ', 777.0421823282961), (b'FSNSC', 287.6456767207707), (b'FSNS', 287.6456760629575), (b'bc_a1_SRF', 221.58333619824262), (b'QRL', 220.00067407167742), (b'TGCLDIWP', 178.0262119461754), (b'ncl_a3_SRF', 144.4378355718746), (b'dst_a1_SRF', 134.04502787376947), (b'num_a2_SRF', 133.59768585294722), (b'num_a1_SRF', 114.73536285732408), (b'QRS', 61.38002207415759), (b'FSNTC', 52.82121615988464), (b'FSNTOAC', 52.82116837309318), (b'TS', 43.230601158076325), (b'DTV', 23.890325598640633), (b'num_a3_SRF', 23.710855497615135), (b'TGCLDLWP', 21.247669139828467), (b'pom_a1_SRF', 18.654890044638314), (b'VD01', 16.757755789754963), (b'U10', 14.257064385019461), (b'H2SO4_SRF', 13.698788418582211), (b'WGUSTD', 13.69239722308186), (b'PBLH', 13.60494497457398), (b'TAUY', 12.299771128420858), (b'FSDSC', 10.633705995788764), (b'DTWR_H2SO4', 10.031753734791826), (b'QREFHT', 9.448450417026466), (b'PRECSL', 7.384330307722539), (b'PRECL', 7.278734297411136), (b'SHFLX', 6.849387002354607), (b'DTWR_H2O2', 5.7230075691333075), (b'DTWR_SO2', 4.925471168874313), (b'QFLX', 3.7075765692867098)]
```
The variables are sorted in descending order by the difference between the ensemble and test medians.  Since the variables are all standardized, the top ranking variables by median difference exhibit the greatest difference according to this technique. Notice that there are many more variables identified by this technique.
