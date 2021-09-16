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


### Details on the ML techniques and boxplot generation
```
docker run -it <image_name> python3 /root/cesm-ruanda/ml-tools/lassovars.py --ensSumUF /root/cesm-ruanda/data/sz750-ufect_T.nc --ensSumYR /root/cesm-ruanda/data/sz300.intel-gnu-pgi-rand1_V6.nc --ensemble /root/cesm-ruanda/data/cheyenne_noavx_dict_10ts_119.pickle --test /root/cesm-ruanda/data/uf_chey_randmt_10ts_119.pickle --timestep 0 --nRuns 119 --uf --regCoef 0.011 --standardize
```
