# cesm-ruanda

### CESM Root caUse Analysis of Numerical DiscrepAncy (CESM-RUANDA)
Welcome to the repository for code developed to locate 
the root causes of numerical discrepancy in the Community 
Earth System Model ([CESM&reg;](https://www.cesm.ucar.edu/)). 
The methods are described in detail in 
[Making Root Cause Analysis Feasible for Large Code Bases: A Solution Approach for a Climate Model](https://doi.org/10.1145/3307681.3325399),
and appear in [Keeping science on keel when software moves](https://doi.org/10.1145/3382037).

### Running RUANDA in Docker
We are updating the code to make it more robust and 
reproducible with Docker. To build the RUANDA image, run: 

```
git/cesm-ruanda$ docker build -t <image_name> -f ruanda.yaml .
```
After successfully building the Docker image, execute the following 
to run the AST analysis (relative to the git/cesm-ruanda directory), e.g., for the RANDMT experiment:

```
$ docker run -it <image_name>:latest python -i cesm-ruanda/ast-processing/ast_graph.py --src cesm-ruanda/ast-processing/randmt/randmtsrc
```
You will be dropped into a Python REPL and further analysis can be performed 
on the CESM metagraph.
