# HyperbolicSteinerTrees
This repository contains the code for the paper called **"HyperSteiner: Computing Heuristic Hyperbolic Steiner Minimal Trees"**.

# Abstract
We propose HyperSteiner—an efficient heuristic algorithm for computing Steiner minimal trees in the hyperbolic space. Our method extends the Euclidean Smith-Lee-Liebman algorithm, which is grounded in a divide-and-conquer approach involving the Delaunay triangulation. Our central idea is rephrasing Steiner tree problems with three terminals as a system of equations in the Klein-Beltrami model. Since hyperbolic geometry is well-suited for representing hierarchical data, we explore applications to hierarchy discovery. Results show that HyperSteiner performs better than the Minimum Spanning Tree and, for large datasets, is more efficient than Neighbor Joining.

# Information
- The code for the method is in `src/utils`. <br />
- The synthetic and real experiments can be reproduced using respectively the notebooks `syntheticExperiments.ipynb` and `realExperiments.ipynb`. <br />
- The folder `Data` contains the hyperbolic representation of the Planaria dataset in the Klein-Beltrami disk and the corresponding groundtruth cell age. This is used for the real-life experiments.
