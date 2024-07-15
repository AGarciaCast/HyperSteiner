# HyperbolicSteinerTrees
This repository contains the code for the paper called **"HyperSteiner: Computing Heuristic Hyperbolic Steiner Minimal Trees"**.

# Abstract
We propose HyperSteiner — an efficient heuristic algorithm for computing Steiner minimal trees in the hyperbolic space. Hypersteiner extends the Euclidean Smith-Lee-Liebman algorithm, which is grounded in a divide-and-conquer approach involving the Delaunay triangulation. The central idea is rephrasing Steiner tree problems with three terminals as a system of equations in the Klein-Beltrami model. Motivated by the fact that hyperbolic geometry is well-suited for representing hierarchies, we explore applications to hierarchy discovery in data. Results show that HyperSteiner infers more realistic hierarchies than the Minimum Spanning Tree and is more scalable to large datasets than Neighbor Joining.

# Information
- The code for the method is in `src/utils`. <br />
- The synthetic and real experiments can be reproduced using respectively the notebooks `syntheticExperiments.ipynb` and `realExperiments.ipynb`. <br />
- The folder `Data` contains the hyperbolic representation of the Planaria dataset in the Klein-Beltrami disk and the corresponding groundtruth cell age. This is used for the real-life experiments.
