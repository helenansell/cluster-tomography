# cluster-tomography

Python script to calculate the nonlinear contribution to the cluster number count in percolation for each of the ten line configurations in a 2d square system, as illustrated in Fig. 1 of Ansell et al., arXiv:2307.04260 (2023).

**percolation_cluster_tomography.py** is the main script, from which the line configuration (1-10), percolation type (site or bond), system size, and number of independent samples can be selected. It returns the mean and standard error of the nonlinear contribution to the cluster number count, averaged across the independent samples. For cases where the line segment can be placed at different angles, the result is returned for $\gamma = \arctan(1/n)$ for integer $n \leq 20$. 

**func_perc_samples.py** contains functions called by the main script that generate 2d site and bond percolation samples, and apply periodic boundary conditions as required.

**func_line_clusters.py** contains functions called by the main script that determine the lattice sites along lines through the system at a given angle, and count the number of clusters along a given line segment.

### Cite as: 
H. S. Ansell, S. J. Frank and I. A. Kovács, Cluster Tomography in Percolation, arXiv:2307.04260 (2023), [doi:10.48550/arXiv.2307.04260](https://doi.org/10.48550/arXiv.2307.04260).

