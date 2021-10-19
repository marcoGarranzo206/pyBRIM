# pyBRIM
Algorithms for community detection on bipartite graphs based on BRIM algorithm

# Overview

Finding communities is a very common and important step in most network analysis. Many community finding algorithms are based on Newmans modularity optimization, initially design for "normal" networks (unweighted, unsigned). Here I implement algorithms for community finding in python for specifically bipartite networks, signed or unsigned, based on the BRIM the algorithm. There is additionally a meta clustering module for combining results from different algorithm runs.

# Example usages:

#### BRIM. Unknown number of communities, unsigned
```python
import networkx as nx
import pyBRIM

g = nx.davis_southern_women_graph()
bs = pyBRIM.BRIM_solver(g)
R_t,S, Q_max = bs.BRIM_bisec()
nodes_to_communities = bs.translate_communities(R_t,S)
```

#### BRIM. Known number of communities, unsigned
```python
import networkx as nx
import pyBRIM

g = nx.davis_southern_women_graph()
bs = pyBRIM.BRIM_solver(g)
R_t,S, Q_max = bs.fit_transform(4)
nodes_to_communities = bs.translate_communities(R_t,S)
```
#### BRIM. Unknown number of communities, signed

#### Meta clustering: HAC

#### Meta clustering: Majprity voting

#### Meta clustering: best meta clusters

# Description

# References

