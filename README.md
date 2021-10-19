# pyBRIM
Algorithms for community detection on bipartite graphs based on BRIM algorithm

# Overview

Finding communities is a very common and important step in most network analysis. Many community finding algorithms are based on Newmans modularity optimization, initially design for "normal" networks (unweighted, unsigned). Here I implement algorithms for community finding in python for specifically bipartite networks, signed or unsigned, based on the BRIM algorithm. There is additionally a meta clustering module for combining results from different algorithm runs.

# Usages and examples:

There are two modules, BRIM and meta.

## BRIM
BRIM module revolves around the BRIM_solver object.Pass it the graph you want to use, the null model (for signed or unsigned networks) and the resolution value (to determine xommunity sizes, more in dexcription).

After creation you can run the fit_transform method to find communities when the number of communities is found. If you dont know, you can use the BRIM_bisec method, which you can pass it a maximum number of communities.

## meta
meta module revolves around the metaCluster object. You pass it an array of cluster runs. Each row corresponds to a cluster run, each column a node. (i,j) element s the community of node j in cluster run i. A co-occurrence matrix is built.

After creation you can do meta clustering based on majority voting or HAC. Run majorityVote method with a threshold t for joining communities (see description). Run HAC nethod for HAC meta clustering with a distance threshold t and distance type distance arguments, as permited by scipys scipy.cluster.hierarchy.linkage

#### BRIM. unknown number of communities, unsigned
```python
import networkx as nx
import pyBRIM

g = nx.davis_southern_women_graph()
bs = pyBRIM.BRIM_solver(g)
R_t,S, Q_max = bs.BRIM_bisec()
nodes_to_communities = bs.translate_communities(R_t,S)
```

#### BRIM. known number of communities, unsigned
```python
import networkx as nx
import pyBRIM

g = nx.davis_southern_women_graph()
bs = pyBRIM.BRIM_solver(g)
R_t,S, Q_max = bs.fit_transform(4)
nodes_to_communities = bs.translate_communities(R_t,S)
```
#### BRIM. unknown number of communities, signed

#### Meta clustering: HAC

#### Meta clustering: majority voting

#### Meta clustering: best meta clusters

# Description

# References

