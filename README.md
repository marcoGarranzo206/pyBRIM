# pyBRIM
Algorithms for community detection on bipartite graphs based on BRIM algorithm

# Overview

Finding communities is a very common and important step in most network analysis. Many community finding algorithms are based on Newmans modularity optimization, initially design for "normal" networks (unweighted, unsigned). Here I implement algorithms for community finding in python for specifically bipartite networks, signed or unsigned, based on the BRIM algorithm. There is additionally a meta clustering module for combining results from different algorithm runs.

# Dependencies

numpy for array operations<br>
scipy for hierarchical clustering<br>
networkx for several graph methods. graphs have to be networkx graphs.<br>
tdqm for progress bar in meta module

# Usages and examples:

There are two modules, BRIM and meta.

## BRIM
BRIM module revolves around the BRIM_solver object.Pass it the graph you want to use, the null model (for signed or unsigned networks) and the resolution value (to determine community sizes, more in description).

After creation you can run the fit_transform method to find communities when the number of communities is found. If you dont know, you can use the BRIM_bisec method, which you can pass it a maximum number of communities.

## meta
meta module revolves around the metaCluster object. You pass it an array of cluster runs. Each row corresponds to a cluster run, each column a node. (i,j) element s the community of node j in cluster run i. A co-occurrence matrix is built.

After creation you can do meta clustering based on majority voting or HAC. Run majorityVote method with a threshold t for joining communities (see description). Run HAC nethod for HAC meta clustering with a distance threshold t and distance type distance arguments, as permited by scipys scipy.cluster.hierarchy.linkage

# Description

# Examples

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
```python
import requests
import networkx as nx
import pyBRIM

url = 'https://raw.githubusercontent.com/DSE-MSU/signed-bipartite-networks/master/data/senate1to10_cikm2019_balance_in_signed_bipartite_networks.txt'
req = requests.get(url)
g = nx.Graph()
for i, (u,v,w) in enumerate(map(lambda x: x.split(), req.text.split("\n"))):
    
    if i == 0:
        
        continue
        
    w = int(w)
    u = "u_" + u
    v = "v_" + v
    g.add_edge(u,v, weight = w)
    
bs = pyBRIM.BRIM_solver(g, "neg")
R_t,S, Q_max = bs.BRIM_bisec()

nodes_to_communities = bs.translate_communities(R_t,S).items()
```

#### Meta clustering: HAC
```python
import requests
import networkx as nx
import pyBRIM
import numpy as np

url = 'https://raw.githubusercontent.com/DSE-MSU/signed-bipartite-networks/master/data/senate1to10_cikm2019_balance_in_signed_bipartite_networks.txt'
req = requests.get(url)
g = nx.Graph()
for i, (u,v,w) in enumerate(map(lambda x: x.split(), req.text.split("\n"))):
    
    if i == 0:
        
        continue
        
    w = int(w)
    u = "u_" + u
    v = "v_" + v
    g.add_edge(u,v, weight = w)
    
bs = pyBRIM.BRIM_solver(g, "neg")
nodes_to_idxs = {n:i for (i,n) in enumerate(g)} #map each node to the same index
n_cluster_runs = 50
cluster_runs = np.zeros((n_cluster_runs, len(g) ))

#run multiple runs of the algorithm
for i in range(n_cluster_runs):
    
    
    R_t,S, Q_max = bs.BRIM_bisec()
    for node,cluster in bs.translate_communities(R_t,S).items():
        
        cluster_runs[i, nodes_to_idxs[node]] = cluster
        
mc = pyBRIM.metaCluster(cluster_runs)
clusters = mc.HAC(0.1, "median")#array. ith position is community of node with index i
```

#### Meta clustering: majority voting
```python
import requests
import networkx as nx
import pyBRIM
import numpy as np

url = 'https://raw.githubusercontent.com/DSE-MSU/signed-bipartite-networks/master/data/senate1to10_cikm2019_balance_in_signed_bipartite_networks.txt'
req = requests.get(url)
g = nx.Graph()
for i, (u,v,w) in enumerate(map(lambda x: x.split(), req.text.split("\n"))):
    
    if i == 0:
        
        continue
        
    w = int(w)
    u = "u_" + u
    v = "v_" + v
    g.add_edge(u,v, weight = w)
    
bs = pyBRIM.BRIM_solver(g, "neg")
nodes_to_idxs = {n:i for (i,n) in enumerate(g)} #map each node to the same index
n_cluster_runs = 50
cluster_runs = np.zeros((n_cluster_runs, len(g) ))

#run multiple runs of the algorithm
for i in range(n_cluster_runs):
    
    
    R_t,S, Q_max = bs.BRIM_bisec()
    for node,cluster in bs.translate_communities(R_t,S).items():
        
        cluster_runs[i, nodes_to_idxs[node]] = cluster
        
mc = pyBRIM.metaCluster(cluster_runs)
clusters = mc.majorityVote(.9)#array. ith position is community of node with index i```
```

# References

