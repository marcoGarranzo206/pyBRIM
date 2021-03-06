# pyBRIM
Algorithms for community detection on bipartite graphs based on BRIM algorithm

# Overview

Finding communities is a very common and important step in most network analysis. Network communities can help find groups of similar nodes, such as users in a social network, groups of products and buyers from sales data,  proteins that belong to similar pathways in a metabolic network, or groups of drugs in a drug interaction network.

Many community finding algorithms are based on Newmans modularity optimization, initially design for "normal" networks (unweighted, unsigned). Here I implement algorithms for community finding in python for specifically bipartite networks, signed or unsigned, based on the BRIM algorithm. There is additionally a meta clustering module for combining results from different algorithm runs.

# Dependencies

numpy for array operations<br>
scipy for hierarchical clustering<br>
networkx for several graph methods. graphs have to be networkx graphs.<br>
tdqm for progress bar in meta module

# Usages and examples:

There are two modules, BRIM and meta.

## BRIM
BRIM module revolves around the BRIM_solver object.Pass it the graph you want to use, the null model (ie signed or unsigned networks) and the resolution value (parameter that determines community sizes, more in description).

After creation you can run the fit_transform method(c) to find c communities. If you dont know the number of communities, you can use the BRIM_bisec method, which you can also pass it a maximum number of communities to search for. The algorithm relies on random initializations, meaning that different runs may (most probably) contain different communities. To solve this issue there is the meta module.

Both methods return a dictionary of node to community assingments and the Q value.

## meta
meta module revolves around the metaCluster object. You pass it an array of cluster runs. Each row corresponds to a cluster run, each column a node. (i,j) element s the community of node j in cluster run i. A co-occurrence matrix is built.

After creation you can do meta clustering based on majority voting or HAC. Run majorityVote method with a threshold t for joining communities. Run HAC nethod for HAC meta clustering with a distance threshold t and distance type distance arguments, as permited by scipys scipy.cluster.hierarchy.linkage

A problem now arises: which is the best meta algorithm. That is up to you to decide. One choice is to use the meta algorithm which yields the best average normalized mutual information score w.r.t to the original cluster runs, since it is on average the most "consistent" with all the different runs. If that is your metric, you may try to optimize directly with general discrete solvers such as genetic algorithms or simulated annealing, but in my experience the meta algorithms work quite well at a fraction of the time.

## Examples
We will be using two networks: davis southern woman graph available at networkx and the US Senate - bills dataset available at the github repo DSE-MSU/signed-bipartite-networks

#### BRIM. unknown number of communities, unsigned
```python
import networkx as nx
import pyBRIM

g = nx.davis_southern_women_graph()
bs = pyBRIM.BRIM_solver(g)
nodes_to_communities, Q_max = bs.BRIM_bisec()
```

#### BRIM. known number of communities, unsigned
```python
import networkx as nx
import pyBRIM

g = nx.davis_southern_women_graph()
bs = pyBRIM.BRIM_solver(g)
nodes_to_communities, Q_max = bs.fit_transform(4)
```
#### BRIM. unknown number of communities, signed, setting resolution value to 2
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
    
bs = pyBRIM.BRIM_solver(g, "neg", 2)
nodes_to_communities, Q_max = bs.BRIM_bisec()
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
    
    
    nodes_to_communities, Q_max = bs.BRIM_bisec()
    for node,cluster in nodes_to_communities.items():
        
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
    
    
    nodes_to_communities, Q_max = bs.BRIM_bisec()
    for node,cluster in nodes_to_communities.items():
        
        cluster_runs[i, nodes_to_idxs[node]] = cluster
        
mc = pyBRIM.metaCluster(cluster_runs)
clusters = mc.majorityVote(.9)#array. ith position is community of node with index i```
```
# Modularity description

Network communities can be loosely defined as groups of nodes that interact a lot with each other and sparsely with other nodes. A rigorous mathematical definition it is not, and Newmann defined the quantity of modularity to judge community assingments of a network based on that principle [1]:


<img src="https://render.githubusercontent.com/render/math?math=Q = \frac{1}{2m}\sum_{i,j}[A_{ij} - \frac{k_ik_j}{2m}]\delta(c_i,c_j)">

Essentially it goes over all pairs of nodes i,j in the graph, and if they are in the same community, sums up the adjacency matrix value (1 if connected, 0 otherwise) and subtracts it the expected value it should have. Here we use the configuration null model, where the expected value of Aij is proprotional to the products of the degrees of i and j.

### Bipartite networks
The problem when using this modularity for bipartite networks is that by defintion, nodes within a set do not interact. If placed on the same community, Aij will be 0, but the expected value will be some positive number under the configuration null model, bringing the modularity down. 

The solution proposed by Barber is simple, have the expected adjacency value of two nodes in the same set be 0. With this in mind and with some clever linear algebra, an iterative algorithm was formed to find communities in bitartite networks [2].

### Signed networks
Signed networks can have positive (Aij = 1) or negative (Aij = -1) interactions between members. These can include:

- Social networks: people can like or dislike each other
- drug-protein interaction network: drugs can inhibit or increase protein activity

Newmans modularity breaks down in these cases. The underlying reason is that it has a probabilistic formulation which cannot function with negative values in the expectation of Aij or in Aij. Gomez et al solved this problem by splitting modularity in two: one for the positive edges and one for the negative edges, with the ultimate goal to maximize Q as [3]:

<img src="https://render.githubusercontent.com/render/math?math=Q^+ = \frac{1}{2w^+}\sum_{i,j}[w_{ij}^+ - \frac{w_i^+w_j^+}{2w^+}]\delta(c_i,c_j)">
<img src="https://render.githubusercontent.com/render/math?math=Q^- = \frac{1}{2w^-}\sum_{i,j}[w_{ij}^- - \frac{w_i^-w_j^-}{2w^-}]\delta(c_i,c_j)">
<img src="https://render.githubusercontent.com/render/math?math=Q = \frac{2w^+}{2w^+ + 2w^-}Q^+ - \frac{2w^-}{2w^+ + 2w^-}Q^-">

With the w+ indicating values for the positive partition of the graph and viceversa. We can optimize Q using the same BRIM algorithm as before.

### Resolution value
Modularity suffers from a resolution limit. As networks get bigger and bigger, the limit to the community sizes it can detect increases. In big sparse networks, any edge will have low probability of ocurring, so modularity will favor much more putting two nodes together in the same community as long as they are connected [4].

This can be mitigated by putting a resolution parameter alpha > 1 to penalize two nodes in the same comunity:
<img src="https://render.githubusercontent.com/render/math?math=Q^+ = \frac{1}{2w^+}\sum_{i,j}[w_{ij}^+ - \alpha\frac{w_i^+w_j^+}{2w^+}]\delta(c_i,c_j)">

For negative networks, once a alpha is chosen, here we put 1/alpha as the resolution for Q-. We can optimize Q using the same BRIM algorithm as before.


# References

[1] Newman M. E. (2006). Modularity and community structure in networks. Proceedings of the National Academy of Sciences of the United States of America, 103(23), 8577???8582. https://doi.org/10.1073/pnas.0601602103

[2] Barber MJ. Modularity and community detection in bipartite networks. Phys Rev E Stat Nonlin Soft Matter Phys. 2007 Dec;76(6 Pt 2):066102. doi: 10.1103/PhysRevE.76.066102. Epub 2007 Dec 7. PMID: 18233893.

[3] Gomez, Sergio & Jensen, Pablo & Arenas, Alex. (2009). Analysis of community structure in networks of correlated data. Physical review. E, Statistical, nonlinear, and soft matter physics. 80. 016114. 10.1103/PhysRevE.80.016114. 

[4] Fortunato, S., & Barth??lemy, M. (2007). Resolution limit in community detection. Proceedings of the National Academy of Sciences of the United States of America, 104(1), 36???41. https://doi.org/10.1073/pnas.0605965104
