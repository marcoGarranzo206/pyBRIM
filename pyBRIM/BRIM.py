import networkx as nx
import numpy as np

def Q(R_T, B, S,m):

    return np.trace((R_T @ B @ S))/m

def BRIM_loop(B,c, S, m):

    """
    Given modularity matrix B (or a rank approximation of it) and a number of modules c,
    run an iteration of BRIM algorithm described in

    There are two sets (bipartite graph)
    Assign modules to the first set (say A) randomly.
    Based on those assingments, assign modules to the other set (say B)
    Based on this new B assignment, reassign the first set to new modules
    """

    Q_old = -1
    while True:

        #assign to R_T based on S
        T = B @ S
        R_T = np.zeros((c, B.shape[0]))
        R_modules = np.argmax(T, axis=1).reshape(1,-1)
        R_T[R_modules , range(B.shape[0])] = 1
        #assign to S based on R_T

        T = R_T @ B
        S = np.zeros((B.shape[1], c))
        S_modules = np.argmax(T, axis=0)
        S[range(B.shape[1]), S_modules ] = 1
        Q_new = Q(R_T,B,S,m)

        if Q_new > Q_old:

            Q_old = Q_new

        else:

            break

    return R_T, S, Q_new


class BRIM_solver:

    def __init__(self, g, null = "config",res = 1):

        self.null = null
        self.resolution = res
        self.g = g

        if null in ("config", "configuration", "configPos", "configuartionPos"):

            self._BNullBipartiteConfigPos()

        elif ("configNeg", "configurationNeg"):

            self._BNullBipartiteConfigNeg()


    def fit_transform(self, c, assingments = None):


        S = np.zeros((self.B.shape[1], c))
        if assingments is None:

            S[range(self.B.shape[1]), np.random.choice(c, self.B.shape[1])] = 1

        else:

            S[range(self.B.shape[1]), assingments] = 1

        return BRIM_loop(self.B,c,S,self.m)

    def translate_communities(self,R_t,S):

        """
        Recieve communty assingments for R_t and S matrix and output
        a dictionary of node-communities
        >>> g = nx.davis_southern_women_graph()
        >>> bs = BRIM_solver(g)
        >>> R_t,S,Q = bs.fit_transform(3) # also works with bs.BRIM_bisec()
        >>> nodes_communites = bs.translate_communities(R_t,S)
        """

        node_to_communites = {}
        #R_t: matrix of (number of communites, len(s2)). each column (node) has a one value
        #in the row corresponding to its communities

        node_to_communites.update({n:c for (n,c) in zip(self.s2, np.argmax(R_t,0))})

        #S: matrix of  (len(s1),number of communites). each row (node) has a one value
        #in the column corresponding to its communities

        node_to_communites.update({n:c for (n,c) in zip(self.s1, np.argmax(S,1))})
        return node_to_communites

    def BRIM_bisec(self, c_max = None):

        """
        Find modules using brim alg without a fixed number of c, using binary search method described
        in BRIM paper (not exactly for now)

        First stage:

            start with low number of communities (2), calculate Q
            double number of communities until Q decreases

        Second stage
            bisection to find optimal c:

                interpolation between number of communities  that caused
                decrease in Q and previous number of communities.

                repeat until convergence
        """

        if c_max is None:

            c_max = self.B.shape[1]
        c = 2
        half = self.B.shape[1]//2
        assingments_s = np.ones(self.B.shape[1]).astype(int)
        Q_old = -100
        prev_c = 1

        while True:


            assingments_s[np.random.choice(range(self.B.shape[1]), half, replace = False )] = np.random.randint(low = 0, high = c, size = half)
            R_T, S, Q_new = self.fit_transform(c,assingments = assingments_s)
            assingments_s = np.where(S != 0 )[1]
            if Q_new > Q_old:

                R_T_max, S_max, Q_max = R_T, S, Q_new
                prev_c = c
                c = min(2*c, c_max)
                Q_old = Q_new

            else:

                break

            if c >= (c_max):
                break

        c_prev = c
        high = c #upper bound of number of communities
        low = prev_c #lower bound of number of communities
        c = (high + low)//2
        Q_upper = Q_new # Q value for number of communities that is the upper bound of number of communities
        Q_lower = Q_old # Q value for number of communities that is the lower bound of number of communities
        #confusingly enough, Q_lower > Q_upper (at least true in first loop, and assumed to be true always)

        while True:

            if high < low or c_prev == c:

                return R_T_max,S_max, Q_max

            assingments_s[np.random.choice(range(self.B.shape[1]), half, replace = False )] = np.random.randint(low = 0, high = c, size = half)
            #cluster number in previus run might be higher than total number of clusters now
            # causing issues when indexing. Assing each cluster to a distinct number
            #between 0 and c
            clusters_to_idxs = {n:min(i,c) for (i,n) in enumerate(set(assingments_s))}
            assingments_s = np.array([ clusters_to_idxs[n] for n in assingments_s])
           
            #print(c, max(assingments_s))
            R_T, S, Q_new = self.fit_transform(c + 1,assingments = assingments_s)
            assingments_s = np.where(S != 0 )[1]


            if Q_new > Q_max:

                R_T_max, S_max, Q_max = R_T, S, Q_new

            if  Q_new < Q_lower:

                #if new Q value is lower than Q value of the lower bound on c
                #search between new c and lower bound

                Q_upper = Q_new
                high = c

            else:

                #if new Q value is higher or equal than Q value of the lower bound on c
                #search between new c and upper bound
                Q_lower = Q_new
                low = c

            c_prev = c
            c = (high + low)//2
    
    def _BNullBipartiteConfigPos(self):

        """
        Computes B_tilde of a bipartite graph following a null configuration model
        Prob of edge (i,j) is 0 if they belong to the same set, degree_i*degree_j/m otherwise
        will end up with a square matrix of two diagonal 0 blocks, and 2 non diagonal blocks
        who are each others transpose

        As such, save space computing only one of the non 0 blocks
        """
        seta, setb = nx.bipartite.sets(self.g)

        if len(seta) > len(setb):

            self.s1 = seta
            self.s2 = setb

        else:

            self.s1 = setb
            self.s2 = seta


        A = nx.bipartite.biadjacency_matrix(self.g,row_order= self.s2, column_order=self.s1)
        ki = np.sum(A, axis=1)
        dj = np.sum(A, axis = 0)
        self.m = np.sum(A)
        self.B = A - self.resolution*((ki@dj)/self.m)

    def _BNullBipartiteConfigNeg(self):

        """
        Computes B_tilde of a bipartite graph following a null configuration model that can include
        both pos and neg edges as described by gomez et al in
        Analysis of community structure in networks of correlated data (2009)

        Prob of any edge (i,j) is 0 if they belong to the same set
        Prob of pos edge (i,j) is pos_degree_i*pos_degree_j/w_pos
        Prob of neg edge (i,j) is neg_degree_i*neg_degree_j/w_pos

        B = A - prob_pos + prob_neg

        will end up with a square matrix of two diagonal 0 blocks, and 2 non diagonal blocks
        who are each others transpose

        As such, save space computing only one of the non 0 blocks
        """

        seta, setb = nx.bipartite.sets(self.g)


        if len(seta) > len(setb):

            self.s1 = seta
            self.s2 = setb

        else:

            self.s1 = setb
            self.s2 = seta

        A = nx.bipartite.biadjacency_matrix(self.g,row_order = self.s2, column_order = self.s1).toarray()
        A_pos = np.zeros(A.shape)
        A_pos[A > 0] = 1

        A_neg = np.zeros(A.shape)
        A_neg[A < 0] = 1

        ki_pos = np.sum(A_pos, axis=0).reshape(-1,1)
        dj_pos = np.sum(A_pos, axis = 1).reshape(1,-1)

        ki_neg = np.sum(A_neg, axis=0).reshape(-1,1)
        dj_neg = np.sum(A_neg, axis = 1).reshape(1,-1)


        self.w_pos = np.sum(A_pos)
        self.w_neg = np.sum(A_neg)
        self.m = self.w_pos + self.w_neg

        self.B = A - self.resolution*((ki_pos@dj_pos).T/self.w_pos) + (1/self.resolution)*((ki_neg@dj_neg).T/self.w_neg)

