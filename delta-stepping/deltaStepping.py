"""
Author: Alexander Holmquist

File: deltaStepping.py

Purpose: sequential implementation of the delta-stepping algorithm.

We do not have to implement something efficient, this is just for
testing purposes.

The algorithm mimics the original idea of delta-stepping, from the
paper:
    U. Meyer, P. Sanders,
    Δ-stepping: a parallelizable shortest path algorithm,
    Journal of Algorithms,
    Volume 49, Issue 1,
    2003,
    Pages 114-152,
    ISSN 0196-6774,
    https://doi.org/10.1016/S0196-6774(03)00076-2
"""

from graphUtils import *

from math import ceil, floor

class deltaStepping:
    def __init__(self, G, sourceNode, delta):
        # Our reusable definition of infinity.
        self.infinity = float("inf")

        # Inputs
        self.G = G
        self.sourceNode = sourceNode
        self.delta = delta

        # Weight position in the edge tuple format given by networkx
        self.wPos = 2

        # Intiialize all kinds of edges at once. The light and heavy
        # edges are simply indices to be used with self.edges
        self.edges, self.edgeWeights, self.lightEdges,\
            self.heavyEdges = self.initEdges()

        # B is our list of buckets
        self.B = self.initBuckets()
        
        # Tentative distance from sourceNode
        self.tent = {vertex: self.infinity for vertex in G.nodes()}

        # Start with only the source node, in the first bucket.
        self.relax(sourceNode, 0)

        # Algorithm results
        result = None

    # This runs the algorithm itself. See Section 2 of the paper.
    def run(self):
        while not self.isBEmpty():
            # First non-empty bucket
            i = self.getMinBucket()
            # R is the set of nodes we remove from the min bucket.
            R = set([])
            
            # Light edges have to be relaxed iteratively, for the
            # relaxation may insert neighbor nodes in the min
            # bucket.
            while self.B[i] != []:
                Req = self.findRequests(self.B[i], 0) # 0 == light
                R = R.union(self.B[i])
                self.B[i] = []
                self.relaxRequests(Req)
                
            # Relax heavy edges all at once, for they cannot be
            # reinserted in the min bucket by the relaxations.
            Req = self.findRequests(R, 1) # 1 == heavy
            self.relaxRequests(Req)

            self.recycleBuckets()

        self.result = self.tent

    ##################################################
    # Algorithm auxiliary procedures
    ##################################################
    def findRequests(self, bucket, relaxType):
        if relaxType == 0:
            return self.findRequestsAux(bucket, self.lightEdges)
        elif relaxType == 1:
            return self.findRequestsAux(bucket, self.heavyEdges)
        else:
            raise ValueError("relaxType should be either light or \
            heavy")

    # This is very, very inneficient, and breaks the
    # work-efficiency of the algorithm (I think it makes it
    # quadratic). However, this implementation is here just to
    # show the working of the algorithm, and its correctness.
    def findRequestsAux(self, bucket, edges):
        Req = []
        for v in bucket:
            for idx in edges:
                w = self.edges[idx][1]
                if (v, w) == self.edges[idx]:
                    Req.append((w, self.tent[v] + 
                                self.edgeWeights[idx]))
    
        return Req

    def relaxRequests(self, Req):
        for (node, dist) in Req:
            self.relax(node, dist)

    def relax(self, node, newDist):
        if self.tent[node] != self.infinity:
            if newDist < self.tent[node]:
                self.B[floor(self.tent[node] / self.delta)]\
                    .remove(node)
                self.B[floor(newDist / self.delta)].append(node)
                self.tent[node] = newDist
        else:
            self.B[floor(newDist / self.delta)].append(node)
            self.tent[node] = newDist

    # Checks if B is empty. This is inneficient.
    def isBEmpty(self):
        for i in list(self.B.keys()):
            if self.B[i] != []:
                return False
        return True

    # Fetches the minimum bucket in B. This is inneficient.
    def getMinBucket(self):
        for i in list(self.B.keys()):
            if self.B[i] != []:
                break

        return i

    # Cleans all the empty buckets, until we find a non-empty one.
    def recycleBuckets(self):
        bKeys = list(self.B.keys())
        bot = bKeys[0]
        top = bKeys[-1:][0]

        for idx in bKeys:
            if self.B[idx] == []:
                # Remove the empty list from the bottom
                self.B.pop(idx)
                # Add a new empty list at the top
                self.B[top + idx - bot + 1] = []
            else:
                break

    ##################################################
    # INITIALIZATION ROUTINES
    ##################################################

    # Checks if an edge is heavy
    def isHeavy(self, edgeWeights, idx):
        return edgeWeights[idx] > self.delta

    def initEdges(self):
        gEdges = self.G.edges.data("weight") # edges from G
        selfEdges = [e[:2] for e in gEdges] # goes to self.edges
        selfEdgeWeights = [e[2] for e in gEdges] # self.edgeWeights

        # Now for the light/heavy indices. This associates the
        # light/heavy information to the edges. We only fetch the
        # indices of the edges that are light / heavy, to keep
        # things uniform, and save a little space.
        lightEdges = []
        heavyEdges = []
        for idx in range(len(selfEdges)):
            if self.isHeavy(selfEdgeWeights, idx):
                heavyEdges.append(idx)
            else:
                lightEdges.append(idx)

        return selfEdges, selfEdgeWeights, lightEdges, heavyEdges

    def initBuckets(self):
        # Fetching the max is a linear operation
        numBuckets = ceil(max([w / self.delta for w in 
                               self.edgeWeights])) + 1

        # Intialize with dummy lists
        buckets = dict()
        for ind in range(numBuckets):
            buckets[ind] = []

        return buckets

def main():
    # Build the entry graph
    inGraphFilePath = "test/binTree.csv"
    outPath = "graph.png"
    gb = graphBuilder(inGraphFilePath)
    G = gb.buildGraphFromFile()
    gb.drawGraph()
    plt.savefig(outPath)
    
    # This can be used to get parameters as input from the user.
    # sourceNode = int(input("Source node: "))
    # delta = float(input("Delta: "))

    sourceNode = 0
    delta = 1000

    alg = deltaStepping(G, sourceNode, delta)
    alg.run()
    print("Results: ", alg.result)


    
if __name__ == '__main__':
    main()
