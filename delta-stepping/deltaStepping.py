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
        self.infinity = float("inf")

        # Inputs
        self.G = G
        self.sourceNode = sourceNode
        self.delta = delta

        # Things reused all over the place
        # 
        # Weight position in the edge tuple format given by networkx
        self.wPos = 2
        self.edges, self.edgeWeights = self.initEdges()

        # These are simply indices
        self.lightEdges, self.heavyEdges = self.getLHEdges()
        
        # Initialization
        #
        self.B = self.initBuckets()
        
        self.tent = {vertex: self.infinity for vertex in G.nodes()}
        self.relax(sourceNode, 0)

        # Algorithm results
        result = None

    # This runs the algorithm itself
    def run(self):
        while not self.isBEmpty():
            i = self.getMinBucket()

            R = set([])
            
            while self.B[i] != []:
                Req = self.findRequests(self.B[i], 0) # 0 == light
                R = R.union(self.B[i])
                self.B[i] = []
                self.relaxRequests(Req)
                
            # Relaxing heavy edges
            Req = self.findRequests(R, 1) # 1 == heavy
            self.relaxRequests(Req)

            self.recycleBuckets()

        self.result = self.tent
    ##################################################
    # Algorithm auxiliary procedures
    ##################################################
    def findRequests(self, bucket, relaxType):
        if relaxType == 0:
            return self.findRequestsLight(bucket)
        elif relaxType == 1:
            return self.findRequestsHeavy(bucket)
        else:
            raise ValueError("relaxType should be either light or \
            heavy")

    def findRequestsLight(self, bucket):
        Req = []
        for v in bucket:
            for idx in self.lightEdges:
                w = self.edges[idx][1]
                if (v, w) == self.edges[idx]:
                    Req.append((w, self.tent[v] + 
                                self.edgeWeights[idx]))
    
        return Req

    def findRequestsHeavy(self, bucket):
        Req = []
        for v in bucket:
            for idx in self.heavyEdges:
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

    # Checks if B is empty
    def isBEmpty(self):
        for i in list(self.B.keys()):
            if self.B[i] != []:
                return False
        return True

    # Fetches the minimum bucket in B
    # TODO: this is very inneficient
    def getMinBucket(self):
        for i in list(self.B.keys()):
            if self.B[i] != []:
                break

        return i

    # Cleans all the empty buckets, until we find a non-empty one
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

    # This associates the light/heavy information to the edges.
    #  0 stands for light, and 1 stands for heavy.        
    def getLHEdges(self):
        lightEdges = []
        heavyEdges = []
        for idx in range(len(self.edges)):
            if self.isHeavy(idx):
                heavyEdges.append(idx)
            else:
                lightEdges.append(idx)

        return lightEdges, heavyEdges

    # Checks if a edge is heavy
    def isHeavy(self, idx):
        return self.edgeWeights[idx] > self.delta

    def initEdges(self):
        gEdges = self.G.edges.data("weight")
        return ([e[:2] for e in gEdges], [e[2] for e in gEdges])

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
    inGraphFilePath = "test/binTree.csv"
    gb = graphBuilder(inGraphFilePath)
    G = gb.buildGraphFromFile()
    gb.drawGraph()
    outPath = "graph.png"
    plt.savefig(outPath)
    
    # sourceNode = int(input("Source node: "))
    # delta = float(input("Delta: "))

    sourceNode = 0
    delta = 1000

    alg = deltaStepping(G, sourceNode, delta)
    alg.run()
    print("Results: ", alg.result)


    
if __name__ == '__main__':
    main()
