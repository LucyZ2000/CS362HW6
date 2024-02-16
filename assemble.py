import sys

class DeBruijnGraph:

    @staticmethod
    def split_into_kmers(str, k):
        '''split a string into segments of length k'''
        kmers = []
        for i in range(len(str) - k + 1):
            kmers.append(str[i:i+k])
        return kmers

    @staticmethod
    def convert_to_list(filename):
        '''convert a textfile into a list of strings'''
        # I wrote this part of the code for testing purposes
        if isinstance(filename, list):
            return filename
        stringList = []
        file = open(filename, "r")
        for line in file:
            if line[-1] == "\n":
                line = line[:-1]
            stringList.append(line)
        return stringList

    class Node:
        def __init__(self, kmer):
            self.kmer = kmer # a string
            self.outdegree = 0
            self.indegree = 0

    def __init__(self, filename, k):
        '''Create the deBruijn Graph object with the set of strings
        in the specified file using parameter k'''
        stringList = self.convert_to_list(filename)

        self.outedges = {} # {node1 : [node1, node2, node3], node2 : []}

        self.inedges = {} # {node1 : [node1, node2, node3], node2 : []}

        self.outedges2 = {}

        # a dictionary with kmers(string) as keys and node objects as values

        self.vertices = {}
        for str in stringList:
            for kmer in self.split_into_kmers(str, k):
                kmerL, kmerR = kmer[:-1], kmer[1:]
                if kmerL in self.vertices:
                    # set node to be the existing node that contains the same kmer
                    nodeL =  self.vertices[kmerL]
                else:
                    # create a new node
                    self.vertices[kmerL] = self.Node(kmerL)
                    nodeL = self.vertices[kmerL]

                if kmerR in self.vertices:
                    nodeR = self.vertices[kmerR]

                else:

                    self.vertices[kmerR] = self.Node(kmerR)
                    nodeR = self.vertices[kmerR]
                    
                nodeL.outdegree += 1
                nodeR.indegree += 1

                self.outedges.setdefault(nodeL, []).append(nodeR)
                self.inedges.setdefault(nodeR, []).append(nodeL)
                self.outedges2.setdefault(nodeL, []).append(nodeR)

    def fix_edges(self):
        '''makes sure the first kmer and last kmer are in self.outedges and self.inedges'''
        for node in self.vertices.values():
            if node not in self.outedges:
                self.outedges[node] = []
            if node not in self.inedges:
                self.inedges[node] = []
    
    def get_start_vertices(self):
        '''returns a list of vertices that start the eulerian paths'''
        startVertices = []
        for node in self.vertices.values():
            if node.outdegree > node.indegree:
                startVertices.append(node)
    
        return startVertices

    def find_eulerian_path(self, vertex):
        ''''''
        countdict = {} # countdict keeps track of the edges that are already visited
        contiglist = []
        contig = vertex.kmer
        for i in range(len(self.outedges[vertex])):

            current = self.outedges[vertex][i]

            countdict.setdefault(current, 0)

            while countdict[current] < len(self.outedges[current]):

                # while the current node still has available edges (edges that are not yet
                # visited) 

                next = self.outedges[current][countdict[current]]

                countdict.setdefault(next, 0)

                contig += current.kmer[-1]

                countdict[current] += 1 

                current = next

                self.outedges.setdefault(current, [])
            
            countdict = {} # initialize the coundict for every eulerian path
            
            contiglist.append(contig)

        return contiglist


    def find_all_eulerian_paths(self):
        '''returns a list of all contigs 
        collected from eulerian paths'''
        contigs = []
        vertices = self.get_start_vertices()
        for vertex in vertices:
            contigs += self.find_eulerian_path(vertex)
        contigs = sorted(contigs, key = len)
        contigs.reverse()
        return contigs[0:10]
    
    def n50l50(self):
        '''returns the N50 and L50 values'''
        contigs = self.find_all_eulerian_paths()
        sortedContigs = sorted(contigs, key=len)
        sortedContigs.reverse()
        totalLength = 0
        n = len(sortedContigs)
        for i in range(n):
            totalLength += len(sortedContigs[i])

        k = 0
        j = 0
        while k <= (totalLength // 2):
            k += len(sortedContigs[j])
            j += 1
    
        return len(sortedContigs[j]), n - j


def write_file(filename, list):
    '''outputs each string(simulated read) in the list as a separate line
    in a file'''
    file = open(filename, "w")
    for string in list:
        file.write(string + "\n")
        
    file.close()


def main():
    args = sys.argv[1:]
    filename = args[0]
    k = int(args[1])
    graph = DeBruijnGraph(filename, k)
    graph.fix_edges()
    n50, l50 = graph.n50l50()
    print("N50: ", n50)
    print("L50: ", l50)

if __name__ == "__main__":
    main()