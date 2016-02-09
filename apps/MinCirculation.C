// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#define WEIGHTED 1
#define DEBUG
#include <vector>

using namespace std;

#include "ligra.h"

struct PR_F {
};

struct PR_Vertex_F {
    int* Visited;

    PR_Vertex_F(int* _Visited) : Visited(_Visited) {
    }

    inline bool operator()(uintE i) {
        Visited[i] = 0;
        return 1;
    }
};

bool isOutAdmissible(graph& GA, uintT node, uintT outEdge) {
    uintT neighbor = GA.V[node].getOutNeighbor(outEdge);
    intE c_p = GA.getOutInfo(node,outEdge).weight + GA.p[neighbor] - GA.p[node];
    intE residual = GA.getOutInfo(node,outEdge).capacity - GA.getOutInfo(node,outEdge).flow;
    cout << node << "->" << neighbor << " admissible: " << (c_p < 0 && residual > 0)
            << "(c_p = " << c_p << ", residual: " << residual << ")\n";
    return (c_p < 0 && residual > 0);
}
bool isInverseAdmissible(graph& GA, uintT node, uintT inEdge) {
    uintT neighbor = GA.V[node].getInNeighbor(inEdge);
    intE c_p = -GA.getInInfo(node,inEdge).weight + GA.p[neighbor] - GA.p[node];
    intE residual = GA.getInInfo(node,inEdge).flow - GA.getInInfo(node,inEdge).lower;
    cout << node << "<-" << neighbor << " admissible: " << (c_p < 0 && residual > 0)
        << "(c_p = " << c_p << ", residual: " << residual << ")\n";
    return (c_p < 0 && residual > 0);
}

void dfs(graph& GA, uintT* blockingFlow, uintT nodeId) {
    cout << "  DFS " << nodeId << ":\n";
    for (uintT i = 0; i < GA.V[nodeId].outDegree; i++) {
        //for each out edge check if it is admissible
        //if so - update distance and run dfs
        uintT neighbor = GA.V[nodeId].getOutNeighbor(i);
        if (isOutAdmissible(GA,nodeId,i) && blockingFlow[neighbor] == INT_MAX) {
            blockingFlow[neighbor] = nodeId;
            dfs(GA, blockingFlow, neighbor);
        }
    }
    //the same for inverse edges
    for (uintT i = 0; i < GA.V[nodeId].inDegree; i++) {
        //for each out edge check if it is admissible
        //if so - update distance and run dfs
        uintT neighbor = GA.V[nodeId].getInNeighbor(i);
        if (isInverseAdmissible(GA,nodeId,i) && blockingFlow[neighbor] == INT_MAX) {
            blockingFlow[neighbor] = nodeId;
            dfs(GA, blockingFlow, neighbor);
        }
    }
}

void calculateExcesses(graph& GA, vertexSubset& excesses, vertexSubset& deficits) {
    cout << "Calculate excesses..." << endl;
    
    uintT totalExcesses = 0;
    uintT totalDeficits = 0;
    for (uintT i = 0; i < GA.n; i++)
    {
        excesses.d[i] = false;
        deficits.d[i] = false;
        intE flow_sum = 0;
        for (uintT j = 0; j < GA.V[i].outDegree; j++) {
            flow_sum -= GA.getOutInfo(i,j).flow;
        }
        for (uint j = 0; j < GA.V[i].inDegree; j++) {
            flow_sum += GA.getInInfo(i,j).flow;
        }
        flow_sum += GA.supply[i];
        if (flow_sum > 0) {
            excesses.d[i] = true;
            totalExcesses++;
        }
        if (flow_sum < 0) {
            deficits.d[i] = true;
            totalDeficits++;
        }
    }
    excesses.m = totalExcesses;
    deficits.m = totalDeficits;
    
    cout << "excesses: " ;
    for (uintT i = 0; i < GA.n; i++) {
        cout << i << ":" << excesses.d[i] << " ";
    }
    cout << endl;
    cout << "deficits: " ;
    for (uintT i = 0; i < GA.n; i++) {
        cout << i << ":" << deficits.d[i] << " ";
    }
    cout << endl;
}

void printGraph(graph GA) {
    for (uintE i = 0; i < GA.n; i++) {
        cout << "#" << i << " (" << GA.p[i] << ", supply: " << GA.supply[i] << "):\t\n";
        cout << "    " << "OUT:\n";
        for (uintE j = 0; j < GA.V[i].outDegree; j++) {
            cout << "        " << GA.V[i].getOutNeighbor(j) 
                    << " (w: " << GA.getOutInfo(i,j).weight << "; "
                    << " f: " << GA.getOutInfo(i,j).flow << "; "
                    << " c: " << GA.getOutInfo(i,j).capacity << "; "
                    << " l: " << GA.getOutInfo(i,j).lower << ")\n";
        }
        cout << "    " << "IN:\n";
        for (uintE j = 0; j < GA.V[i].inDegree; j++) {
            cout << "        " << GA.V[i].getInNeighbor(j) 
                    << " (w: " << GA.getInInfo(i,j).weight << "; "
                    << " f: " << GA.getInInfo(i,j).flow << "; "
                    << " c: " << GA.getInInfo(i,j).capacity << "; "
                    << " l: " << GA.getInInfo(i,j).lower << ")\n";
        }
    }
}

void refine(graph& GA, const double epsilon) {
    
    cout << "New iteration: e = " << epsilon << endl;
    printGraph(GA);
    
    //loop through arcs, saturate those with c_p less than 0 and add excesses
    for (uintT i = 0; i < GA.n; i++)
    {
        for (uintT j = 0; j < GA.V[i].getOutDegree(); j++) {
            uintT neighbour = GA.V[i].getOutNeighbor(j);
            intE c_p = GA.getOutInfo(i,j).weight + GA.p[neighbour] - GA.p[i];
            if (c_p < 0) {
                GA.getOutInfo(i,j).flow = GA.getOutInfo(i,j).capacity;
            }
            cout << "c_p("<< i << "," << neighbour << ") = " << c_p << endl;
        }
    }
    
    vertexSubset excesses(GA.n);
    vertexSubset deficits(GA.n);
    excesses.toDense();
    deficits.toDense();
    calculateExcesses(GA,excesses,deficits);
    
    while (excesses.isEmpty() == false) {
        
        cout << "Raise potentials..." << endl;
        
        //calculate shortest distances
        uintE* ShortestDistance = newA(uintT,GA.n);
        vertexSubset Frontier(GA.n);
        Frontier.toDense();
        for(uintT i = 0; i < GA.n; i++) Frontier.d[i] = excesses.d[i];
        Frontier.m = excesses.m;
        
        for (uintT i = 0; i < GA.n; i++) {
            ShortestDistance[i] = INT_MAX;
        }
        
        //calculate shortest path in Residual (!) graph
        uintE length_to_deficit;
        while (!Frontier.isEmpty()) {
            vertexSubset nextIterSubset(GA.n);
            nextIterSubset.toDense();
            for (uintT i = 0; i < GA.n; i++) {
                //if vertex in the active subset
                if (Frontier.d[i]) {
                    //if source - zero distance by definition
                    if (excesses.d[i]) {
                        ShortestDistance[i] = 0;
                    }
                    
                    //if a vertex is a deficit - terminate
                    if (deficits.d[i]) {
                        nextIterSubset.del();
                        length_to_deficit = ShortestDistance[i];
                        break;
                    }
                    
                    //iterate through residual edges to both directions!
                    for (uintT j = 0; j < GA.V[i].outDegree; j++) {
                        if (GA.getOutInfo(i,j).flow < GA.getOutInfo(i,j).capacity) {
                            uintT neighbour = GA.V[i].getOutNeighbor(j);
                            cout << "got new forward residual edge: " << i << "->" << neighbour << endl;
                            if (ShortestDistance[neighbour] > ShortestDistance[i] + 1) {
                                ShortestDistance[neighbour] = ShortestDistance[i] + 1;
                                if (!nextIterSubset.d[neighbour]) {
                                    nextIterSubset.d[neighbour] = true;
                                    nextIterSubset.m++;
                                }
                            }
                        }
                    }
                    for (uintT j = 0; j < GA.V[i].inDegree; j++) {
                        if (GA.getInInfo(i,j).flow > GA.getInInfo(i,j).lower) {
                            uintT neighbour = GA.V[i].getInNeighbor(j);
                            cout << "got new backward residual edge: " << i << "->" 
                                    << j << "[original <-]" << endl;
                            if (ShortestDistance[neighbour] > ShortestDistance[i] + 1) {
                                ShortestDistance[neighbour] = ShortestDistance[i] + 1;
                                if (!nextIterSubset.d[neighbour]) {
                                    nextIterSubset.d[neighbour] = true;
                                    nextIterSubset.m++;
                                }
                            }
                        }
                    }
                }
            }
            Frontier.del();
            Frontier = nextIterSubset;
        }
        cout << "Shortest path done. MaxLength = " << length_to_deficit << endl;
        for (int i = 0; i < GA.n; i++) {
            cout << i << ":" << ShortestDistance[i] << " ";
        }
        cout << endl;
        
        //raise potentials
        for (uintT i = 0; i < GA.n; i++) {
            if (ShortestDistance[i] < length_to_deficit) {
                GA.p[i] += (length_to_deficit - ShortestDistance[i])*epsilon;
            }
        }
        cout << "new potentials: " << endl;
        printGraph(GA);
        
        //calculate blocking flow via depth-first search in admissible graph
        // depth-first search: creating an array of pointers to parents
        // can be extended to the case of non-unit capacity by saving 
        // a pointer to a parent together with total flow
        //terminate when a deficit reached
        uintT* blockingSearch = newA(uintT, GA.n);
        for (uintT i = 0; i < GA.n; i++)    blockingSearch[i] = INT_MAX;
        for (uintT i = 0; i < GA.n; i++) {
            if (excesses.d[i]) {
                dfs(GA,blockingSearch,i);
            }
        }
        cout << "Blocking flow array: ";
        for (int i = 0; i < GA.n; i++) {
            cout << i << ":" << blockingSearch[i] << " ";
        }
        cout << endl;
        cout << "blocking flow done." << endl;
        
        cout << "Paths: \n";
        //raise flow: back propagation from each deficit
        for (uintT i = 0; i < GA.n; i++) {
            if (deficits.d[i] && blockingSearch[i] < INT_MAX) {
                uintT curNode = i;
                while (excesses.d[curNode] != true) {
                    cout << curNode << "<-";
                    curNode = blockingSearch[curNode];
                }
                cout << curNode << endl;
            }
        }
        
        //now raise flows
        cout << "Raising flows... Path:" << endl;
        for (uintT i = 0; i < GA.n; i++) {
            //iteration only through deficits only
            if (deficits.d[i] && blockingSearch[i] < INT_MAX) {
                uintT curNode = i;
                uintT nextNode;
                //go back through edges until excess is found (back propagation))
                //works only in case of unit flows possible
                while (excesses.d[curNode] != true) {
                    cout << "Node: " << curNode << ":" << endl;
                    nextNode = blockingSearch[curNode];
                    //we don't know what edge is correct - so iterate 
                    //through accessible input and inverse input
                    int j = 0;
                    while (j < GA.V[nextNode].outDegree &&
                            (!isOutAdmissible(GA,nextNode,j) || 
                            GA.V[nextNode].getOutNeighbor(j) != curNode)) {
                        j++;
                    }
                    if (j == GA.V[nextNode].outDegree) {
                        j = 0;
                        while (j < GA.V[nextNode].inDegree &&
                                (!isInverseAdmissible(GA,nextNode,j) || 
                                GA.V[nextNode].getInNeighbor(j) != curNode)) {
                            j++;
                        }
                        if (j == GA.V[curNode].inDegree) {
                            cout << "Error: missing edge in raising flows" << endl;
                            exit(1);
                        }
                    }
                    //now we have ID of an edge
                    GA.getInInfo(curNode,j).flow++;
                    curNode = nextNode;
                }
                cout << "Iteration finished. Final excess: " << curNode << endl;
                printGraph(GA);
            }
            cout << endl;
        }
    
        //recalculate excesses and deficits
        calculateExcesses(GA,excesses,deficits);
    }
    
    deficits.del();
    excesses.del();
}

void Compute(graph& GA, commandLine P) {
    
    long source = P.getOptionLongValue("-s",0);
    
    long n = GA.n;
    long m = GA.m;
    
    long target = P.getOptionLongValue("-t",n-1);
    
    //initialize ShortestPathLen to "infinity"
    {parallel_for(long i = 0; i<m; i++) {
            GA.E[i].capacity = 1;
            GA.E[i].flow = 0;
            GA.E[i].lower = 0;
    }}
    
    intE* p = newA(intE, n);
    intE* supply = newA(intE, n);
    {parallel_for(long i = 0; i < n; i++) {
            p[i] = 0;
            supply[i] = 0;
    }}
    supply[source] = 1;
    supply[target] = -1;
    
    GA.p = p;
    GA.supply = supply;
    
    //getMaxWeight
    long maxWeight = 0;
    for(long i = 0; i < GA.n; i++) {
        for (long j = 0; j < GA.V[i].getOutDegree(); j++) {
            maxWeight = (maxWeight < GA.getOutInfo(i,j).weight) ? GA.getOutInfo(i,j).weight : maxWeight;
        }
    }
    
    double epsilon = maxWeight;
    
    while (epsilon > (double)1./n) {
        refine(GA, epsilon);
        epsilon = epsilon / 2.;
    }
    
    cout << "Done." << endl;
    
    free(p);
}

