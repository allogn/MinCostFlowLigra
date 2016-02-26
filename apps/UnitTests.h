/* 
 * File:   UnitTests.h
 * Author: alvis
 *
 * Created on February 16, 2016, 6:04 PM
 */

#ifndef UNITTESTS_H
#define	UNITTESTS_H

bool is_feasible(double threshold, graph& GA) {
    for (uintT i = 0; i < GA.n; i++) {
        for (uintT j = 0; j < GA.V[i].outDegree; j++) {
            if (GA.getOutInfo(i,j).flow < GA.getOutInfo(i,j).capacity) {
                intT c_p = GA.getOutInfo(i,j).weight - GA.p[i] + GA.p[GA.V[i].getOutNeighbor(j)];
                if (c_p < threshold) {
                    if (verbose) log_file << "c_p for forward edge " << i << " -> " << GA.V[i].getOutNeighbor(j)
                            << " is negative" << endl;
                    return false;
                }
            }
        }
    }
    for (uintT i = 0; i < GA.n; i++) {
        for (uintT j = 0; j < GA.V[i].inDegree; j++) {
            if (GA.getInInfo(i,j).flow > 0) {
                intT c_p = - GA.getInInfo(i,j).weight - GA.p[i] + GA.p[GA.V[i].getInNeighbor(j)];
                if (c_p < threshold) {
                    if (verbose) log_file << "c_p for backward edge " << i << " -> " << GA.V[i].getInNeighbor(j)
                            << " is negative" << endl;
                    return false;
                }
            }
        }
    }
    return true;
}

bool is_flow(graph& GA) {
    intT flow = 0;
    for (uintT i = 0; i < GA.n; i++) {
        flow = 0;
        for (uintT j = 0; j < GA.V[i].outDegree; j++) {
            flow -= GA.getOutInfo(i,j).flow;
        }
        for (uintT j = 0; j < GA.V[i].inDegree; j++) {
            flow += GA.getInInfo(i,j).flow;
        }
        flow += GA.supply[i];
        if (flow != 0) {
            if (verbose) log_file << "error: not a flow in node #" << i << " flow: " << flow << endl;
            return false;
        }
    }
    return true;
}

bool test_shortest_path(graph& GA, NodeList& buckets, const double epsilon, 
    const vertexSubset& excesses, const vertexSubset& deficits) {
    //calculate shortest distances
    uintT* ShortestDistance = newA(unsigned long,GA.n);
    vertexSubset Frontier(GA.n);
    Frontier.toDense();
    for(uintT i = 0; i < GA.n; i++) Frontier.d[i] = excesses.d[i];
    Frontier.m = excesses.m;

    for (uintT i = 0; i < GA.n; i++) {
        ShortestDistance[i] = UINT_T_MAX;
    }

    //calculate shortest path in Residual (!) graph
    uintT length_to_deficit = UINT_T_MAX;
    while (!Frontier.isEmpty()) {
        vertexSubset nextIterSubset(GA.n);
        nextIterSubset.m = 0;
        nextIterSubset.toDense();
        bool gotDeficit = false;
        for (uintT i = 0; i < GA.n; i++) {
            //if vertex in the active subset
            if (Frontier.d[i]) {
                //if source - zero distance by definition
                if (excesses.d[i]) {
                    ShortestDistance[i] = 0;
                }

                if (deficits.d[i]) {
                    length_to_deficit = (length_to_deficit > ShortestDistance[i])?
                        ShortestDistance[i] : length_to_deficit;
                }

                //iterate through residual edges to both directions!
                for (uintT j = 0; j < GA.V[i].outDegree; j++) {
                    if (GA.getOutInfo(i,j).flow < GA.getOutInfo(i,j).capacity) {
                        uintT neighbour = GA.V[i].getOutNeighbor(j);

                        //length function: l = floor(c_p/epsilon) + 1
                        double c_p = GA.getOutInfo(i,j).weight + GA.p[neighbour] - GA.p[i];
                        uintT length = floor(c_p/epsilon) + 1;
                        
                        uintT newDist = ShortestDistance[i] + length;
                        if (newDist > 3*GA.n) { newDist = 3*GA.n; }
                        if (ShortestDistance[neighbour] > newDist) {
                            ShortestDistance[neighbour] = newDist;
                            if (!nextIterSubset.d[neighbour]) {
                                nextIterSubset.d[neighbour] = true;
                                nextIterSubset.m++;
                            }
                        }
                    }
                }
                //inverted edges (see G_f condition)
                for (uintT j = 0; j < GA.V[i].inDegree; j++) {
                    if (GA.getInInfo(i,j).flow > GA.getInInfo(i,j).lower) {
                        uintT neighbour = GA.V[i].getInNeighbor(j);

                        //length function: l = floor(c_p/epsilon) + 1
                        double c_p = - GA.getInInfo(i,j).weight + GA.p[neighbour] - GA.p[i];
                        uintT length = floor(c_p/epsilon) + 1;

                        uintT newDist = ShortestDistance[i] + length;
                        if (newDist > 3*GA.n) { newDist = 3*GA.n; }
                        if (ShortestDistance[neighbour] > newDist) {
                            ShortestDistance[neighbour] = newDist;
                            if (!nextIterSubset.d[neighbour]) {
                                nextIterSubset.d[neighbour] = true;
                                nextIterSubset.m++;
                            }
                        }
                    }
                }
            }
        }
        if (gotDeficit) break;
        Frontier.del();
        Frontier = nextIterSubset;
    }
    
    Node* node = buckets.visited;
    bool check[GA.n];
    for (uintT i = 0; i < GA.n; i++) {
        check[i] = false;
    }
    while (node != NULL) {
        check[node->nodeId] = true;
        if (ShortestDistance[node->nodeId] != node->dist) {
            return false;
        }
        node = node->next;
    }
    for(uintT i = 0; i < GA.n; i++) {
        if (ShortestDistance[i] < length_to_deficit && check[i] == false) {
            return false;
        }
    }
    return true;
}

#endif	/* UNITTESTS_H */

