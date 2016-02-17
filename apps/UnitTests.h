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

#endif	/* UNITTESTS_H */

