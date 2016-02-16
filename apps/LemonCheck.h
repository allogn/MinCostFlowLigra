// Alvis Logins (c) 2015

#include "MinCirculation.h"

bool flow_check(graph& GA, string& filename) {
    
    //build lemon graph
    ListDigraph g;
    ListDigraph::ArcMap<int> weight(g);
    ListDigraph::ArcMap<int> flow(g);
    for (uintT node = 0; node < GA.n; node++) {
        g.addNode();
    }
    for (uintT node = 0; node < GA.n; node++) {
        for (uintT edge = 0; edge < GA.V[node].outDegree; edge++) {
            ListDigraph::Arc e = g.addArc(g.nodeFromId(node),
                    g.nodeFromId(GA.V[node].getOutNeighbor(edge)));
            weight[e] = GA.getOutInfo(node,edge).weight;
            flow[e] = GA.getOutInfo(node,edge).flow;
        }
    }
    
    //load lemon graph with results
    log_file << "Loading answers from " << filename << endl;
    int lemon_result;
    long totalCost;
    ListDigraph result;
    ListDigraph::ArcMap<int> rflow(result);
    digraphReader(result, filename)
            .arcMap("flow", rflow)
            .attribute("totalCost",totalCost)
            .attribute("result",lemon_result)
            .run();
    log_file << "Answers loaded" << endl;
    
    if (lemon_result != 1) {
        log_file << "Error: wrong lemon answer. Result = " << lemon_result << endl;
        log_file << "Probably wrong input file for lemon" << endl;
        return false;
    }
    
    //calculate total cost in internal graph
    intT self_totalcost = 0;
    for (uintT node = 0; node < GA.n; node++) {
        for (uintT edge = 0; edge < GA.V[node].outDegree; edge++) {
            edgeInfo ei = GA.getOutInfo(node,edge);
            self_totalcost += ei.weight*ei.flow;
        }
    }
    
    if (self_totalcost != totalCost) {
        log_file << "Different total costs: internal=" << self_totalcost 
                << " lemon=" << totalCost << endl;
        return false;
    }
    return true;
}