// Alvis Logins (c) 2015

#include "MinCirculation.h"

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
    double c_p = GA.getOutInfo(node,outEdge).weight + GA.p[neighbor] - GA.p[node];
    intE residual = GA.getOutInfo(node,outEdge).capacity - GA.getOutInfo(node,outEdge).flow;
    if (testing && verbose)
        log_file << node << "->" << neighbor << " admissible: " << (c_p < 0 && residual > 0)
            << "(c_p = " << c_p << ", residual: " << residual << ")\n";
    return (c_p < 0 && residual > 0);
}
bool isInverseAdmissible(graph& GA, uintT node, uintT inEdge) {
    uintT neighbor = GA.V[node].getInNeighbor(inEdge);
    double c_p = -GA.getInInfo(node,inEdge).weight + GA.p[neighbor] - GA.p[node];
    intE residual = GA.getInInfo(node,inEdge).flow - GA.getInInfo(node,inEdge).lower; //unit case!
    if (testing && verbose) log_file << node << "<-" << neighbor << " admissible: " << (c_p < 0 && residual > 0)
        << "(c_p = " << c_p << ", residual: " << residual << ")\n";
    return (c_p < 0 && residual > 0);
}

void dfs(graph& GA, uintT* blockingFlow, uintT nodeId,
        const vertexSubset& excesses, const vertexSubset& deficits) {
    if (testing && verbose) log_file << "  DFS " << nodeId << ":\n";
    
    queue<uintT> node_queue;
    node_queue.push(nodeId);
    
    while(!node_queue.empty()) {
        uintT node = node_queue.front();
        node_queue.pop();
        
        for (uintT i = 0; i < GA.V[node].outDegree; i++) {
            //for each out edge check if it is admissible
            //if so - update distance and run dfs
            uintT neighbor = GA.V[node].getOutNeighbor(i);
            if (isOutAdmissible(GA,node,i) && blockingFlow[neighbor] == UINT_T_MAX/2
                    && !excesses.d[neighbor]) {
                blockingFlow[neighbor] = node;
                node_queue.push(neighbor);
                
                if (deficits.d[neighbor]) {
                    return;
                }
            }
        }
        //the same for inverse edges
        for (uintT i = 0; i < GA.V[node].inDegree; i++) {
            //for each out edge check if it is admissible
            //if so - update distance and run dfs
            uintT neighbor = GA.V[node].getInNeighbor(i);
            if (isInverseAdmissible(GA,node,i) && blockingFlow[neighbor] == UINT_T_MAX/2
                    && !excesses.d[neighbor]) {
                blockingFlow[neighbor] = node;
                node_queue.push(neighbor);
                
                if (deficits.d[neighbor]) {
                    return;
                }
            }
        }
    }
}

void calculateExcesses(graph& GA, vertexSubset& excesses, vertexSubset& deficits) {
    if (testing && verbose) log_file << "Calculate excesses..." << endl;
    
    uintT totalExcesses = 0;
    uintT totalDeficits = 0;
    for (uintT i = 0; i < GA.n; i++)
    {
        excesses.d[i] = false;
        deficits.d[i] = false;
        intT flow_sum = 0;
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
    
    if (testing && verbose) log_file << "excesses: " ;
    for (uintT i = 0; i < GA.n; i++) {
        if (testing) log_file << i << ":" << excesses.d[i] << " ";
    }
    if (testing && verbose) log_file << endl;
    if (testing && verbose) log_file << "deficits: " ;
    for (uintT i = 0; i < GA.n; i++) {
        if (testing && verbose) log_file << i << ":" << deficits.d[i] << " ";
    }
    if (testing && verbose) log_file << endl;
}

void printGraph(graph GA) {
    for (uintE i = 0; i < GA.n; i++) {
        log_file << "#" << i << " (" << GA.p[i] << ", supply: " << GA.supply[i] << "):\t\n";
        log_file << "    " << "OUT:\n";
        for (uintE j = 0; j < GA.V[i].outDegree; j++) {
            log_file << "        " << GA.V[i].getOutNeighbor(j) 
                    << " (w: " << GA.getOutInfo(i,j).weight << "; "
                    << " f: " << GA.getOutInfo(i,j).flow << "; "
                    << " c: " << GA.getOutInfo(i,j).capacity << "; "
                    << " l: " << GA.getOutInfo(i,j).lower << ")\n";
        }
        log_file << "    " << "IN:\n";
        for (uintE j = 0; j < GA.V[i].inDegree; j++) {
            log_file << "        " << GA.V[i].getInNeighbor(j) 
                    << " (w: " << GA.getInInfo(i,j).weight << "; "
                    << " f: " << GA.getInInfo(i,j).flow << "; "
                    << " c: " << GA.getInInfo(i,j).capacity << "; "
                    << " l: " << GA.getInInfo(i,j).lower << ")\n";
        }
    }
}

uintT raise_potentials(graph& GA, const double epsilon, 
    const vertexSubset& excesses, const vertexSubset& deficits) {
    
    if (testing && verbose) log_file << "Raise potentials..." << endl;

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
                        if (newDist >= length_to_deficit) { 
                            newDist = length_to_deficit; 
                        } else if (ShortestDistance[neighbour] > newDist) {
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
                        if (newDist >= length_to_deficit) { 
                            newDist = length_to_deficit; 
                        } else if (ShortestDistance[neighbour] > newDist) {
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
    if (testing && verbose) {
        log_file << "Shortest path done. MaxLength = " << length_to_deficit << endl;
//        for (uintT i = 0; i < GA.n; i++) {
//            log_file << i << ":" << ShortestDistance[i] << endl;
//        }
        log_file << endl;
    }

    //raise potentials
    for (uintT i = 0; i < GA.n; i++) {
        if (ShortestDistance[i] < length_to_deficit) {
            GA.p[i] += (length_to_deficit - ShortestDistance[i])*epsilon;
        }
    }

    free(ShortestDistance);
}

void raise_flows(graph& GA, const double epsilon, 
    const vertexSubset& excesses, const vertexSubset& deficits) {
    //calculate blocking flow via depth-first search in admissible graph
    // depth-first search: creating an array of pointers to parents
    // can be extended to the case of non-unit capacity by saving 
    // a pointer to a parent together with total flow
    //terminate when a deficit reached
    uintT* blockingSearch = newA(uintT, GA.n);
    for (uintT i = 0; i < GA.n; i++) blockingSearch[i] = UINT_T_MAX/2;
    for (uintT i = 0; i < GA.n; i++) {
        if (excesses.d[i]) {
            dfs(GA, blockingSearch, i, excesses, deficits);
        }
    }
    if (testing && verbose) {
        log_file << "Blocking flow array: ";
        for (uintT i = 0; i < GA.n; i++) {
            log_file << i << ":" << blockingSearch[i] << " ";
        }
        log_file << endl;
        log_file << "blocking flow done." << endl;
        log_file << "Paths: \n";
    }
    
    //raise flow: back propagation from each deficit
    for (uintT i = 0; i < GA.n; i++) {
        if (deficits.d[i] && blockingSearch[i] < UINT_T_MAX/2) {
            uintT curNode = i;
            while (excesses.d[curNode] != true) {
                if (testing) log_file << curNode << "<-";
                curNode = blockingSearch[curNode];
            }
            if (testing) log_file << curNode << endl;
        }
    }

    //now raise flows
    if (testing && verbose) log_file << "Raising flows... Path:" << endl;
    for (uintT i = 0; i < GA.n; i++) {
        //iteration only through deficits only
        if (deficits.d[i] && blockingSearch[i] < UINT_T_MAX/2) {
            uintT curNode = i;
            uintT nextNode;
            //go back through edges until excess is found (back propagation))
            //works only in case of unit flows possible
            while (excesses.d[curNode] != true) {
                if (testing && verbose) log_file << "Node: " << curNode << ":" << endl;
                nextNode = blockingSearch[curNode];
                //we don't know what edge is correct - so iterate 
                //through accessible input and inverse input
                uintT j = 0;
                while (j < GA.V[nextNode].outDegree &&
                        (!isOutAdmissible(GA, nextNode, j) ||
                        GA.V[nextNode].getOutNeighbor(j) != curNode)) {
                    j++;
                }
                if (j == GA.V[nextNode].outDegree) {
                    j = 0;
                    while (j < GA.V[nextNode].inDegree &&
                            (!isInverseAdmissible(GA, nextNode, j) ||
                            GA.V[nextNode].getInNeighbor(j) != curNode)) {
                        j++;
                    }
                    if (j == GA.V[nextNode].inDegree) {
                        if (testing && verbose) log_file << "Error: missing edge in raising flows" << endl;
                        exit(1);
                    }
                    GA.getInInfo(nextNode, j).flow--;
                } else {
                    GA.getOutInfo(nextNode, j).flow++;
                }

                curNode = nextNode;
            }
            if (testing && verbose) {
                log_file << "Iteration finished. Final excess: " << curNode << endl;
                printGraph(GA);
            }
        }
    }
        
    free(blockingSearch);
}

void refine(graph& GA, const double epsilon) {
    
    if (testing) {
        if (verbose) log_file << "New iteration: e = " << epsilon << endl;
        printGraph(GA);
        
        //check feasibility of c_p
        if (!is_feasible(-2*epsilon,GA)) {
            cout << "Error: not a feasible flow" << endl;
            exit(1);
        }
    }
    
    //loop through arcs, saturate those with c_p less than 0 and add excesses
    for (uintT i = 0; i < GA.n; i++)
    {
        for (uintT j = 0; j < GA.V[i].getOutDegree(); j++) {
            uintT neighbour = GA.V[i].getOutNeighbor(j);
            double c_p = GA.getOutInfo(i,j).weight + GA.p[neighbour] - GA.p[i];
            if (isOutAdmissible(GA,i,j)) {
                GA.getOutInfo(i,j).flow = GA.getOutInfo(i,j).capacity;
            }
            if (testing && verbose) log_file << "c_p("<< i << "," << neighbour << ") = " << c_p << endl;
        }
    }
    //iterate through inverted edges
    if (testing && verbose) log_file << "Inverted:" << endl;
    for (uintT i = 0; i < GA.n; i++)
    {
        for (uintT j = 0; j < GA.V[i].inDegree; j++) {
            uintT neighbour = GA.V[i].getInNeighbor(j);
            double c_p = - GA.getInInfo(i,j).weight + GA.p[neighbour] - GA.p[i];
            if (isInverseAdmissible(GA,i,j)) {
                GA.getInInfo(i,j).flow = GA.getInInfo(i,j).lower;
            }
            if (testing && verbose) log_file << "c_p("<< neighbour << "," << i << ") = " << c_p << endl;
        }
    }
    
    if (testing && !is_feasible(-epsilon,GA)) {
        cout << "Error: not a feasible flow after saturating" << endl;
        exit(1);
    }
    
    vertexSubset excesses(GA.n);
    vertexSubset deficits(GA.n);
    excesses.toDense();
    deficits.toDense();
    calculateExcesses(GA,excesses,deficits);
    
    while (excesses.isEmpty() == false) {
        
        uintE target_id = raise_potentials(GA,epsilon,excesses,deficits);
        
        if (testing && !is_feasible(-epsilon,GA)) {
            cout << "Error: not a feasible flow after raising potentials" << endl;
            exit(1);
        }
        
        raise_flows(GA,epsilon,excesses, deficits);
    
        //recalculate excesses and deficits
        calculateExcesses(GA,excesses,deficits);
    }
    
    deficits.del();
    excesses.del();
    
    if (testing) {
        if (!is_flow(GA)) {
            cout << "Error: not a flow" << endl;
            exit(1);
        };
    }
}

void Compute(graph& GA, commandLine P) {
    
    uintT source = P.getOptionLongValue("-s",0);
    char* iFile = P.getArgument(0);
    
    testing = P.getOptionIntValue("-test",0);
    verbose = P.getOptionIntValue("-verbose",0);
    string log_filename = P.getOptionValue("-log","log.txt");
    string profiling_file = P.getOptionValue("-profile","profile.txt");
    if (testing) {
        log_file.open(log_filename.c_str());
    }
    
    Timers timers;
    
    uintT n = GA.n;
    uintT m = GA.m;
    
    uintT target = P.getOptionLongValue("-t",n-1);
    
    //initialize ShortestPathLen to "infinity"
    {parallel_for(long i = 0; i<m; i++) {
            GA.E[i].capacity = 1;
            GA.E[i].flow = 0;
            GA.E[i].lower = 0;
    }}
    
    double* p = newA(double, n);
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
    uintT maxWeight = 0;
    for(uintT i = 0; i < GA.n; i++) {
        for (uintT j = 0; j < GA.V[i].getOutDegree(); j++) {
            maxWeight = (maxWeight < GA.getOutInfo(i,j).weight) ? GA.getOutInfo(i,j).weight : maxWeight;
        }
    }
    
    double epsilon = maxWeight;
    
    timers.total->start();
    while (epsilon > (double)1./n) {
        refine(GA, epsilon);
        epsilon = epsilon / 2.;
    }
    timers.total->stop();
    
    if (testing) printGraph(GA);
    bool checked = true;
    if (testing) {
        //flow check
        string sol = P.getOptionValue("-sol", "");
        if (sol == "") {
            cout << "Error: Solution source must be provided" << endl;
            exit(1);
        }
        log_file << "Running unit tests" << endl;
        if (!flow_check(GA,sol)) {
            log_file << "Flow check failed" << endl;
            checked = false;
        } else {
            log_file << "Flow check successful" << endl;
        }
        
        if (checked) cout << "Testing done." << endl;
    }
    
    cout << "Done." << endl;
    cout << "Total time: " << timers.total->realTime() << endl;
    
    intT self_totalcost = 0;
    for (uintT node = 0; node < GA.n; node++) {
        for (uintT edge = 0; edge < GA.V[node].outDegree; edge++) {
            edgeInfo ei = GA.getOutInfo(node,edge);
            self_totalcost += ei.weight*ei.flow;
        }
    }
    cout << "Total cost: " << self_totalcost << endl;
    
    timers.print(profiling_file, GA);
    log_file.close();
}

