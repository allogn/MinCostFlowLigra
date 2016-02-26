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

inline bool isOutAdmissible(graph& GA, uintT node, uintT outEdge) {
    uintT neighbor = GA.V[node].getOutNeighbor(outEdge);
    edgeInfo ei = GA.getOutInfo(node,outEdge);
    double c_p = ei.weight + GA.p[neighbor] - GA.p[node];
    intE residual = ei.capacity - ei.flow;
    return (c_p < 0 && residual > 0);
}
inline bool isInverseAdmissible(graph& GA, uintT node, uintT inEdge) {
    uintT neighbor = GA.V[node].getInNeighbor(inEdge);
    edgeInfo ei = GA.getInInfo(node,inEdge);
    double c_p = -ei.weight + GA.p[neighbor] - GA.p[node];
    intE residual = ei.flow - ei.lower; //unit case!
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

inline int ftoi_sse1(float f)
{
    return _mm_cvtt_ss2si(_mm_load_ss(&f));     // SSE1 instructions for float->int
}

uintT raise_potentials(graph& GA, const double epsilon, 
    const vertexSubset& excesses, const vertexSubset& deficits) {
    
    if (testing && verbose) log_file << "Raise potentials..." << endl;

    // boolean visited array
    bool visited[GA.n];
    
    //calculate shortest distances
    NodeList buckets(GA.n);
    for (uintT i = 0; i < GA.n; i++) {
        visited[i] = false;
        if (excesses.d[i]) {
            buckets.addToBucket(0, i);
        }
    }

    //calculate shortest path in Residual (!) graph
    uintT length_to_deficit;
    while (true) {
        uintT i = buckets.getNearest();
        
        if (visited[i] == true) {
            buckets.popNearest();
            continue;
        }
        visited[i] = true;
        buckets.saveAndPopNearest();
        
        if (deficits.d[i]) {
            length_to_deficit = buckets.smallest_dist;
            break;
        }

        //iterate through residual edges in both directions!
        for (uintT j = 0; j < GA.V[i].outDegree; j++) {
            uintT neighbour = GA.V[i].getOutNeighbor(j);
            if (! visited[neighbour]) { //@todo change if statements (with edge) and maybe faster
                edgeInfo ei = GA.getOutInfo(i,j);
                if (ei.flow < ei.capacity) {
                    //length function: l = floor(c_p/epsilon) + 1
                    double c_p = ei.weight + GA.p[neighbour] - GA.p[i];
                    uintT length = floor(c_p/epsilon) + 1;
                    uintT newDist = buckets.smallest_dist + length;
                    buckets.addToBucket(newDist, neighbour);
                }
            }
        }
        //inverted edges (see G_f condition)
        for (uintT j = 0; j < GA.V[i].inDegree; j++) {
            uintT neighbour = GA.V[i].getInNeighbor(j);
            if (!visited[neighbour]) {
                edgeInfo ei = GA.getInInfo(i,j);
                if (ei.flow > ei.lower) {
                    //length function: l = floor(c_p/epsilon) + 1
                    double c_p = - ei.weight + GA.p[neighbour] - GA.p[i];
                    uintT length = floor(c_p/epsilon) + 1;
                    uintT newDist = buckets.smallest_dist + length;
                    buckets.addToBucket(newDist, neighbour);
                }
            }
        }
    }
    if (testing && verbose) {
        log_file << "Shortest path done. MaxLength = " << length_to_deficit << endl;
//        for (uintT i = 0; i < GA.n; i++) {
//            log_file << i << ":" << ShortestDistance[i] << endl;
//        }
        log_file << endl;
        
        bool t = test_shortest_path(GA, buckets, epsilon, excesses, deficits);
        log_file << "Testing Shortest Path: " << t << endl;
        if (!t) {
            cout << "Shortest Path test fail" << endl;
            exit(1);
        }
    }

    //delete everything else from buckets without moving to "visited"
    buckets.clearBuckets(visited);
    
    //raise potentials for all visited vertices
    while (buckets.visited != NULL) {
        uintT i = buckets.getIdVisited();
        uintT dist = buckets.getDistVisited();
        buckets.popVisited();
        GA.p[i] += (length_to_deficit - dist)*epsilon;
    }
    
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

