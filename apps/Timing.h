/* 
 * File:   timing.h
 * Author: alvis
 *
 * Created on February 17, 2016, 6:13 PM
 */

#ifndef TIMING_H
#define	TIMING_H

#include <lemon/time_measure.h>
#include <fstream>

struct AllTimes {
    
    ofstream f;
    Timer* total;
    
    AllTimes() {
        total = new Timer(false);
    }
    
    void print(string filename, graph& GA) {
        f.open(filename.c_str());
        
        f << "Nodes: " << GA.n << endl;
        f << "Edges: " << GA.m << endl;
        
        f << "Total time: " << total->realTime() << endl;
        f.close();
    }
    
    ~AllTimes() {
        delete total;
    }
};

typedef struct AllTimes Timers;

#endif	/* TIMING_H */

