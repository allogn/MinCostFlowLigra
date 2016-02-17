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
        total = new TimeReport("Total time:",f,false);
    }
    
    void print(string filename) {
        f.open(filename.c_str());
        f << "Total time: " << total->realTime() << endl;
        f.close();
    }
    
    ~AllTimes() {
        delete total;
    }
};

typedef struct AllTimes Timers;

#endif	/* TIMING_H */

