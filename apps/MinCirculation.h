/* 
 * File:   MinCirculation.h
 * Author: alvis
 *
 * Created on February 15, 2016, 4:00 PM
 * 
 * Usage example:
 * -rounds 0 -test 1 -verbose 1 ../../data/tests/small_4.adj
 */

#ifndef MINCIRCULATION_H
#define	MINCIRCULATION_H

#define WEIGHTED 1
#define DEBUG
#include <vector>
#include <fstream>

using namespace std;

#include "ligra.h"
#include <queue>

#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/cost_scaling.h>
#include <emmintrin.h>
using namespace std;
using namespace lemon;

bool testing;
ofstream log_file;
bool verbose;

#include "Timing.h"
#include "NodeList.h"
#include "LemonCheck.h"
#include "UnitTests.h"

#endif	/* MINCIRCULATION_H */

