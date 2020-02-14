#ifndef TRIANGLE_DIST_RUN_HPP
#define TRIANGLE_DIST_RUN_HPP

#include "base_struct.hpp"
#include "source.hpp"
#include "worker.hpp"
#include <vector>
#include <iostream>
#include <unordered_map>
#include <sys/time.h>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <limits>

// run actual mpi implementation
void run_mpi(const char* filename, MPIIO &hIO, int workerNum, int memSize, int lenBuf, unsigned int seed, long double * moments, double &srcCompCost, double &workerCompCostMax, double &workerCompCostSum);

// run experiment with mpi
void run_exp (const char* input, const char* outPath, MPIIO &hIO, int workerNum, int memSize, int repeat, int bufLen=1000);

void print_cnt(const char* outPath, const long double * moments, int id = 0);

#endif //TRIANGLE_DIST_RUN_HPP
