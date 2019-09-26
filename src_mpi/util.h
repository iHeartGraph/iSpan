#ifndef __UTIL_H__
#define __UTIL_H__
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h>
#include <algorithm>
#include <mpi.h>

#define LOCK(vert, lock) while(!__sync_bool_compare_and_swap(lock+vert,0,-1))
#define UNLOCK(vert, lock) lock[vert]=0

///change to int for SCC
//typedef long index_t;
//typedef long vertex_t;
//typedef double path_t;
//typedef long depth_t;
//
//typedef int index_t;
//typedef int vertex_t;
//typedef double path_t;
//typedef signed char depth_t;
//typedef int color_t;

typedef int index_t;
typedef int vertex_t;
typedef double path_t;
typedef signed char depth_t;
typedef int color_t;
typedef unsigned int long_t;

#define INFTY (float)10000000 
#define NEGATIVE (int)-1
#define ORPHAN	(unsigned char)254
#define UNVIS		(long)-1
#define DEBUG 0
#define VERBOSE 0 
#define OUTPUT_TIME 1
//#define ALPHA 30 
//#define BETA 200 
//#define GAMMA 10
//#define THETA 0.01
#define TRIM_TIMES 3
inline off_t fsize(const char *filename) {
	struct stat st; 
	if (stat(filename, &st) == 0)
		return st.st_size;
	return -1; 
}

#endif
