#ifndef __GRAPH_H__
#define __GRAPH_H__
#include "util.h"
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "wtime.h"
class graph
{
	public:
		index_t *fw_beg_pos;
		vertex_t *fw_csr;
		index_t *bw_beg_pos;
		vertex_t *bw_csr;
		path_t *weight;
		index_t vert_count;
		index_t edge_count;
		vertex_t *src_list;
		index_t src_count;

	public:
		graph(){};
		~graph(){};
		graph(const char *fw_beg_file, 
				const char *fw_csr_file,
                const char *bw_beg_file,
                const char *bw_csr_file);
		void gen_src(){};
		void groupby(){};
};
#endif
