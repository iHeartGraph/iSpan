#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <omp.h>
#include <algorithm>
#include "wtime.h"
#include "graph.h"
//#include "frontier_queue.h"
#include "scc_common.h"
//#include "trim_1.h"
#include "trim_1_gfq.h"
#include "trim_2_3.h"
#include "color_propagation.h"
#include "fw_bw.h"
#include "openmp_wcc.hpp"
#define INF -1

//0 trim, 1 largest SCC, 2 small SCC, 3 total time
//4 trim_size_1, 5 trim_size_2, 6 pivot_selection, 7 fw_bfs, 8 bw_bfs, 9 color propagation, 10 color identify, 11 color_init 

void scc_detection(
        const graph *g,
        const int alpha, 
        const int beta,
        const int gamma,
        const double theta,
        const index_t thread_count,
        double *avg_time
        )
{
    const index_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;
    if(DEBUG)
        printf("vert_count = %d, edge_count = %d, avg_degree = %.3lf\n", vert_count, edge_count, avg_degree);
    //step 0: initialization
//    return;

    index_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    index_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    if(VERBOSE)
    {
        for(int i=fw_beg_pos[vert_count]; i<fw_beg_pos[vert_count+1]; ++i)
            printf("%d\n", fw_csr[i]);
    }

    index_t *scc_id = new index_t[vert_count + 1];
   
    index_t *color = new index_t[vert_count + 1];
//    index_t *color_times = new index_t[vert_count + 1];

    index_t *max_pivot_list = new index_t[thread_count];
    index_t *max_degree_list = new index_t[thread_count];
    

//    index_t *mul_degree = new index_t[vert_count + 1];
//    index_t *degree_prop = new index_t[vert_count + 1];
	
    depth_t *fw_sa;
	depth_t *bw_sa;

//    fw_sa = new depth_t[vert_count+1];
//    bw_sa = new depth_t[vert_count+1];

	if(posix_memalign((void **)&fw_sa,getpagesize(),
		sizeof(depth_t)*(vert_count+1)))
		perror("posix_memalign");
	
    if(posix_memalign((void **)&bw_sa,getpagesize(),
		sizeof(depth_t)*(vert_count+1)))
		perror("posix_memalign");
    
    index_t *small_queue = new index_t[vert_count + 1];
    index_t *temp_queue = new index_t[vert_count + 1];
    index_t *inter_queue = new index_t[vert_count + 1];
    index_t *wcc_fq= new index_t[vert_count + 1];
    
    index_t *thread_bin = new index_t[thread_count];
    index_t *prefix_sum = new index_t[thread_count];
	
    index_t *front_comm=new index_t[thread_count];	
	index_t *work_comm=new index_t[thread_count];
    bool *color_change = new bool[thread_count];
    memset(color_change, 0, sizeof(bool) * thread_count);

    //WCC + FW-BW
//    index_t *wcc_color = new index_t[vert_count + 1];
//    index_t *color_redirect = new index_t[vert_count + 1];
//    bool *is_redirect = new bool[thread_count];
//    color_t *global_color = new color_t[1];
//    global_color[0] = 0;

    vertex_t wcc_fq_size = 0;

//    memset(fw_sa, -1, sizeof(depth_t)*vert_count);
//	memset(bw_sa, -1, sizeof(depth_t)*vert_count);
//    memset(scc_id, 0, sizeof(index_t) * (vert_count + 1));

//Initialization
    
    if(DEBUG)
    {
        printf("Initialization\n");
    }

//    #pragma omp parallel for
    for(index_t i=0; i<vert_count + 1; ++i)
    {
        color[i] = i;
//            color_times[i] = 0;
        fw_sa[i] = -1;
        bw_sa[i] = -1;
        scc_id[i] = 0;
//            wcc_color[i] = -1;
//            color_redirect[i] = i;
    }

    //step 1: trim size-1
    //step 2: largest_scc, asynchronize, direction optimized
    //step 3: trim size-1 && size-2 && size-3
    /// generate new FQ
    //step 4: small sccs, graph coloring
    
    if(DEBUG)
    {
        printf("Parallel starts\n");
    }
    vertex_t vertex_fw = 0;
    vertex_t vertex_bw = 0;
    index_t size_3_1 = 0;
    index_t size_3_2 = 0;
    bool changed = false;
    double end_time;
    double start_time = wtime();
    #pragma omp parallel \
    num_threads(thread_count) 
//    \shared(global_color)
    {
        const index_t tid = omp_get_thread_num();
        index_t step = vert_count / thread_count;
        index_t vert_beg = tid * step;
        index_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
        double time_size_1_first;
        double time_size_1;
        double time_fw;
        double time_bw; 
        double time_size_2;
        double time_size_3;
        double time_gfq;
        double time_color_1;
        double time_color_2;
        double time_color;
        double pivot_time;
        double time_color_init;
        double time_wcc;
        double time_mice_fw_bw; 
        const vertex_t upper_bound = vert_count / thread_count * 5;
        vertex_t *thread_queue = new vertex_t[upper_bound];
        
        
        double time = wtime();
// not using frontier queue for trim
//
        index_t trim_times = 1;

// change to control by trimmed vertices
        
        trim_1_first(scc_id,
                fw_beg_pos,
                bw_beg_pos,
                vert_beg,
                vert_end);
        trim_times ++;
        #pragma omp barrier
        if(tid == 0)
        {
            time_size_1_first = wtime() - time;
        }
//            trim_1_normal_gfq(scc_id,
//                    fw_beg_pos,
//                    bw_beg_pos,
//                    vert_beg,
//                    vert_end,
//                    fw_csr,
//                    bw_csr,
//                    thread_count,
//                    small_queue,
//                    thread_bin,
//                    prefix_sum,
//                    tid);
    
        
//        time = wtime();
//        while(trim_times < 10)
//        {
//            trim_1_normal_only_size(scc_id,
//                    fw_beg_pos,
//                    bw_beg_pos,
//                    vert_beg,
//                    vert_end,
//                    fw_csr,
//                    bw_csr,
//                    thread_count,
//                    thread_bin,
//                    prefix_sum,
//                    tid);
//            #pragma omp barrier
//            vertex_t trimmed_v = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
////FISC plan            
//            if(trimmed_v < theta * vert_count)
//                break;
//
////            //Hong
////            if(trimmed_v < 100)
////                break;
//            trim_times++;
//            #pragma omp barrier
//        }
//        #pragma omp barrier
//        if(tid == 0)
//        {
//            time_size_1 = wtime() - time;
//            if(DEBUG)
//            {
//                printf("trim_1_times, %d\n", trim_times);
//            }
//        }
//
//        trim_1_first_gfq(scc_id,
//                fw_beg_pos,
//                bw_beg_pos,
//                vert_beg,
//                vert_end,
//                thread_count,
//                small_queue,
//                thread_bin,
//                prefix_sum,
//                tid);
/*        
        vertex_t fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
        step = fq_size / thread_count;
        vert_beg = tid * step;
        vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
        
        if(DEBUG)
        {
            if(tid == 0)
            {
                printf("trim_1, 1st, %.3lf, fq_size, %d\n", time_size_1_first * 1000, fq_size);
            }
        }
        vertex_t prev_fq_size = fq_size;
        
        #pragma omp barrier

        //step 1.2: trim size_1
        while(prev_fq_size > 0 && trim_times < TRIM_TIMES)
        {
            time = wtime();
//            trim_1_normal_gfq(scc_id,
//                    fw_beg_pos,
//                    bw_beg_pos,
//                    vert_beg,
//                    vert_end,
//                    fw_csr,
//                    bw_csr,
//                    thread_count,
//                    small_queue,
//                    thread_bin,
//                    prefix_sum,
//                    tid);
            
            trim_1_from_fq_gfq(scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    thread_count,
                    small_queue,
                    thread_bin,
                    prefix_sum,
                    tid,
                    temp_queue);
            
            #pragma omp barrier
            if(tid == 0)
            {
                time_size_1 += wtime() - time;
            }
            fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            step = fq_size / thread_count;
            vert_beg = tid * step;
            vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
            
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("trim_1, 2nd, %.3lf, trimmed vertices, %d\n", time_size_1 * 1000, prev_fq_size - fq_size);
                }
            }
            if(prev_fq_size - fq_size < THETA * prev_fq_size)
                break;
            prev_fq_size = fq_size;
            trim_times++;
            #pragma omp barrier
        }

//        if(tid == 0)
//        {
//            std::sort(small_queue, small_queue+fq_size);
//        }
//        #pragma omp barrier

        if(DEBUG)
        {
            if(tid == 0)
                printf("trim-1 times befor FW-BW, %d\nfq_size, %d\n", trim_times - 1, fq_size);
        }
 
// not using FW-BW to detect the Elephant SCC
//        printf("before fw\n");
        time = wtime();
        //step 2: largest_scc, asynchronize, direction optimized
        //step 2.1: pivot selection
//        vertex_t root = pivot_selection(scc_id,
//                        fw_beg_pos,
//                        bw_beg_pos,
//                        vert_beg,
//                        vert_end,
//                        fw_csr,
//                        bw_csr,
//                        max_pivot_list,
//                        max_degree_list,
//                        tid,
//                        thread_count);

        */

        time = wtime();

        generate_frontier_queue(vert_count,
                scc_id,
                thread_count,
                small_queue,
                thread_bin,
                prefix_sum,
                vert_beg,
                vert_end,
                tid);
        #pragma omp barrier

        if(tid == 0)
        {
            time_size_1 = wtime() - time;
        }
        index_t fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
        step = fq_size / thread_count;
        vert_beg = tid * step;
        vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
        
        if(DEBUG)
        {
            if(tid == 0)
            {
                printf("fq_size, %d\n", fq_size);
            }
        }
        
        time = wtime();
        vertex_t root = pivot_selection_from_fq(scc_id,
                        fw_beg_pos,
                        bw_beg_pos,
                        vert_beg,
                        vert_end,
                        fw_csr,
                        bw_csr,
                        max_pivot_list,
                        max_degree_list,
                        tid,
                        thread_count,
                        small_queue);
//        ,
//                        mul_degree,
//                        degree_prop);
        #pragma omp barrier
        pivot_time = wtime() - time;

        
        //step 2.2: fw bfs
//        time = wtime();
//        fw_bfs(scc_id,
//                fw_beg_pos,
//                bw_beg_pos,
//                vert_beg,
//                vert_end,
//                fw_csr,
//                bw_csr,
//                fw_sa,
//                front_comm,
//                work_comm,
//                root,
//                tid,
//                thread_count,
//                alpha,
//                beta,
//                edge_count,
//                vert_count);
//        #pragma omp barrier
//        time_fw = wtime() - time;
//
//        ///step 2.3: bw bfs 
//        time = wtime();
//        bw_bfs(scc_id,
//                fw_beg_pos,
//                bw_beg_pos,
//                vert_beg,
//                vert_end,
//                fw_csr,
//                bw_csr,
//                fw_sa,
//                bw_sa,
//                front_comm,
//                work_comm,
//                root,
//                tid,
//                thread_count,
//                alpha,
//                beta,
//                edge_count,
//                vert_count);
//        #pragma omp barrier
//        time_bw = wtime() - time;

        time = wtime();
        fw_bfs_fq_queue(scc_id,
                fw_beg_pos,
                bw_beg_pos,
                vert_beg,
                vert_end,
                fw_csr,
                bw_csr,
                fw_sa,
                front_comm,
                work_comm,
                root,
                tid,
                thread_count,
                alpha,
                beta,
                gamma,
                small_queue,
                fq_size,
                avg_degree,
                vertex_fw,
                temp_queue,
                prefix_sum,
                upper_bound,
                thread_queue);
        #pragma omp barrier
        time_fw = wtime() - time;

        
        time = wtime();
        bw_gfq_from_fw(fw_sa, 
                thread_count,
                small_queue,
                thread_bin,
                prefix_sum,
                vert_beg,
                vert_end,
                tid,
                temp_queue);
        #pragma omp barrier
        time_bw = wtime() - time;

        index_t temp_fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
        index_t temp_step = temp_fq_size / thread_count;
        index_t temp_vert_beg = tid * temp_step;
        index_t temp_vert_end = (tid == thread_count - 1 ? temp_fq_size : temp_vert_beg + temp_step);
        
        if(DEBUG)
        {
            if(tid == 0)
            {
                printf("gfq bw, %.3lf, fq_size, %d\n", time_bw * 1000, temp_fq_size);
            }
        }
        #pragma omp barrier
        
//        printf("before fw\n");
        time = wtime();

        ///step 2.3: bw bfs 
        time = wtime();
        bw_bfs_fq_queue(scc_id,
                fw_beg_pos,
                bw_beg_pos,
                temp_vert_beg,
                temp_vert_end,
                fw_csr,
                bw_csr,
                fw_sa,
                bw_sa,
                front_comm,
                work_comm,
                root,
                tid,
                thread_count,
                alpha,
                beta,
                gamma,
                temp_queue,
                temp_fq_size,
                avg_degree,
                vertex_bw,
                inter_queue,
                prefix_sum,
                upper_bound,
                thread_queue);
//        bw_bfs_fq(scc_id,
//                fw_beg_pos,
//                bw_beg_pos,
//                temp_vert_beg,
//                temp_vert_end,
//                fw_csr,
//                bw_csr,
//                fw_sa,
//                bw_sa,
//                front_comm,
//                work_comm,
//                root,
//                tid,
//                thread_count,
//                alpha,
//                beta,
//                temp_queue,
//                temp_fq_size,
//                avg_degree,
//                vertex_bw);
        trim_times = 0;
        #pragma omp barrier
        time_bw += wtime() - time;
        delete[] thread_queue;
    
        vertex_t prev_fq_size = fq_size;

        //step 3: trim size-1 && size-2
        //step 3.1: trim size_1 && gfq
        
        while(prev_fq_size > 0 && trim_times < TRIM_TIMES)
        {
            time = wtime();
            trim_1_from_fq_gfq(scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    thread_count,
                    small_queue,
                    thread_bin,
                    prefix_sum,
                    tid,
                    temp_queue);
            #pragma omp barrier
            
            if(tid == 0)
            {
                time_size_1 += wtime() - time;
            }

            fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("trim_1, 3rd, %.3lf, trimmed vertices, %d\n", time_size_1 * 1000, prev_fq_size - fq_size);
                }
            }
            step = fq_size / thread_count;
            vert_beg = tid * step;
            vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
            if(prev_fq_size - fq_size < prev_fq_size * theta)
                break;
            prev_fq_size = fq_size;
            trim_times += 1;
            #pragma omp barrier
        }
        if(DEBUG)
        {
            if(tid == 0)
                printf("trim-1 times before trim-2-3, %d\n", trim_times);
        }

        if(fq_size > 0)
        {
     // not using trim-2
            fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            step = fq_size / thread_count;
            vert_beg = tid * step;
            vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);

            if(DEBUG)
            {
                if(tid == 0)
                    printf("fq_size, %d\n", fq_size);
            }
            time = wtime();
            //step 3.2: trim size_2
    //        trim_2_from_graph(scc_id,
    //                fw_beg_pos,
    //                bw_beg_pos,
    //                vert_beg,
    //                vert_end,
    //                fw_csr,
    //                bw_csr);
            trim_2_from_fq(scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    small_queue);
            #pragma omp barrier
            if(tid == 0)
                time_size_2 = wtime() - time;
            
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("time_size_2, %.3lf, fq_size, %d\n", time_size_2 * 1000, fq_size);
                }
            }
//note! just for record            
            //time_size_1 = 0;
            trim_times = 0;
            #pragma omp barrier
            
            
// not using trim-3
            //step 3.3: trim size-3
            //size-3: type_1, A --> B --> C --> A
            time = wtime();
    //        trim_3_1_from_graph(scc_id,
    //                fw_beg_pos,
    //                bw_beg_pos,
    //                vert_beg,
    //                vert_end,
    //                fw_csr,
    //                bw_csr);

            trim_3_1_from_fq(scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    small_queue);
            #pragma omp barrier
            double time_size_3_1 = wtime() - time;

            if(DEBUG)
            {
                if(tid == 0)
                    printf("time_size_3_1, %.3lf\n", time_size_3_1 * 1000);
            }

            //step 3.3: trim size-3
            //size-3: type_1, A --> B --> A --> C --> A 
            //starting from hub vertex A
            
            time = wtime();
    //        trim_3_2_from_graph(scc_id,
    //                fw_beg_pos,
    //                bw_beg_pos,
    //                vert_beg,
    //                vert_end,
    //                fw_csr,
    //                bw_csr);
            trim_3_2_from_fq(scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    small_queue);
            
            trim_times = 0;
            
/// note!
//            time_size_1 = 0;
            #pragma omp barrier
            double time_size_3_2 = wtime() - time;

            if(DEBUG)
            {
                if(tid == 0)
                    printf("time_size_3_2, %.3lf\n", time_size_3_2 * 1000);
            }
            time_size_3 = time_size_3_1 + time_size_3_2;
            
            
            //step 3.4: trim size_1
            time = wtime();
    //        trim_1_normal(scc_id,
    //                fw_beg_pos,
    //                bw_beg_pos,
    //                vert_beg,
    //                vert_end,
    //                fw_csr,
    //                bw_csr);
    //        trim_1_from_fq(scc_id,
    //                fw_beg_pos,
    //                bw_beg_pos,
    //                vert_beg,
    //                vert_end,
    //                fw_csr,
    //                bw_csr,
    //                small_queue);
    //        #pragma omp barrier
    //        time_size_1 += wtime() - time;
    //        
    //        ///step 4.0 generate frontier queue
    //        time = wtime();
    //        gfq_from_queue(vert_count,
    //                scc_id,
    //                thread_count,
    //                small_queue,
    //                thread_bin,
    //                prefix_sum,
    //                vert_beg,
    //                vert_end,
    //                tid,
    //                temp_queue);
    //        double gfq_time = wtime() - time;
    //        time_gfq += gfq_time;
                
            prev_fq_size = fq_size;
            while(prev_fq_size > 0)
            {
                time = wtime();
                trim_1_from_fq_gfq(scc_id,
                        fw_beg_pos,
                        bw_beg_pos,
                        vert_beg,
                        vert_end,
                        fw_csr,
                        bw_csr,
                        thread_count,
                        small_queue,
                        thread_bin,
                        prefix_sum,
                        tid,
                        temp_queue);
                #pragma omp barrier
                if(tid == 0)
                    time_size_1 += wtime() - time;

                fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
                step = fq_size / thread_count;
                vert_beg = tid * step;
                vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
                if(DEBUG)
                {
                    if(tid == 0)
                    {
                        printf("trim_1, 4th, %.3lf, fq_size, %d, trimmed vertices, %d\n", time_size_1 * 1000, fq_size, prev_fq_size - fq_size);
                    }
                }
                if(prev_fq_size - fq_size < prev_fq_size * theta)
                    break;
                prev_fq_size = fq_size;
                #pragma omp barrier
                
                ///trim-1 only once
                break;
            }
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("trim_1 times before Color, %d, fq_size, %d\n", trim_times, fq_size);
                }
            }
          
//            generate_frontier_queue(vert_count,
//                    scc_id,
//                    thread_count,
//                    small_queue,
//                    thread_bin,
//                    prefix_sum,
//                    vert_beg,
//                    vert_end,
//                    tid);
//            #pragma omp barrier
//            index_t fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
//            step = fq_size / thread_count;
//            vert_beg = tid * step;
//            vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
//

/*
            time = wtime();
            gfq_from_queue(vert_count,
                    scc_id,
                    thread_count,
                    small_queue,
                    thread_bin,
                    prefix_sum,
                    vert_beg,
                    vert_end,
                    tid,
                    temp_queue);
            #pragma omp barrier
            index_t fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            step = fq_size / thread_count;
            vert_beg = tid * step;
            vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
*/
//            if(tid == 0)
//                wcc_color[small_queue[vert_beg]] = 0;
            #pragma omp barrier
//            index_t queue_size;
//            if(fq_size * thread_count > vert_count)
//            {
//                queue_size = vert_count / thread_count;
//            }
//            else
//            {
//                queue_size = fq_size;
//            }
//    
//            index_t *bw_beg_pos = g->bw_beg_pos;
//    vertex_t *bw_csr = g->bw_csr;

// try using WCC + FW-BW

            time = wtime();
// color WCC, with both FW & BW
            coloring_wcc(fq_size,
                    scc_id,
                    thread_count,
                    small_queue,
                    vert_beg,
                    vert_end,
                    tid,
                    color,
                    color_change,
                    fw_beg_pos,
                    fw_csr,
                    bw_beg_pos,
                    bw_csr);
/*
// Hang's label propragation WCC
            openmp_wcc<vertex_t, index_t, color_t>
                (fw_csr,//change
                 fw_beg_pos,//change
                 wcc_color,//wcc_color[]
                 color_redirect,
                 is_redirect,
                 color_change,
                 global_color,
                 tid,
                 thread_count,
                 vert_beg,
                 vert_end,
                 vert_count,
                 edge_count,
                 scc_id,
                 small_queue);
                 */
            #pragma omp barrier
            if(tid == 0)
            {
                time_wcc = wtime() - time;
                //printf("WCC time (ms), %.3lf\n", time_wcc * 1000);
//                printf("global color, %d\n", global_color[0]);
//                for(int i=0; i<fq_size; ++i)
//                {
//                    printf("%d, %d\n", small_queue[i], wcc_color[small_queue[i]]);
//                }
            }

            time = wtime();
// using fw-bw to detect mice SCCs
            init_fw_sa(vert_beg,
                    vert_end,
                    fw_sa,
                    small_queue,
                    fq_size,
                    wcc_fq,
                    tid,
                    color,
                    wcc_fq_size);
            #pragma omp barrier
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("wcc_fq, %d\n", wcc_fq_size);
                }
            }
//            if(tid == 0)
//            {
            mice_fw_bw(color,
                    scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    fw_csr,
                    bw_csr,
                    fw_sa,
                    tid,
                    thread_count,
                    small_queue,
                    fq_size,
                    wcc_fq,
                    wcc_fq_size
                    ); 
//            }
            #pragma omp barrier
            if(tid == 0)
            {
                time_mice_fw_bw = wtime() - time;
//                printf("FW-BW mice (ms), %.3lf\n", time_mice_fw_bw * 1000);
            }
            

// try sorting the frontier queue
//            if(tid == 0)
//                printf("queue_size, %d\n", queue_size);
//
//            double time_sort = wtime();
//            if(tid == 0)
//            {
//                std::sort(small_queue, small_queue + queue_size);
//            }
//
//            #pragma omp barrier
//            time_sort = wtime() - time_sort;
//            if(DEBUG && tid == 0)
//            {
//                printf("sort time, %g\n", time_sort * 1000);
//            }

/** coloring algorithm
            index_t *q = new index_t[queue_size];
            double time_rank = wtime();
            
            degree_rank(fq_size,
                scc_id,
                thread_count,
                small_queue,
                vert_beg,
                vert_end,
                tid,
                mul_degree,
                degree_prop,
                fw_beg_pos,
                fw_csr,
                bw_beg_pos,
                bw_csr);
            #pragma omp barrier
            time_rank = wtime() - time_rank;
            if(DEBUG && tid == 0)
            {
                printf("degree rank time, %g\n", time_rank * 1000);
            }

            time_color = wtime();
            time_color_1 = 0;
            time_color_2 = 0;
            time_color_init = 0;
            
            graph_color(fq_size,
                    scc_id,
                    thread_count,
                    small_queue,
                    vert_beg,
                    vert_end,
                    tid,
                    color,
                    color_change,
                    time_color_1,
                    time_color_2,
                    time_color_init,
                    bw_beg_pos,
                    bw_csr,
                    q,
                    mul_degree,
                    degree_prop,
                    color_times);
                    
//            printf("tid, %d\n", tid);
            #pragma omp barrier
            time_color = wtime() - time_color;
            if(DEBUG)
            {
                if(tid == 0)
                    printf("tid, %d, color time, %.3lf\n", tid, time_color * 1000);
            }
//            if(DEBUG)
//            {
//                if(tid == 0)
//                {
//                    vertex_t color_total = 0;
//                    for(index_t i=0; i<vert_count; ++i)
//                    {
//                        if(color_times[i] > 0)
//                            color_total += color_times[i];
//                    }
//                    printf("total color times, %d\n", color_total);
//                }
//            }
//            printf("q[0] = %d\n", q[0]);
            delete[] q;
            */
        }
        if(tid == 0)
        {
            avg_time[0] += time_size_1_first + time_size_1 + time_size_2 + time_size_3;
//            avg_time[0] += time_size_1;
            avg_time[1] += time_fw + time_bw + pivot_time;
            avg_time[2] += time_wcc + time_mice_fw_bw;
            avg_time[4] += time_size_1_first + time_size_1;
//            avg_time[4] += time_size_1;
            avg_time[5] += time_size_2;
            avg_time[6] += pivot_time;
            avg_time[7] += time_fw;
            avg_time[8] += time_bw;
            avg_time[9] += time_wcc;
            avg_time[10] += time_mice_fw_bw;
//            avg_time[11] += time_color_init;
            
            avg_time[13] += time_size_3;
            avg_time[14] += time_gfq;
            
//            if(DEBUG)
//            {
//                for(int i=0; i< 15; ++i)
//                    printf("%.3lf\n", avg_time[i]);
//                printf("algorithm run finishes\n");
//            }
        }
        if(OUTPUT_TIME)
        {
            if(tid == 0)
            {
                printf("\ntime size_1_first, %.3lf\ntime size_1, %.3lf\ntime pivot, %.3lf\nlargest fw, %.3lf\nlargest bw, %.3lf\nlargest fw/bw, %.3lf\ntrim size_2, %.3lf\ntrim size_3, %.3lf\nwcc time, %.3lf\nmice fw-bw time, %.3lf\nmice scc time, %.3lf\ntotal time, %.3lf\n", time_size_1_first * 1000, time_size_1 * 1000, pivot_time * 1000, time_fw * 1000, time_bw * 1000, (pivot_time + time_fw + time_bw) * 1000, time_size_2 * 1000, time_size_3 * 1000, time_wcc * 1000, time_mice_fw_bw * 1000, (time_wcc + time_mice_fw_bw) * 1000, (time_size_1_first + time_size_1 + pivot_time + time_fw + time_bw + time_size_2 + time_size_3 + time_wcc + time_mice_fw_bw) * 1000);
            }
        }
        #pragma omp barrier
    }
    end_time = wtime() - start_time;
    avg_time[3] += end_time;
    if(DEBUG)
        printf("total time, %.3lf\n", end_time * 1000);

    get_scc_result(scc_id,
            vert_count);

    delete[] scc_id;
    delete[] color;
    delete[] max_pivot_list;
    delete[] max_degree_list;
    delete[] small_queue;
    delete[] temp_queue;
    delete[] thread_bin;
    delete[] prefix_sum;
    delete[] front_comm;
	delete[] work_comm;
    delete[] color_change;
//    delete[] wcc_color;
//    delete[] color_redirect;
//    delete[] is_redirect;

}


