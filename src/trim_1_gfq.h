#ifndef TRIM_1_GFQ_H
#define TRIM_1_GFQ_H

#include "scc_common.h"
#include "util.h"
#include "wtime.h"
#include <omp.h>

inline void trim_1_first(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end
        )
{
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id] == 0)
        {
            scc_id[vert_id] = -1;
            continue;
        }
        else
        {
            if(bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id] == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
        }
    }
}

inline void trim_1_first_gfq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid
        )
{
    //step 1: get thread_bin_size during trim size-1
    index_t thread_bin_size = 0;
    
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id] == 0)
        {
            scc_id[vert_id] = -1;
            continue;
        }
        else
            if(bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id] == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            else
            {
                thread_bin_size ++;
            }
    }
    
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }
    
    #pragma omp barrier

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
}

//step 1.2: trim size_1
inline void trim_1_normal(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr
        )
{
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;

            }
        }
    }
}

inline void trim_1_from_fq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        index_t *small_queue
        )
{
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;

            }
        }
    }
}

// return the number of trimmed vertices
inline void trim_1_normal_only_size(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        const index_t thread_count,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid
        )
{
    //step 1: get thread_bin_size during trim size-1
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                thread_bin_size ++;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;
                thread_bin_size ++;
                continue;
            }
            
        }
    }
    
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

//    //step 3: write the vertices into fq
//    vertex_t start_pos = prefix_sum[tid];
//    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
//    {
//        if(scc_id[vert_id] == 0)
//        {
//            frontier_queue[start_pos++] = vert_id;
//        }
//    }
}

//step 1.2: trim size_1
inline void trim_1_normal_gfq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid
        )
{
    //step 1: get thread_bin_size during trim size-1
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            
            thread_bin_size ++;
        }
    }
    
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
}

inline void trim_1_from_fq_gfq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    out_degree = 1;
                    break;
                }
            }
            if(out_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    in_degree = 1;
                    break;
                }
            }
            if(in_degree == 0)
            {
                scc_id[vert_id] = -1;
                continue;
            }
            thread_bin_size++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = frontier_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
    #pragma omp barrier
    if(DEBUG)
    {
        if(tid == 0)
        {
            printf("In normal trim, thread bin size, %d\n", prefix_sum[55]);
        }
    }
    //step 4: write back to small_queue
    for(index_t i=prefix_sum[tid]; i<prefix_sum[tid]+thread_bin[tid]; ++i)
    {
        frontier_queue[i] = temp_queue[i];
    }
}

inline static void get_queue(
        vertex_t *thread_queue,
        vertex_t *thread_bin,
        index_t *prefix_sum,
        index_t tid,
        vertex_t *temp_queue
        )
{
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }
    #pragma omp barrier
    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = start_pos; vert_id < start_pos + thread_bin[tid]; ++vert_id)
    {
        temp_queue[vert_id] = thread_queue[vert_id-start_pos];
    }
}

// Using prefix sum to generate frontier queue 
inline static void generate_frontier_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
//    #pragma omp barrier
//    return prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

}

inline static void gfq_from_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
    #pragma omp barrier
    //step 4: write back to small_queue
    for(index_t i=prefix_sum[tid]; i<prefix_sum[tid]+thread_bin[tid]; ++i)
    {
        small_queue[i] = temp_queue[i];
    }

}

inline static void bw_gfq_from_fw(
        index_t *fw_sa,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(fw_sa[vert_id] != -1)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    #pragma omp barrier
    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(fw_sa[vert_id] != -1)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
//    #pragma omp barrier
//    //step 4: write back to small_queue
//    for(index_t i=prefix_sum[tid]; i<prefix_sum[tid]+thread_bin[tid]; ++i)
//    {
//        small_queue[i] = temp_queue[i];
//    }

}

inline static void gfq_fw_bw_from_queue(
        index_t *sa,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(sa[vert_id] == -1)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(sa[vert_id] == -1)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
}

#endif


