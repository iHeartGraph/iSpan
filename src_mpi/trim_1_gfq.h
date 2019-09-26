#ifndef TRIM_1_GFQ_H
#define TRIM_1_GFQ_H

#include "scc_common.h"
#include "util.h"
#include "wtime.h"
#include <omp.h>

inline void trim_1_first(
        index_t *scc_id,
        long_t *fw_beg_pos,
        long_t *bw_beg_pos,
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
        long_t *fw_beg_pos,
        long_t *bw_beg_pos,
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
            long_t my_beg = fw_beg_pos[vert_id];
            long_t my_end = fw_beg_pos[vert_id+1];
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

inline static void gfq_origin(
        const index_t vert_count,
        index_t *scc_id,
        index_t *frontier_queue,
        vertex_t vert_beg,
        vertex_t vert_end,
        long_t *fw_beg_pos,
        vertex_t *fw_csr,
        long_t *bw_beg_pos,
        vertex_t *bw_csr,
        vertex_t *sub_fw_beg,
        vertex_t *sub_fw_csr,
        vertex_t *sub_bw_beg,
        vertex_t *sub_bw_csr,
        vertex_t *front_comm,
        long_t *work_comm,
        vertex_t world_rank,
        vertex_t world_size,
        vertex_t *vert_map
        )
{
    index_t index = 0;
    index_t fw_edge_num = 0;
    index_t bw_edge_num = 0;
    for(vertex_t vert_id = 0; vert_id < vert_count; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            frontier_queue[index] = vert_id;
            vert_map[vert_id] = index;
            index++;
        }
    }

    for(vertex_t v_index = 0; v_index < index; ++v_index)
    {
        vertex_t vert_id = frontier_queue[v_index];

//    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
//    {
//        vertex_t vert_id = frontier_queue[v_index];
        //Forward CSR
        long_t my_beg = fw_beg_pos[vert_id];
        long_t my_end = fw_beg_pos[vert_id + 1];
        for(long_t w_index = my_beg; w_index < my_end; ++w_index)
        {
            vertex_t w = fw_csr[w_index];
            if(scc_id[w] == 0)
            {
                sub_fw_csr[fw_edge_num] = vert_map[w];
                fw_edge_num++;
            }
        }
        sub_fw_beg[v_index + 1] = fw_edge_num;

        //Backward CSR
        my_beg = bw_beg_pos[vert_id];
        my_end = bw_beg_pos[vert_id + 1];
        for(long_t w_index = my_beg; w_index < my_end; ++w_index)
        {
            vertex_t w = bw_csr[w_index];
            if(scc_id[w] == 0)
            {
                sub_bw_csr[bw_edge_num] = vert_map[w];
                bw_edge_num++;
            }
        }
        sub_bw_beg[v_index + 1] = bw_edge_num;
    //    }
    }
    
    front_comm[world_rank] = index; //Vertex count
    work_comm[world_rank] = fw_edge_num; //Edge count
    std::cout<<"sub v_count, "<<front_comm[world_rank]<<", sub e_count, "<<work_comm[world_rank]<<","<<bw_edge_num<<"\n";

//    return front_comm[world_rank];
}
//inline static void gfq_distribute(
//        const index_t vert_count,
//        index_t *scc_id,
//        index_t *frontier_queue,
//        vertex_t vert_beg,
//        vertex_t vert_end,
//        index_t *fw_beg_pos,
//        vertex_t *fw_csr,
//        index_t *bw_beg_pos,
//        vertex_t *bw_csr,
//        vertex_t *sub_fw_beg,
//        vertex_t *sub_fw_csr,
//        vertex_t *sub_bw_beg,
//        vertex_t *sub_bw_csr,
//        vertex_t *buf_fw_beg,
//        vertex_t *buf_fw_csr,
//        vertex_t *buf_bw_beg,
//        vertex_t *buf_bw_csr,
//        int *size_fw_beg,
//        int *size_fw_csr,
//        int *size_bw_beg,
//        int *size_bw_csr,
//        vertex_t world_rank,
//        vertex_t world_size,
//        vertex_t *vert_map
//        )
//{
//    index_t index = 0;
//    index_t fw_edge_num = 0;
//    index_t bw_edge_num = 0;
//    for(vertex_t vert_id = 0; vert_id < vert_count; ++vert_id)
//    {
//        if(scc_id[vert_id] == 0)
//        {
//            frontier_queue[index] = vert_id;
//            vert_map[vert_id] = index;
//            index++;
//        }
//    }
//
////    front_comm[world_rank] = index; //Vertex count
////    work_comm[world_rank] = fw_edge_num; //Edge count
////    work_comm_bw[world_rank] = bw_edge_num;
//
////    for(vertex_t v_index = 0; v_index < index; ++v_index)
////        vertex_t vert_id = frontier_queue[v_index];
//    
//    vertex_t num_fw_beg = 0;
////    vertex_t virtual_id = 0;
//    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
//    {
////        vertex_t vert_id = frontier_queue[v_index];
//        //Forward CSR
//        if(scc_id[vert_id] != 0)
//            continue;
//        num_fw_beg += 1;
//        index_t my_beg = fw_beg_pos[vert_id];
//        index_t my_end = fw_beg_pos[vert_id + 1];
//        for(vertex_t w_index = my_beg; w_index < my_end; ++w_index)
//        {
//            vertex_t w = fw_csr[w_index];
//            if(scc_id[w] == 0)
//            {
//                buf_fw_csr[fw_edge_num] = vert_map[w];
//                fw_edge_num++;
//            }
//        }
//        buf_fw_beg[num_fw_beg] = fw_edge_num;
////        buf_fw_beg[vert_map[vert_id] + 1] = fw_edge_num;
//
//        //Backward CSR
//        my_beg = bw_beg_pos[vert_id];
//        my_end = bw_beg_pos[vert_id + 1];
//        for(vertex_t w_index = my_beg; w_index < my_end; ++w_index)
//        {
//            vertex_t w = bw_csr[w_index];
//            if(scc_id[w] == 0)
//            {
//                buf_bw_csr[bw_edge_num] = vert_map[w];
//                bw_edge_num++;
//            }
//        }
////        buf_bw_beg[vert_map[vert_id] + 1] = bw_edge_num;
//        buf_bw_beg[num_fw_beg] = bw_edge_num;
//    }
//    size_fw_beg[world_rank] = num_fw_beg;
//    size_bw_beg[world_rank] = num_fw_beg;
//    size_fw_csr[world_rank] = fw_edge_num;
//    size_bw_csr[world_rank] = bw_edge_num;
//    
//    double temp_time = wtime();
//
//    MPI_Allgather(&size_fw_beg[world_rank],
//        1,
//        MPI_LONG,
//        size_fw_beg,
//        1,
//        MPI_LONG,
//        MPI_COMM_WORLD);
//
//    MPI_Allgather(&size_bw_beg[world_rank],
//        1,
//        MPI_LONG,
//        size_bw_beg,
//        1,
//        MPI_LONG,
//        MPI_COMM_WORLD);
//
//    MPI_Allgather(&size_fw_csr[world_rank],
//        1,
//        MPI_LONG,
//        size_fw_csr,
//        1,
//        MPI_LONG,
//        MPI_COMM_WORLD);
//
//    MPI_Allgather(&size_bw_csr[world_rank],
//        1,
//        MPI_LONG,
//        size_bw_csr,
//        1,
//        MPI_LONG,
//        MPI_COMM_WORLD);
//
//    double comm_time = wtime() - temp_time;
//    
////    if(world_rank == 0)
////    vertex_t e_num_fw = 0;
//    vertex_t e_num = 0;
//    vertex_t v_num = 0;
//    vertex_t *beg_fw_beg = (vertex_t *)calloc(world_size, sizeof(int));
//    vertex_t *beg_bw_beg = (vertex_t *)calloc(world_size, sizeof(int));
//    vertex_t *beg_fw_csr = (vertex_t *)calloc(world_size, sizeof(int));
//    vertex_t *beg_bw_csr = (vertex_t *)calloc(world_size, sizeof(int));
//
//    for(int i = 0; i < world_size - 1; ++i)
//    {
//        v_num += size_fw_beg[i];
//        e_num += size_fw_csr[i];
//        beg_fw_beg[i + 1] = beg_fw_beg[i] + size_fw_beg[i];
//        beg_bw_beg[i + 1] = beg_bw_beg[i] + size_bw_beg[i];
//        beg_fw_csr[i + 1] = beg_fw_csr[i] + size_fw_csr[i];
//        beg_bw_csr[i + 1] = beg_bw_csr[i] + size_bw_csr[i];
//    }
//    v_num += size_fw_beg[world_size - 1];
//    e_num += size_fw_csr[world_size - 1];
//
//
////    std::cout<<"e_num_fw,"<<e_num_fw<<","<<"e_num_bw,"<<e_num_bw<<"\n";
//    
//    // Accumulating the begin position
////    vertex_t part_id = 0;
////    vertex_t fw_edge_start = 0;
//
//    for(int i = 0; i < num_fw_beg; ++i)
//    {
//        buf_fw_beg[i] += beg_fw_csr[world_rank];
//        buf_bw_beg[i] += beg_bw_csr[world_rank];
//    }
////    if(world_rank == 0)
////    {
////        for(int i = 0; i < num_fw_beg; ++i)
////        {
////            std::cout<<" "<<i<<":"<<buf_bw_beg[i]<<",";
////        }
//
////    }
////    std::cout<<"\n";
////    for(int i = 0; i < world_size; ++i)
////    {
//////        std::cout<<"size array,"<<size_fw_beg[i]<<","<<size_fw_csr[i]<<","<<size_bw_csr[i]<<"\n";
////        std::cout<<"sub_bw,"<<i<<","<<sub_bw_beg[i]<<",size_bw_beg,"<<size_bw_beg[i]<<",beg_bw_beg,"<<beg_bw_beg[i]<<"\n";
////    }
//    
//
////    vertex_t * displs = (vertex_t *)calloc(world_size, sizeof(int));
//    
//    temp_time = wtime();
//    
//    MPI_Allgatherv(buf_fw_beg,
//        num_fw_beg,
//        MPI_INT,
//        sub_fw_beg,
//        size_fw_beg,
//        beg_fw_beg,
//        MPI_INT,
//        MPI_COMM_WORLD);
//
//    MPI_Allgatherv(buf_fw_csr,
//        fw_edge_num,
//        MPI_INT,
//        sub_fw_csr,
//        size_fw_csr,
//        beg_fw_csr,
//        MPI_INT,
//        MPI_COMM_WORLD);
//
//    MPI_Allgatherv(buf_bw_beg,
//        num_fw_beg,
//        MPI_INT,
//        sub_bw_beg,
//        size_bw_beg,
//        beg_bw_beg,
//        MPI_INT,
//        MPI_COMM_WORLD);
//
//    MPI_Allgatherv(buf_bw_csr,
//        bw_edge_num,
//        MPI_INT,
//        sub_bw_csr,
//        size_bw_csr,
//        beg_bw_csr,
//        MPI_INT,
//        MPI_COMM_WORLD);
//
//    comm_time += wtime() - temp_time;
//
//    std::cout<<"gfq_comm time," << comm_time * 1000 << " (ms)\n";
//    
////    if(world_rank == 0)
////    {
//////        for(int i = 0; i < v_num; ++i)
//////        {
//////            std::cout<<" "<<i<<":"<<sub_fw_beg[i]<<",";
//////        }
//////        std::cout<<"\n";
//////        for(int i = 0; i < e_num; ++i)
//////        {
//////            std::cout<<" "<<i<<":"<<sub_fw_csr[i]<<",";
//////        }
//////        std::cout<<"\n";
//////
////        for(int i = 0; i < v_num; ++i)
////        {
////            std::cout<<" "<<i<<":"<<sub_bw_beg[i]<<",";
////        }
////        std::cout<<"\n";
//////
//////        for(int i = 0; i < e_num; ++i)
//////        {
//////            std::cout<<" "<<i<<":"<<sub_bw_csr[i]<<",";
//////        }
//////        std::cout<<"\n";
////    }
//
//    
////    std::cout<<"sub v_count, "<<front_comm[world_rank]<<", sub e_count, "<<work_comm[world_rank]<<","<<bw_edge_num<<"\n";
//
//}

#endif


