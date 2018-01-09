#ifndef COLOR_PROPAGATION_H
#define COLOR_PROPAGATION_H

#include "graph.h"
#include "util.h"
#include "wtime.h"
#include <set>

inline static void degree_rank(
        const index_t fq_size,
        index_t *scc_id,
        const index_t thread_count,
        index_t *small_queue,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *mul_degree,
        index_t *degree_prop,
        index_t *fw_beg_pos,
        index_t *fw_csr,
        index_t *bw_beg_pos,
        index_t *bw_csr
        )
{
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++ fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            /// fw, out degree
            index_t out_degree = 0;
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];

            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0)
                    out_degree ++;
            }
            
            /// bw, in degree
            index_t in_degree = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
//            printf("%d, %d, %d ", vert_id, my_beg, my_end);
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0)
                    in_degree ++;
            }
//            printf("%d, %d, %d\n", vert_id, out_degree, in_degree);
//          fisc used
//            mul_degree[vert_id] = out_degree * (in_degree + 1);
            mul_degree[vert_id] = (out_degree + 5) * (in_degree);
//            mul_degree[vert_id] = 0;
            degree_prop[vert_id] = mul_degree[vert_id];
        }
    }

}

inline static void color_propagation(
        const index_t fq_size,
        index_t *scc_id,
        const index_t thread_count,
        index_t *small_queue,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *color,
        bool *color_change,
//        index_t *fw_beg_pos,
//        index_t *fw_csr,
        index_t *bw_beg_pos,
        index_t *bw_csr,
        index_t *mul_degree,
        index_t *degree_prop,
        index_t *color_times
        )
{
    index_t depth = 0;
    while(true)
    {
        #pragma omp barrier
        if(DEBUG)
        {
            depth += 1;
        }
        
        bool color_changed = false;
        //option 1: bottom up
        for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++ fq_vert_id)
        {
            vertex_t vert_id = small_queue[fq_vert_id];
            if(scc_id[vert_id] == 0)
            {
                index_t my_beg = bw_beg_pos[vert_id];
                index_t my_end = bw_beg_pos[vert_id+1];

                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    if(vert_id == w)
                        continue;
//                    if(scc_id[w] == 0 && color[w] > color[vert_id])
                    if(scc_id[w] == 0)
                    {
                        if(degree_prop[vert_id] < degree_prop[w])
                        {
                            degree_prop[vert_id] = degree_prop[w];
                            color[vert_id] = color[w];
                            color_times[vert_id] += 1;
                            if(!color_changed)
                                color_changed = true;
                        }
                        else
                            if(degree_prop[vert_id] == degree_prop[w])
                            {
                                if(color[vert_id] < color[w])
                                {
                                    color[vert_id] = color[w];
                                    color_times[vert_id] += 1;
                                    if(!color_changed)
                                        color_changed = true;
                                }
                            }
                    }
                }
            }
        }
        #pragma omp barrier
        
//        printf("depth, %d, tid, %d\n", depth, tid);
        // path compression                
        if(color_changed)
        {
            for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
            {
                vertex_t vert_id = small_queue[fq_vert_id];
//                if(DEBUG)
//                {
//                    printf("tid, %d, vert_id, %d, color, %d, degree, %d\n", tid, vert_id, color[vert_id], degree_prop[vert_id]);
//                }
//                printf("tid, %d, vert_id, %d\n", tid, vert_id);
                if(scc_id[vert_id] == 0 && color[vert_id] != vert_id)
                {
                    index_t root = color[vert_id];
//                    index_t root = vert_id;
                    index_t c_depth = 0;
                    while(color[root] != root && c_depth < 100)// && degree_prop[root] <= degree_prop[color[root]])
                    {
//                        printf("tid, %d, vert_id, %d, root, %d, color[root], %d, degree, %d\n", tid, vert_id, root, color[root], degree_prop[root]);
                        root = color[root];
                        c_depth ++;
                    }
                    index_t v_id = vert_id;
                    while(v_id != root && color[v_id] != root)// && degree_prop[v_id] < degree_prop[root])
                    {
                        index_t prev = v_id;
                        v_id = color[v_id];
                        color[prev] = root;
                        degree_prop[prev] = degree_prop[root];
                        color_times[prev] += 1;
//                                depth ++;
                    }
//                            printf("%d\n", depth);
//                            if(root != color[root])
//                            {
//                                color[vert_id] = color[root];
//                            }
                }
            }
            color_change[tid] = true;
        }
        else
        {
            color_change[tid] = false;
        }
        #pragma omp barrier
//        printf("after path compression, %d\n", tid);

        bool final_color_change = false;
        for(index_t i=0; i<thread_count; ++i)
        {
            if(color_change[i])
            {
                final_color_change = true;
                break;
            }
        }
        #pragma omp barrier
        if(final_color_change == false)
        {
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("iteration_depth, %d, ", depth);
                }
            }
            return;
        }
    }
}


//setp 4.2: color identify, first select root, then top down bfs
inline static void color_identify(
        const index_t fq_size,
        index_t *scc_id,
        const index_t thread_count,
        index_t *small_queue,
        index_t vert_beg,
        index_t vert_end,
        index_t *color,
        index_t *bw_beg_pos,
        index_t *bw_csr,
        index_t *q,
        index_t *mul_degree,
        index_t *degree_prop
        )
{
    
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++ fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0 && color[vert_id] == vert_id)
        {
            scc_id[vert_id] = vert_id;
            ///top down bfs, queue based, no need to use status array
            index_t head = 0;
            index_t tail = 0;
            q[tail++] = vert_id;
            if(tail == fq_size)
                tail = 0;

            while(head != tail)
            {
                vertex_t temp_v = q[head++];
                if(head == fq_size)
                    head = 0;
                index_t my_beg = bw_beg_pos[temp_v];
                index_t my_end = bw_beg_pos[temp_v+1];

                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    
                    if(scc_id[w] == 0 && color[w] == vert_id)
                    {
                        scc_id[w] = vert_id;
                        q[tail++] = w;
                        if(tail == fq_size)
                            tail = 0;
                    }
                }
            }
        }
    }
}

inline static void color_init(
        index_t *scc_id,
        index_t *small_queue,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *color,
        bool *color_change,
        index_t *mul_degree,
        index_t *degree_prop,
        index_t *color_times
        )
{
    bool color_changed = false;
   
//    index_t sum = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++ fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            color[vert_id] = vert_id;
            degree_prop[vert_id] = mul_degree[vert_id];
//            sum ++;
            if(!color_changed)
                color_changed = true;
            color_times[vert_id] += 1;
        }
    }
//    #pragma omp barrier
    color_change[tid] = color_changed;
//    printf("%d, %d\n", tid, sum);
}

inline static void color_statistic(index_t *scc_id,
                index_t *small_queue,
                index_t vert_beg,
                index_t vert_end,
                index_t tid,
                index_t *color,
                index_t *color_times,
                index_t fq_size)
{
    index_t scc_num = 0;
    std::set<int> statis_set;
    vertex_t color_total = 0;
    if(tid == 0)
    {
        for(vertex_t fq_vert_id = 0; fq_vert_id < fq_size; ++ fq_vert_id)
        {
            vertex_t vert_id = small_queue[fq_vert_id];
            if(scc_id[vert_id] != 0)
            {
                scc_num++;
                if(statis_set.find(color[vert_id]) == statis_set.end())
                {
                    statis_set.insert(color[vert_id]);
                }
            }
        }
        for(index_t i=0; i<fq_size; ++i)
        {
            vertex_t vert_id = small_queue[i];
            if(color_times[vert_id] > 0)
                color_total += color_times[vert_id];
        }
        printf("detected vertex number, %d, scc_num, %d, color times, %d\n", scc_num, statis_set.size(), color_total);
    }
}

inline static void graph_color(
        const index_t fq_size,
        index_t *scc_id,
        const index_t thread_count,
        index_t *small_queue,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *color,
        bool *color_change,
        double &time_color_1,
        double &time_color_2,
        double &time_color_init,
//        index_t *fw_beg_pos,
//        index_t *fw_csr,
        index_t *bw_beg_pos,
        index_t *bw_csr,
        index_t *q,
        index_t *mul_degree,
        index_t *degree_prop,
        index_t *color_times
        )
{
    index_t round_num = 0;
    while(true)
    {
        #pragma omp barrier
        double time = wtime();
        color_propagation(fq_size,
                scc_id,
                thread_count,
                small_queue,
                vert_beg,
                vert_end,
                tid,
                color,
                color_change,
//                fw_beg_pos,
//                fw_csr,
                bw_beg_pos,
                bw_csr,
                mul_degree,
                degree_prop,
                color_times);
        #pragma omp barrier
//        if(DEBUG && tid == 0)
//        {
//            printf("color propagation finishes\n");
//        }
        time_color_1 += wtime() - time;

        time = wtime();
        color_identify(fq_size,
                scc_id,
                thread_count,
                small_queue,
                vert_beg,
                vert_end,
                color,
                bw_beg_pos,
                bw_csr,
                q,
                mul_degree,
                degree_prop);
        #pragma omp barrier
        time_color_2 += wtime() - time;
        
//        if(DEBUG && tid == 0)
//        {
//            printf("color identify finishes\n");
//        }
        if(DEBUG)
        {
            color_statistic(scc_id,
                    small_queue,
                    vert_beg,
                    vert_end,
                    tid,
                    color,
                    color_times,
                    fq_size);
        }
        #pragma omp barrier
       
        time = wtime();
        color_init(scc_id,
                small_queue,
                vert_beg,
                vert_end,
                tid,
                color,
                color_change,
                mul_degree,
                degree_prop,
                color_times);
        #pragma omp barrier
        time_color_init += wtime() - time;

//        if(DEBUG && tid == 0)
//        {
//            printf("color init finishes\n");
//        }
        //check whether it finishes
        bool final_color_change = false;
        for(index_t i=0; i<thread_count; ++i)
        {
            if(color_change[i])
            {
//                if(DEBUG && tid == 0)
//                    printf("tid, %d\n", i);
                final_color_change = true;
                break;
            }
        }
        #pragma omp barrier
        if(DEBUG)
        {
            round_num ++;
//            if(tid == 0)
//                printf("%d\n", round_num);
        }
        #pragma omp barrier
        if(final_color_change == false)
        {
            if(DEBUG)
            {
                if(tid == 0)
                    printf("round num, %d\n", round_num);
            }
            return;
        }
    }
//    #pragma omp barrier
}
inline static void coloring_wcc(
        const index_t fq_size,
        index_t *scc_id,
        const index_t thread_count,
        index_t *small_queue,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *color,
        bool *color_change,
        index_t *fw_beg_pos,
        index_t *fw_csr,
        index_t *bw_beg_pos,
        index_t *bw_csr
        )
{
    index_t depth = 0;
    while(true)
    {
        #pragma omp barrier
        if(DEBUG)
        {
            depth += 1;
        }
        
        bool color_changed = false;
        //option 1: bottom up
        for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++ fq_vert_id)
        {
            vertex_t vert_id = small_queue[fq_vert_id];
            if(scc_id[vert_id] == 0)
            {

                //1. using in edge
                index_t my_beg = bw_beg_pos[vert_id];
                index_t my_end = bw_beg_pos[vert_id+1];

                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    if(vert_id == w)
                        continue;
//                    if(scc_id[w] == 0 && color[w] > color[vert_id])
                    if(scc_id[w] == 0)
                    {
                        if(color[vert_id] < color[w])
                        {
                            color[vert_id] = color[w];
                            if(!color_changed)
                                color_changed = true;
                        }
                    }
                }
                //2. using out edge
                my_beg = fw_beg_pos[vert_id];
                my_end = fw_beg_pos[vert_id+1];

                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = fw_csr[my_beg];
                    if(vert_id == w)
                        continue;
//                    if(scc_id[w] == 0 && color[w] > color[vert_id])
                    if(scc_id[w] == 0)
                    {
                        if(color[vert_id] < color[w])
                        {
                            color[vert_id] = color[w];
                            if(!color_changed)
                                color_changed = true;
                        }
                    }
                }
            }
        }
        #pragma omp barrier
        
        // path compression                
        if(color_changed)
        {
            for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
            {
                vertex_t vert_id = small_queue[fq_vert_id];
                if(scc_id[vert_id] == 0 && color[vert_id] != vert_id)
                {
                    index_t root = color[vert_id];
                    index_t c_depth = 0;
                    while(color[root] != root && c_depth < 100)// && degree_prop[root] <= degree_prop[color[root]])
                    {
                        root = color[root];
                        c_depth ++;
                    }
                    index_t v_id = vert_id;
                    while(v_id != root && color[v_id] != root)// && degree_prop[v_id] < degree_prop[root])
                    {
                        index_t prev = v_id;
                        v_id = color[v_id];
                        color[prev] = root;
                    }
                }
            }
            color_change[tid] = true;
        }
        else
        {
            color_change[tid] = false;
        }
        #pragma omp barrier

        bool final_color_change = false;
        for(index_t i=0; i<thread_count; ++i)
        {
            if(color_change[i])
            {
                final_color_change = true;
                break;
            }
        }
        #pragma omp barrier
        if(final_color_change == false)
        {
            return;
        }
    }

}
#endif
