#ifndef TRIM_2_3_H
#define TRIM_2_3_H

#include "scc_common.h"
#include "util.h"
#include "wtime.h"

//function 1: trim size_2 from scanning the whole graph
inline void trim_2_from_graph(
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
            index_t out_neighbor = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(out_degree == 0)
                    {
                        out_neighbor = w;
                        out_degree = 1;
                    }
                    else
                        if(w != out_neighbor)
                        {
                            out_degree = 2;
                            break;
                        }
                }
            } 
            ///add trim size-1

            //out_degree is 1, check the out_neighbor
            if(out_degree == 1)
            {
                my_beg = fw_beg_pos[out_neighbor];
                my_end = fw_beg_pos[out_neighbor + 1];

                index_t neighbor_out_degree = 0;
                index_t neighbor_out_neighbor = 0;
                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = fw_csr[my_beg];
                    if(scc_id[w] == 0 && w != out_neighbor)
                    {
                        if(neighbor_out_degree == 0)
                        {
                            neighbor_out_degree = 1;
                            neighbor_out_neighbor = w;
                        }
                        else
                            if(w != neighbor_out_neighbor)
                            {
                                neighbor_out_degree = 2;
                                break;
                            }
                    }
                }
                ///add trim size-1: neighbor_out_degree == 0
                if(neighbor_out_degree == 1 && neighbor_out_neighbor == vert_id)
                {
                    scc_id[vert_id] = -2;
                    scc_id[out_neighbor] = -2;
                    continue;
                }
            }
            
            index_t in_degree = 0;
            index_t in_neighbor = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(in_degree == 0)
                    {
                        in_degree = 1; 
                        in_neighbor = w;
                    }
                    else
                        if(w != in_neighbor)
                        {
                            in_degree = 2;
                            break;
                        }
                }
            }
            /// add trim size-1
            //in_degree is 1, check the in_neighbor
            if(in_degree == 1)
            {
                my_beg = bw_beg_pos[in_neighbor];
                my_end = bw_beg_pos[in_neighbor + 1];

                index_t neighbor_in_degree = 0;
                index_t neighbor_in_neighbor = 0;
                
                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    if(scc_id[w] == 0 && w != in_neighbor)
                    {
                        if(neighbor_in_degree == 0)
                        {
                            neighbor_in_degree = 1;
                            neighbor_in_neighbor = w;
                        }
                        else
                            if(w != neighbor_in_neighbor)
                            {
                                neighbor_in_degree = 2;
                                break;
                            }
                    }
                }
                ///add trim size-1: neighbor_out_degree == 0
                if(neighbor_in_degree == 1 && neighbor_in_neighbor == vert_id)
                {
                    scc_id[vert_id] = -2;
                    scc_id[in_neighbor] = -2;
                    continue;
                }
            }
        }
                    
    }
}
//function 2: trim 2 from last frontier queue
inline void trim_2_from_fq(
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
//        if(scc_id[vert_id] == 0)
//        {
            
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id+1];
            index_t out_degree = 0;
            index_t out_neighbor = 0;
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(out_degree == 0)
                    {
                        out_neighbor = w;
                        out_degree = 1;
                    }
                    else
                        if(w != out_neighbor)
                        {
                            out_degree = 2;
                            break;
                        }
                }
            } 
            ///add trim size-1

            //out_degree is 1, check the out_neighbor
            if(out_degree == 1)
            {
                my_beg = fw_beg_pos[out_neighbor];
                my_end = fw_beg_pos[out_neighbor + 1];

                index_t neighbor_out_degree = 0;
                index_t neighbor_out_neighbor = 0;
                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = fw_csr[my_beg];
                    if(scc_id[w] == 0 && w != out_neighbor)
                    {
                        if(neighbor_out_degree == 0)
                        {
                            neighbor_out_degree = 1;
                            neighbor_out_neighbor = w;
                        }
                        else
                            if(w != neighbor_out_neighbor)
                            {
                                neighbor_out_degree = 2;
                                break;
                            }
                    }
                }
                ///add trim size-1: neighbor_out_degree == 0
                if(neighbor_out_degree == 1 && neighbor_out_neighbor == vert_id)
                {
                    scc_id[vert_id] = -2;
                    scc_id[out_neighbor] = -2;
                    continue;
                }
            }
            
            index_t in_degree = 0;
            index_t in_neighbor = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(in_degree == 0)
                    {
                        in_degree = 1; 
                        in_neighbor = w;
                    }
                    else
                        if(w != in_neighbor)
                        {
                            in_degree = 2;
                            break;
                        }
                }
            }
            /// add trim size-1
            //in_degree is 1, check the in_neighbor
            if(in_degree == 1)
            {
                my_beg = bw_beg_pos[in_neighbor];
                my_end = bw_beg_pos[in_neighbor + 1];

                index_t neighbor_in_degree = 0;
                index_t neighbor_in_neighbor = 0;
                
                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    if(scc_id[w] == 0 && w != in_neighbor)
                    {
                        if(neighbor_in_degree == 0)
                        {
                            neighbor_in_degree = 1;
                            neighbor_in_neighbor = w;
                        }
                        else
                            if(w != neighbor_in_neighbor)
                            {
                                neighbor_in_degree = 2;
                                break;
                            }
                    }
                }
                ///add trim size-1: neighbor_out_degree == 0
                if(neighbor_in_degree == 1 && neighbor_in_neighbor == vert_id)
                {
                    scc_id[vert_id] = -2;
                    scc_id[in_neighbor] = -2;
                    continue;
                }
            }
//        }
                    
    }
}

//function 3: trim size_3_1 from scanning the graph
//size-3: type_1, A --> B --> C --> A
inline void trim_3_1_from_graph(
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
            index_t out_neighbor = 0;
            
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                ///self-cycle
                if(scc_id[w] == 0 && vert_id != w)
                {
                    ///redundant edge
                    if(out_degree == 0)
                    {
                        out_degree = 1;
                        out_neighbor = w;
                    }
                    else
                        if(w != out_neighbor)
                        {
                            out_degree = 2;
                            break;
                        }
                }
            }

            //A out_degree is 1, check the out_neighbor
            if(out_degree == 1)
            {
                my_beg = fw_beg_pos[out_neighbor];
                my_end = fw_beg_pos[out_neighbor + 1];

                index_t B_out_degree = 0;
                index_t B_out_neighbor = 0;
                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = fw_csr[my_beg];
                    if(scc_id[w] == 0 && out_neighbor != w && w != vert_id)
                    {
                        if(B_out_degree == 0)
                        {
                            B_out_degree = 1;
                            B_out_neighbor = w;
                        }
                        else
                            if(w != B_out_neighbor)
                            {
                                B_out_degree = 2;
                                break;
                            }
                    }
                }
                ///add trim size-1: neighbor_out_degree == 0
                if(B_out_degree == 1)
                {
                    //B out_degree is 1, check C
                    my_beg = fw_beg_pos[B_out_neighbor];
                    my_end = fw_beg_pos[B_out_neighbor + 1];
                    index_t C_out_degree = 0;
                    index_t C_out_neighbor = 0;

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        if(scc_id[w] == 0 && w != B_out_neighbor && w != out_neighbor)
                        {
                            if(C_out_degree == 0)
                            {
                                C_out_degree = 1;
                                C_out_neighbor = w;
                            }
                            else
                                if(w != C_out_neighbor)
                                {
                                    C_out_degree = 2;
                                    break;
                                }
                        }
                    }
                    if(C_out_degree == 1 && C_out_neighbor == vert_id)
                    {
                        scc_id[vert_id] = -3;
                        scc_id[out_neighbor] = -3;
                        scc_id[B_out_neighbor] = -3;
//                            printf("Out_A = %d, B = %d, C = %d\n", vert_id, out_neighbor, B_out_neighbor);
//                        size_3_type_1++;
                    }
                }
            }
            
            index_t in_degree = 0;
            index_t in_neighbor = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(in_degree == 0)
                    {
                        in_neighbor = w;
                        in_degree = 1;
                    }
                    else
                        if(w != in_neighbor)
                        {
                            in_degree = 2;
                            break;
                        }
                }
            }
            /// add trim size-1
            //in_degree is 1, check the in_neighbor
            if(in_degree == 1 && in_neighbor != vert_id)
            {
                my_beg = bw_beg_pos[in_neighbor];
                my_end = bw_beg_pos[in_neighbor + 1];

                index_t B_in_degree = 0;
                index_t B_in_neighbor = 0;
                
                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    if(scc_id[w] == 0 && in_neighbor != w && w != vert_id)
                    {
                        if(B_in_degree == 0)
                        {
                            B_in_degree = 1;
                            B_in_neighbor = w;
                        }
                        else
                            if(w != B_in_neighbor)
                            {
                                B_in_degree = 2;
                                break;
                            }
                    }
                }
                ///add trim size-1: neighbor_in_degree == 0
                if(B_in_degree == 1)
                {
                    //B in_degree is 1, check C
                    my_beg = bw_beg_pos[B_in_neighbor];
                    my_end = bw_beg_pos[B_in_neighbor + 1];
                    index_t C_in_degree = 0;
                    index_t C_in_neighbor = 0;

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = bw_csr[my_beg];
                        if(scc_id[w] == 0 && w != B_in_neighbor && w != in_neighbor)
                        {
                            if(C_in_degree == 0)
                            {
                                C_in_degree = 1;
                                C_in_neighbor = w;
                            }
                            else
                                if(w != C_in_neighbor)
                                {
                                    C_in_degree = 2;
                                    break;
                                }
                        }
                    }
                    if(C_in_degree == 1 && C_in_neighbor == vert_id)
                    {
                        scc_id[vert_id] = -3;
                        scc_id[in_neighbor] = -3;
                        scc_id[B_in_neighbor] = -3;
//                        size_3_type_1++;
                    }
                }
            }
        }
                    
    }
}

inline void trim_3_1_from_fq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *small_queue
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
            index_t out_neighbor = 0;
            
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                ///self-cycle
                if(scc_id[w] == 0 && vert_id != w)
                {
                    ///redundant edge
                    if(out_degree == 0)
                    {
                        out_degree = 1;
                        out_neighbor = w;
                    }
                    else
                        if(w != out_neighbor)
                        {
                            out_degree = 2;
                            break;
                        }
                }
            }

            //A out_degree is 1, check the out_neighbor
            if(out_degree == 1)
            {
                my_beg = fw_beg_pos[out_neighbor];
                my_end = fw_beg_pos[out_neighbor + 1];

                index_t B_out_degree = 0;
                index_t B_out_neighbor = 0;
                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = fw_csr[my_beg];
                    if(scc_id[w] == 0 && out_neighbor != w && w != vert_id)
                    {
                        if(B_out_degree == 0)
                        {
                            B_out_degree = 1;
                            B_out_neighbor = w;
                        }
                        else
                            if(w != B_out_neighbor)
                            {
                                B_out_degree = 2;
                                break;
                            }
                    }
                }
                ///add trim size-1: neighbor_out_degree == 0
                if(B_out_degree == 1)
                {
                    //B out_degree is 1, check C
                    my_beg = fw_beg_pos[B_out_neighbor];
                    my_end = fw_beg_pos[B_out_neighbor + 1];
                    index_t C_out_degree = 0;
                    index_t C_out_neighbor = 0;

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        if(scc_id[w] == 0 && w != B_out_neighbor && w != out_neighbor)
                        {
                            if(C_out_degree == 0)
                            {
                                C_out_degree = 1;
                                C_out_neighbor = w;
                            }
                            else
                                if(w != C_out_neighbor)
                                {
                                    C_out_degree = 2;
                                    break;
                                }
                        }
                    }
                    if(C_out_degree == 1 && C_out_neighbor == vert_id)
                    {
                        scc_id[vert_id] = -3;
                        scc_id[out_neighbor] = -3;
                        scc_id[B_out_neighbor] = -3;
//                            printf("Out_A = %d, B = %d, C = %d\n", vert_id, out_neighbor, B_out_neighbor);
//                        size_3_type_1++;
                    }
                }
            }
            
            index_t in_degree = 0;
            index_t in_neighbor = 0;
            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id+1];
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(in_degree == 0)
                    {
                        in_neighbor = w;
                        in_degree = 1;
                    }
                    else
                        if(w != in_neighbor)
                        {
                            in_degree = 2;
                            break;
                        }
                }
            }
            /// add trim size-1
            //in_degree is 1, check the in_neighbor
            if(in_degree == 1 && in_neighbor != vert_id)
            {
                my_beg = bw_beg_pos[in_neighbor];
                my_end = bw_beg_pos[in_neighbor + 1];

                index_t B_in_degree = 0;
                index_t B_in_neighbor = 0;
                
                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    if(scc_id[w] == 0 && in_neighbor != w && w != vert_id)
                    {
                        if(B_in_degree == 0)
                        {
                            B_in_degree = 1;
                            B_in_neighbor = w;
                        }
                        else
                            if(w != B_in_neighbor)
                            {
                                B_in_degree = 2;
                                break;
                            }
                    }
                }
                ///add trim size-1: neighbor_in_degree == 0
                if(B_in_degree == 1)
                {
                    //B in_degree is 1, check C
                    my_beg = bw_beg_pos[B_in_neighbor];
                    my_end = bw_beg_pos[B_in_neighbor + 1];
                    index_t C_in_degree = 0;
                    index_t C_in_neighbor = 0;

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = bw_csr[my_beg];
                        if(scc_id[w] == 0 && w != B_in_neighbor && w != in_neighbor)
                        {
                            if(C_in_degree == 0)
                            {
                                C_in_degree = 1;
                                C_in_neighbor = w;
                            }
                            else
                                if(w != C_in_neighbor)
                                {
                                    C_in_degree = 2;
                                    break;
                                }
                        }
                    }
                    if(C_in_degree == 1 && C_in_neighbor == vert_id)
                    {
                        scc_id[vert_id] = -3;
                        scc_id[in_neighbor] = -3;
                        scc_id[B_in_neighbor] = -3;
//                        size_3_type_1++;
                    }
                }
            }
        }
                    
    }
}

//function 4: trim size_3_2 from scanning the graph
//size-3: type_1, A --> B --> A --> C --> A 
//starting from hub vertex A
inline void trim_3_2_from_graph(
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
            index_t my_end = fw_beg_pos[vert_id + 1];
            index_t out_degree = 0;
            index_t out_neighbor_B = -1;
            index_t out_neighbor_C = -1;
            
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(out_degree == 0)
                    {
                        out_neighbor_B = w;
                        out_degree = 1;
                    }
                    else
                        if(out_degree == 1 && out_neighbor_B != w)
                        {
                            out_degree = 2;
                            out_neighbor_C = w;
                        }
                        else
                            if(out_degree == 2 && out_neighbor_B != w && out_neighbor_C != w)
                            {
                                out_degree = 3;
                                break;
                            }
                }
            }
            if(out_degree == 2)
            {
                //check whether B/C can point to A

                my_beg = fw_beg_pos[out_neighbor_B];
                my_end = fw_beg_pos[out_neighbor_B + 1];
                index_t B_out_degree = 0;
                index_t B_out_neighbor = 0;

                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = fw_csr[my_beg];
                    if(scc_id[w] == 0 && w != out_neighbor_B && w != out_neighbor_C)
                    {
                        if(B_out_degree == 0)
                        {
                            B_out_degree = 1;
                            B_out_neighbor = w;
                        }
                        else
                        {
                            B_out_degree = 2;
                            break;
                        }
                    }
                }
                if(B_out_degree == 1 && B_out_neighbor == vert_id)
                {
                    my_beg = fw_beg_pos[out_neighbor_C];
                    my_end = fw_beg_pos[out_neighbor_C + 1];
                    index_t C_out_degree = 0;
                    index_t C_out_neighbor = 0;

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        if(scc_id[w] == 0 && w != out_neighbor_C && w != out_neighbor_B)
                        {
                            if(C_out_degree == 0)
                            {
                                C_out_degree = 1;
                                C_out_neighbor = w;
                            }
                            else
                            {
                                C_out_degree = 2;
                                break;
                            }
                        }
                    }
                    if(C_out_degree == 1 && C_out_neighbor == vert_id)
                    {
                        scc_id[vert_id] = -4;
                        scc_id[out_neighbor_B] = -4;
                        scc_id[out_neighbor_C] = -4;
//                        size_3_type_2 ++;
                    }
                }
            }

            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id + 1];
            index_t in_degree = 0;
            index_t in_neighbor_B = -1;
            index_t in_neighbor_C = -1;
            
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(in_degree == 0)
                    {
                        in_neighbor_B = w;
                        in_degree = 1;
                    }
                    else
                        if(in_degree == 1 && in_neighbor_B != w)
                        {
                            in_degree = 2;
                            in_neighbor_C = w;
                        }
                        else
                            if(in_degree == 2 && in_neighbor_B != w && in_neighbor_C != w)
                            {
                                in_degree = 3;
                                break;
                            }
                }
            }
            if(in_degree == 2)
            {
                //check whether B/C can point to A

                my_beg = bw_beg_pos[in_neighbor_B];
                my_end = bw_beg_pos[in_neighbor_B + 1];
                index_t B_in_degree = 0;
                index_t B_in_neighbor = 0;

                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    if(scc_id[w] == 0 && w != in_neighbor_B && w != in_neighbor_C)
                    {
                        if(B_in_degree == 0)
                        {
                            B_in_degree = 1;
                            B_in_neighbor = w;
                        }
                        else
                        {
                            B_in_degree = 2;
                            break;
                        }
                    }
                }
                if(B_in_degree == 1 && B_in_neighbor == vert_id)
                {
                    my_beg = bw_beg_pos[in_neighbor_C];
                    my_end = bw_beg_pos[in_neighbor_C + 1];
                    index_t C_in_degree = 0;
                    index_t C_in_neighbor = 0;

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = bw_csr[my_beg];
                        if(scc_id[w] == 0 && w != in_neighbor_C && w != in_neighbor_B)
                        {
                            if(C_in_degree == 0)
                            {
                                C_in_degree = 1;
                                C_in_neighbor = w;
                            }
                            else
                            {
                                C_in_degree = 2;
                                break;
                            }
                        }
                    }
                    if(C_in_degree == 1 && C_in_neighbor == vert_id)
                    {
                        scc_id[vert_id] = -4;
                        scc_id[in_neighbor_B] = -4;
                        scc_id[in_neighbor_C] = -4;
//                        size_3_type_2 ++;
                    }
                }
            }


        }
    }
}

inline void trim_3_2_from_fq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *small_queue
        )
{
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            index_t my_beg = fw_beg_pos[vert_id];
            index_t my_end = fw_beg_pos[vert_id + 1];
            index_t out_degree = 0;
            index_t out_neighbor_B = -1;
            index_t out_neighbor_C = -1;
            
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = fw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(out_degree == 0)
                    {
                        out_neighbor_B = w;
                        out_degree = 1;
                    }
                    else
                        if(out_degree == 1 && out_neighbor_B != w)
                        {
                            out_degree = 2;
                            out_neighbor_C = w;
                        }
                        else
                            if(out_degree == 2 && out_neighbor_B != w && out_neighbor_C != w)
                            {
                                out_degree = 3;
                                break;
                            }
                }
            }
            if(out_degree == 2)
            {
                //check whether B/C can point to A

                my_beg = fw_beg_pos[out_neighbor_B];
                my_end = fw_beg_pos[out_neighbor_B + 1];
                index_t B_out_degree = 0;
                index_t B_out_neighbor = 0;

                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = fw_csr[my_beg];
                    if(scc_id[w] == 0 && w != out_neighbor_B && w != out_neighbor_C)
                    {
                        if(B_out_degree == 0)
                        {
                            B_out_degree = 1;
                            B_out_neighbor = w;
                        }
                        else
                        {
                            B_out_degree = 2;
                            break;
                        }
                    }
                }
                if(B_out_degree == 1 && B_out_neighbor == vert_id)
                {
                    my_beg = fw_beg_pos[out_neighbor_C];
                    my_end = fw_beg_pos[out_neighbor_C + 1];
                    index_t C_out_degree = 0;
                    index_t C_out_neighbor = 0;

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        if(scc_id[w] == 0 && w != out_neighbor_C && w != out_neighbor_B)
                        {
                            if(C_out_degree == 0)
                            {
                                C_out_degree = 1;
                                C_out_neighbor = w;
                            }
                            else
                            {
                                C_out_degree = 2;
                                break;
                            }
                        }
                    }
                    if(C_out_degree == 1 && C_out_neighbor == vert_id)
                    {
                        scc_id[vert_id] = -4;
                        scc_id[out_neighbor_B] = -4;
                        scc_id[out_neighbor_C] = -4;
//                        size_3_type_2 ++;
                    }
                }
            }

            my_beg = bw_beg_pos[vert_id];
            my_end = bw_beg_pos[vert_id + 1];
            index_t in_degree = 0;
            index_t in_neighbor_B = -1;
            index_t in_neighbor_C = -1;
            
            for(; my_beg < my_end; ++my_beg)
            {
                index_t w = bw_csr[my_beg];
                if(scc_id[w] == 0 && w != vert_id)
                {
                    if(in_degree == 0)
                    {
                        in_neighbor_B = w;
                        in_degree = 1;
                    }
                    else
                        if(in_degree == 1 && in_neighbor_B != w)
                        {
                            in_degree = 2;
                            in_neighbor_C = w;
                        }
                        else
                            if(in_degree == 2 && in_neighbor_B != w && in_neighbor_C != w)
                            {
                                in_degree = 3;
                                break;
                            }
                }
            }
            if(in_degree == 2)
            {
                //check whether B/C can point to A

                my_beg = bw_beg_pos[in_neighbor_B];
                my_end = bw_beg_pos[in_neighbor_B + 1];
                index_t B_in_degree = 0;
                index_t B_in_neighbor = 0;

                for(; my_beg < my_end; ++my_beg)
                {
                    index_t w = bw_csr[my_beg];
                    if(scc_id[w] == 0 && w != in_neighbor_B && w != in_neighbor_C)
                    {
                        if(B_in_degree == 0)
                        {
                            B_in_degree = 1;
                            B_in_neighbor = w;
                        }
                        else
                        {
                            B_in_degree = 2;
                            break;
                        }
                    }
                }
                if(B_in_degree == 1 && B_in_neighbor == vert_id)
                {
                    my_beg = bw_beg_pos[in_neighbor_C];
                    my_end = bw_beg_pos[in_neighbor_C + 1];
                    index_t C_in_degree = 0;
                    index_t C_in_neighbor = 0;

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = bw_csr[my_beg];
                        if(scc_id[w] == 0 && w != in_neighbor_C && w != in_neighbor_B)
                        {
                            if(C_in_degree == 0)
                            {
                                C_in_degree = 1;
                                C_in_neighbor = w;
                            }
                            else
                            {
                                C_in_degree = 2;
                                break;
                            }
                        }
                    }
                    if(C_in_degree == 1 && C_in_neighbor == vert_id)
                    {
                        scc_id[vert_id] = -4;
                        scc_id[in_neighbor_B] = -4;
                        scc_id[in_neighbor_C] = -4;
//                        size_3_type_2 ++;
                    }
                }
            }


        }
    }
}

#endif
