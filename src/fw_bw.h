#ifndef FW_BW_H
#define FW_BW_H

#include "wtime.h"
#include "util.h"
#include "graph.h"
#include <set>
#include <iostream>

inline void fw_bfs(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *fw_sa,
        index_t *front_comm,
        index_t *work_comm,
        vertex_t root,
        index_t tid,
        index_t thread_count,
        double alpha,
        double beta,
        vertex_t edge_count,
        vertex_t vert_count
        )
{
    depth_t level = 0;
    fw_sa[root] = 0;
    bool is_top_down = true;	
    bool is_top_down_queue = false;
    index_t queue_size = vert_count / thread_count;
    while(true)
    {
        double ltm= wtime();
        index_t front_count=0;
        index_t my_work_next=0;
        index_t my_work_curr=0;

        if(is_top_down)
        {
            for(vertex_t vert_id=vert_beg; vert_id<vert_end; vert_id++)
            {
                if(scc_id[vert_id] == 0 && fw_sa[vert_id]==level)
                {
                    index_t my_beg = fw_beg_pos[vert_id];
                    index_t my_end = fw_beg_pos[vert_id+1];

                    for(; my_beg<my_end; my_beg++)
                    {
                        vertex_t nebr=fw_csr[my_beg];
                        if(scc_id[nebr] == 0 && fw_sa[nebr] == -1)
                        {
                            fw_sa[nebr] = level+1;
                            my_work_next+=fw_beg_pos[nebr+1]-fw_beg_pos[nebr];

                            front_count++;
                        }
                    }
                }
            }
            work_comm[tid]=my_work_next;
        }
        else
            if(!is_top_down_queue)
            {
                for(vertex_t vert_id=vert_beg; vert_id<vert_end; vert_id++)
                {
                    if(scc_id[vert_id] == 0 && fw_sa[vert_id]== -1)
                    {
                        index_t my_beg = bw_beg_pos[vert_id];
                        index_t my_end = bw_beg_pos[vert_id+1];
                        my_work_curr+=my_end-my_beg;

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=bw_csr[my_beg];
                            if(scc_id[vert_id] == 0 && fw_sa[nebr] != -1)
                            {
                                fw_sa[vert_id] = level+1;
                                front_count++;
                                break;
                            }
                        }
                    }
                }
                work_comm[tid]=my_work_curr;
                
            }
            else
            {
                index_t *q = new index_t[queue_size];
                index_t head = 0;
                index_t tail = 0;
                //std::queue<index_t> q;

                //Option 1: put current level vertices into fq
                for(vertex_t vert_id=vert_beg; vert_id<vert_end; vert_id++)
                {
                    if(scc_id[vert_id] == 0 && fw_sa[vert_id] == level)
                    {
                        q[tail++] = vert_id;
                    }

                }
                while(head != tail)
                {
                    vertex_t temp_v = q[head++];
                    if(head == queue_size)
                        head = 0;
                    index_t my_beg = fw_beg_pos[temp_v];
                    index_t my_end = fw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        
                        if(scc_id[w] == 0 && fw_sa[w] == -1)
                        {
                            q[tail++] = w;
                            if(tail == queue_size)
                                tail = 0;
                            fw_sa[w] = level + 1;
                        }
                    }
                }
                delete[] q;
                break;
            }

        front_comm[tid]=front_count;

        #pragma omp barrier
        front_count=0;
        my_work_next=0;

        for(index_t i=0;i<thread_count;++i)
        {
            front_count += front_comm[i];
            my_work_next += work_comm[i];
        }
            
        if(front_count == 0) break;
        
        if(is_top_down && my_work_next>(alpha*edge_count)) 
        {
            is_top_down=false;
    //		if(tid==0)
    //			std::cout<<"--->Switch to bottom up"<<my_work_next
    //				<<" "<<edge_count<<"<----\n";
    //		if(tid==0)
    //		{
    //			long todo=0;
    //			for(int i=0;i<vert_count;i++)
    //				if(fw_sa[i]==-1)
    //					todo+=fw_beg_pos[i+1]-fw_beg_pos[i];

    //			std::cout<<"Edges connected to unvisited: "<<todo<<"\n";
    //		}

        }	
        if(!is_top_down && my_work_next < (beta * edge_count))
        {
            is_top_down = true;
//				if(tid==0)
//					std::cout<<"--->Switch to top down"<<my_work_next
//						<<" "<<edge_count<<"<----\n";
//				if(tid==0)
//				{
//					long todo=0;
//					for(int i=0;i<vert_count;i++)
//						if(fw_sa[i]==-1)
//							todo+=fw_beg_pos[i+1]-fw_beg_pos[i];
//
//					std::cout<<"Edges connected to unvisited: "<<todo<<"\n";
//				}
        }
        
        #pragma omp barrier
    
//			if(tid==0) std::cout<<"Level-"<<(int)level
//				<<"-frontier-time-futurework(currwork in btup): "
//				<<front_count<<" "
//				<<(wtime() - ltm) * 1000<<" ms "
//				<<my_work_next<<"\n";
//			#pragma omp barrier

        level ++;
    }
}

inline void bw_bfs(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *fw_sa,
        vertex_t *bw_sa,
        index_t *front_comm,
        index_t *work_comm,
        vertex_t root,
        index_t tid,
        index_t thread_count,
        double alpha,
        double beta,
        vertex_t edge_count,
        vertex_t vert_count
        )
{
    bw_sa[root] = 0;
    bool is_top_down = true;
    bool is_top_down_queue = false;
    index_t level = 0;
    scc_id[root] = 1;
    index_t queue_size = vert_count / thread_count;
    while(true)
    {
        double ltm= wtime();
        index_t front_count=0;
        index_t my_work_next=0;
        index_t my_work_curr=0;

        if(is_top_down)
        {
            for(vertex_t vert_id=vert_beg; vert_id<vert_end; vert_id++)
            {
                if(scc_id[vert_id] == 1 && bw_sa[vert_id]==level)
                {
                    index_t my_beg = bw_beg_pos[vert_id];
                    index_t my_end = bw_beg_pos[vert_id+1];

//                        printf("my_beg = %d, my_end = %d\n", my_beg, my_end);
                    for(; my_beg<my_end; my_beg++)
                    {
                        vertex_t nebr=bw_csr[my_beg];
                        if(scc_id[nebr] == 0 && bw_sa[nebr] == -1 && fw_sa[nebr] != -1)
                        {
                            bw_sa[nebr] = level+1;
                            my_work_next+=bw_beg_pos[nebr+1]-bw_beg_pos[nebr];
                            front_count++;
                            scc_id[nebr] = 1;	
                        }
                    }
                }
            }
            work_comm[tid]=my_work_next;
        }
        else
            if(!is_top_down_queue)
            {
                for(vertex_t vert_id=vert_beg; vert_id<vert_end; vert_id++)
                {
                    if(scc_id[vert_id] == 0 && bw_sa[vert_id] == -1 && fw_sa[vert_id] != -1)
                    {
                        index_t my_beg = fw_beg_pos[vert_id];
                        index_t my_end = fw_beg_pos[vert_id+1];
                        my_work_curr+=my_end-my_beg;

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=fw_csr[my_beg];
                            if(scc_id[nebr] == 1 && bw_sa[nebr] != -1)//fw_sa[nebr] != -1)
                            {
                                bw_sa[vert_id] = level+1;
                                front_count++;
                                scc_id[vert_id] = 1;
                                break;
                            }
                        }
                    }
                }
                work_comm[tid]=my_work_curr;
            }
            else
            {
//                std::queue<index_t> q;
                index_t *q = new index_t[queue_size];
                index_t head = 0;
                index_t tail = 0;
                
				for(vertex_t vert_id=vert_beg; vert_id<vert_end; vert_id++)
                {
                    if(scc_id[vert_id] == 1 && bw_sa[vert_id] == level)
                    {
//                        q.push(vert_id);
                        q[tail++] = vert_id;
                        //impossible here
//                        if(tail == queue_size)
//                            tail = 0;
                    }
                }

//                printf("queue_size = %d, tail = %d\n", queue_size, tail);
//                printf("q[%d].size = %d\n", tid, q.size());
                while(head != tail)
                {
//                    index_t temp_v = q.front();
//                    q.pop();
                    vertex_t temp_v = q[head++];
                    if(head == queue_size)
                        head = 0;
                    index_t my_beg = bw_beg_pos[temp_v];
                    index_t my_end = bw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = bw_csr[my_beg];
                        
                        if(scc_id[w] == 0 && bw_sa[w] == -1 && fw_sa[w] != -1)
                        {
//                            q.push(w);
                            q[tail++] = w;
                            if(tail == queue_size)
                                tail = 0;
                            scc_id[w] = 1;
                            bw_sa[w] = level + 1;
                        }
                    }
//                    if(tid == 0)
//                        printf("head, %d; tail, %d\n", head, tail);
                }
//                break;
//                printf("head = %d\n", head);
                delete[] q;
            }

        front_comm[tid]=front_count;

        #pragma omp barrier
        front_count=0;
        my_work_next=0;

        for(index_t i=0;i<thread_count;++i)
        {
            front_count += front_comm[i];
            my_work_next += work_comm[i];
        }
            
        if(front_count == 0) 
        {
            break;
        }
        
        if(is_top_down && my_work_next>(alpha*edge_count)) 
        {
            is_top_down=false;
//				if(tid==0)
//					std::cout<<"--->Switch to bottom up"<<my_work_next
//						<<" "<<edge_count<<"<----\n";
//				if(tid==0)
//				{
//					long todo=0;
//					for(int i=0;i<vert_count;i++)
//						if(bw_sa[i]==-1)
//							todo+=bw_beg_pos[i+1]-bw_beg_pos[i];
//
//					std::cout<<"Edges connected to unvisited: "<<todo<<"\n";
//				}

        }	
        if(!is_top_down && my_work_next < (beta * edge_count))
        {
            is_top_down = true;
//				if(tid==0)
//					std::cout<<"--->Switch to top down"<<my_work_next
//						<<" "<<edge_count<<"<----\n";
//				if(tid==0)
//				{
//					long todo=0;
//					for(int i=0;i<vert_count;i++)
//						if(bw_sa[i]==-1)
//							todo+=bw_beg_pos[i+1]-bw_beg_pos[i];
//
//					std::cout<<"Edges connected to unvisited: "<<todo<<"\n";
//				}
        }
        
        #pragma omp barrier
//		
//			if(tid==0) std::cout<<"Level-"<<(int)level
//				<<"-frontier-time-futurework(currwork in btup): "
//				<<front_count<<" "
//				<<(wtime() - ltm) * 1000<<" ms "
//				<<my_work_next<<"\n";
//			#pragma omp barrier

        level ++;
    }
}

inline void fw_bfs_fq_queue(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *fw_sa,
        index_t *vertex_cur,
        index_t *vertex_front,
        vertex_t root,
        index_t tid,
        index_t thread_count,
        const int alpha,
        const int beta,
        const int gamma,
        vertex_t *frontier_queue,
        vertex_t fq_size,
        const double avg_degree,
        vertex_t vertex_visited,
        vertex_t *temp_queue,
        index_t *prefix_sum,
        vertex_t upper_bound,
        vertex_t *thread_queue
        )
{
    depth_t level = 0;
    fw_sa[root] = 0;
    temp_queue[0] = root;
    vertex_t queue_size = 1;
//    vertex_t upper_bound = vert_count / thread_count;
//    vertex_t *thread_queue = new vertex_t[upper_bound];
//    vertex_t root_out_degree = fw_beg_pos[root+1] - fw_beg_pos[root];
    bool is_top_down = true;
//    bool is_top_down_async = false;
//    if(VERBOSE)
//    {
//        if(tid == 0)
//        {
//            printf("out_degree, %d, limit, %.3lf\n", root_out_degree, alpha * beta * fq_size);
//        }
//    }
//    if(root_out_degree < alpha * beta * fq_size)
//    {
//        is_top_down_async = true;
//    }
    bool is_top_down_queue = false;
//    index_t queue_size = fq_size / thread_count;
    #pragma omp barrier
    while(true)
    {
        double ltm= wtime();
        vertex_t vertex_frontier = 0;
        
//        printf("tid, %d, %d, %d, %d\n", tid, step, queue_beg, queue_end);
        #pragma omp barrier
        if(is_top_down)
        {
            vertex_t step = queue_size / thread_count;
            vertex_t queue_beg = tid * step;
            vertex_t queue_end = (tid == thread_count - 1 ? queue_size: queue_beg + step);
            for(vertex_t q_vert_id=queue_beg; q_vert_id<queue_end; q_vert_id++)
            {
                vertex_t vert_id = temp_queue[q_vert_id];
                //in fq, scc_id[vert_id] is always not 0
                if(scc_id[vert_id] == 0 && fw_sa[vert_id] != -1)
                {
                    index_t my_beg = fw_beg_pos[vert_id];
                    index_t my_end = fw_beg_pos[vert_id+1];
                    for(; my_beg<my_end; my_beg++)
                    {
                        vertex_t nebr=fw_csr[my_beg];
                        if(scc_id[nebr] == 0 && fw_sa[nebr] == -1)
                        {
                            fw_sa[nebr] = level+1;
                            thread_queue[vertex_frontier] = nebr;
                            vertex_frontier++;
                        }
                    }
                }
            }
        }
        else
            if(!is_top_down_queue)
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 0 && fw_sa[vert_id] == -1)
                    {
                        index_t my_beg = bw_beg_pos[vert_id];
                        index_t my_end = bw_beg_pos[vert_id+1];
//                        my_work_curr+=my_end-my_beg;

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=bw_csr[my_beg];
                            if(scc_id[nebr] == 0 && fw_sa[nebr] != -1)
//                            if(scc_id[vert_id] == 0 && fw_sa[nebr] == level)
                            {
                                fw_sa[vert_id] = level+1;
                                vertex_frontier++;
//                                front_count++;
                                break;
                            }
                        }
                    }
                }
            }
            else
            {
                vertex_t end_queue = upper_bound;
                index_t head = 0;
                index_t tail = 0;
                //std::queue<index_t> q;
                vertex_t step = queue_size / thread_count;
                vertex_t queue_beg = tid * step;
                vertex_t queue_end = (tid == thread_count - 1 ? queue_size: queue_beg + step);

                //Option 1: put current level vertices into fq
                for(vertex_t q_vert_id=queue_beg; q_vert_id<queue_end; q_vert_id++)
                {
                    thread_queue[tail] = temp_queue[q_vert_id];
                    tail ++;
                }
                while(head != tail)
                {
                    vertex_t temp_v = thread_queue[head++];
//                    front_count ++;
                    if(head == end_queue)
                        head = 0;
                    index_t my_beg = fw_beg_pos[temp_v];
                    index_t my_end = fw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        
                        if(scc_id[w] == 0 && fw_sa[w] == -1)
                        {
                            thread_queue[tail++] = w;
                            if(tail == end_queue)
                                tail = 0;
                            fw_sa[w] = level + 1;
                        }
                    }
                }
            }

        
        vertex_front[tid] = vertex_frontier;
        #pragma omp barrier
        vertex_frontier = 0;

        for(index_t i=0; i<thread_count; ++i)
        {
            vertex_frontier += vertex_front[i];
        }
        vertex_visited += vertex_frontier;
//        #pragma omp barrier
        if(VERBOSE)
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
			if(tid==0 && level < 50) 
                std::cout<<"Level-"<<(int)level<<" "
//				<<"-frontier-time-visited:"
				<<vertex_frontier<<" "
                <<fq_size<<" "
                <<(double)(fq_size)/vertex_frontier<<" "
				<<(wtime() - ltm) * 1000<<"ms "
				<<vertex_visited<<" "
                <<edge_frontier<<" "
                <<edge_remaider<<" "
                <<edge_remaider/edge_frontier<<"\n";
        }
        
        if(vertex_frontier == 0) break;
        
        if(is_top_down) 
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;
//            printf("edge_remainder/alpha = %g, edge_froniter = %g\n", edge_remainder / alpha, edge_frontier);
            if(!is_top_down_queue && (edge_remainder / alpha) < edge_frontier)
            {
                is_top_down = false;
                if(VERBOSE)
                {
                    if(tid==0)
                    {
//                        double Nf = vertex_frontier;
//                        double Nu = fq_size - vertex_visited;
//                        double Mf = Nf * Nf / Nu + avg_degree * (Nu - Nf);
//                        printf("mf=%.0lf, mu=%.0lf, alpha=%d, Mf=%.0lf\n", edge_frontier, edge_remainder, ALPHA, Mf);
                        std::cout<<"--->Switch to bottom up\n";
                    }
                }
            }

        }
        else
            if((!is_top_down && !is_top_down_queue && (fq_size*1.0/beta) > vertex_frontier) || (!is_top_down && !is_top_down_queue && level > gamma))
//            if(level > 10)
            {
                //if(!is_top_down_queue)
                
                vertex_frontier = 0;
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 0 && fw_sa[vert_id] == level + 1)
                    {
                        thread_queue[vertex_frontier] = vert_id;
                        vertex_frontier++;
                    }
                }
                vertex_front[tid] = vertex_frontier;
                
                is_top_down = false;
                is_top_down_queue = true;

                if(VERBOSE)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to top down queue\n";
                }
            }
//                is_top_down = true;
//                is_top_down_queue = true;
//                if(VERBOSE && tid==0)
//                    std::cout<<"--->Switch back to top down\n";
//                vertex_frontier = 0;
//                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
//                {
//                    vertex_t vert_id = frontier_queue[fq_vert_id];
//                    if(scc_id[vert_id] == 0 && fw_sa[vert_id] == level + 1)
//                    {
//                        thread_queue[vertex_frontier] = vert_id;
//                        vertex_frontier++;
//                    }
//                }
//                vertex_front[tid] = vertex_frontier;
//            }
        
//            if(is_top_down && !is_top_down_queue && (fq_size*1.0/BETA) > vertex_frontier)
// switch to async top down queue

//        if(is_top_down && level > gamma)
//        {
//            // get the frontier queue from bottom up
//            if(!is_top_down_queue)
//            {
//                vertex_frontier = 0;
//                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
//                {
//                    vertex_t vert_id = frontier_queue[fq_vert_id];
//                    if(scc_id[vert_id] == 0 && fw_sa[vert_id] == level + 1)
//                    {
//                        thread_queue[vertex_frontier] = vert_id;
//                        vertex_frontier++;
//                    }
//                }
//                vertex_front[tid] = vertex_frontier;
//            }
//            
//            is_top_down = false;
//            is_top_down_queue = true;
//
//            if(VERBOSE)
//            {
//                if(tid==0)
//                    std::cout<<"--->Switch to top down queue\n";
//            }
//        }
        
        #pragma omp barrier

        if(is_top_down || is_top_down_queue)
        {
            get_queue(thread_queue,
                    vertex_front,
                    prefix_sum,
                    tid,
                    temp_queue);
            queue_size = prefix_sum[thread_count-1] + vertex_front[thread_count-1];
//            if(VERBOSE)
//            {
//                if(tid == 0)
//                    printf("queue_size, %d\n", queue_size);
//            }
        }
        #pragma omp barrier
        
        level ++;
    }
//    if(tid == 0)
//    {
//        printf("fw_level, %d\n", level);
//    }
//    delete[] thread_queue;
}

inline void fw_bfs_fq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *fw_sa,
        index_t *vertex_cur,
        index_t *vertex_front,
        vertex_t root,
        index_t tid,
        index_t thread_count,
        double alpha,
        double beta,
        vertex_t *frontier_queue,
        vertex_t fq_size,
        const double avg_degree,
        vertex_t vertex_visited
        )
{
    depth_t level = 0;
    fw_sa[root] = 0;
    vertex_t root_out_degree = fw_beg_pos[root+1] - fw_beg_pos[root];
    bool is_top_down = true;
    bool is_top_down_async = false;
    if(VERBOSE)
    {
        if(tid == 0)
        {
            printf("out_degree, %d, limit, %.3lf\n", root_out_degree, alpha * beta * fq_size);
        }
    }
    if(root_out_degree < alpha * beta * fq_size)
    {
        is_top_down_async = true;
    }
    bool is_top_down_queue = false;
    index_t queue_size = fq_size / thread_count;
    #pragma omp barrier
    
    while(true)
    {
//        #pragma omp barrier
        double ltm= wtime();
//        index_t front_count=0;
//        index_t my_work_next=0;
//        index_t my_work_curr=0;
        
        vertex_t vertex_frontier = 0;
//        vertex_t vertex_current = 0;

        if(is_top_down)
        {
            if(is_top_down_async)
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    //in fq, scc_id[vert_id] is always not 0
                    if(scc_id[vert_id] == 0 && (fw_sa[vert_id]==level || fw_sa[vert_id]==level+1))
                    {
                        index_t my_beg = fw_beg_pos[vert_id];
                        index_t my_end = fw_beg_pos[vert_id+1];

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=fw_csr[my_beg];
                            if(scc_id[nebr] == 0 && fw_sa[nebr] == -1)
                            {
                                fw_sa[nebr] = level+1;
//                                edge_frontier += fw_beg_pos[nebr+1] - fw_beg_pos[nebr];
                                vertex_frontier ++;
//                                my_work_next+=fw_beg_pos[nebr+1]-fw_beg_pos[nebr];
//
//                                front_count++;
                            }
                        }
                    }
                }
            }
            else
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    //in fq, scc_id[vert_id] is always not 0
                    if(scc_id[vert_id] == 0 && fw_sa[vert_id]==level)
    //                if(scc_id[vert_id] == 0 && (fw_sa[vert_id]==level || fw_sa[vert_id]==level+1))
                    {
                        index_t my_beg = fw_beg_pos[vert_id];
                        index_t my_end = fw_beg_pos[vert_id+1];

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=fw_csr[my_beg];
                            if(scc_id[nebr] == 0 && fw_sa[nebr] == -1)
                            {
                                fw_sa[nebr] = level+1;
//                                my_work_next+=fw_beg_pos[nebr+1]-fw_beg_pos[nebr];
//                                front_count++;
                                vertex_frontier++;
                            }
                        }
                    }
                }
            }
//            vertex_front[tid] = vertex_frontier;
        }
        else
            if(!is_top_down_queue)
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 0 && fw_sa[vert_id]== -1)
                    {
                        index_t my_beg = bw_beg_pos[vert_id];
                        index_t my_end = bw_beg_pos[vert_id+1];
//                        my_work_curr+=my_end-my_beg;

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=bw_csr[my_beg];
                            if(scc_id[vert_id] == 0 && fw_sa[nebr] != -1)
                            {
                                fw_sa[vert_id] = level+1;
                                vertex_frontier++;
//                                front_count++;
                                break;
                            }
                        }
                    }
                }
            }
            else
            {
                index_t *q = new index_t[queue_size];
                index_t head = 0;
                index_t tail = 0;
                //std::queue<index_t> q;

                //Option 1: put current level vertices into fq
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 0 && fw_sa[vert_id] == level)
                    {
                        q[tail++] = vert_id;
//                        front_count ++;
                    }

                }
                while(head != tail)
                {
                    vertex_t temp_v = q[head++];
//                    front_count ++;
                    if(head == queue_size)
                        head = 0;
                    index_t my_beg = fw_beg_pos[temp_v];
                    index_t my_end = fw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        
                        if(scc_id[w] == 0 && fw_sa[w] == -1)
                        {
                            q[tail++] = w;
                            if(tail == queue_size)
                                tail = 0;
                            fw_sa[w] = level + 1;
                        }
                    }
                }
                delete[] q;
//                if(!DEBUG)
//                {
//                    break;
//                }
            }
        
        vertex_front[tid] = vertex_frontier;

        #pragma omp barrier
        vertex_frontier = 0;

        for(index_t i=0; i<thread_count; ++i)
        {
            vertex_frontier += vertex_front[i];
        }
        vertex_visited += vertex_frontier;

        if(VERBOSE)
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
			if(tid==0) 
                std::cout<<"Level-"<<(int)level<<" "
//				<<"-frontier-time-visited:"
				<<vertex_frontier<<" "
                <<fq_size<<" "
                <<(double)(fq_size)/vertex_frontier<<" "
				<<(wtime() - ltm) * 1000<<"ms "
				<<vertex_visited<<" "
                <<edge_frontier<<" "
                <<edge_remaider<<" "
                <<edge_remaider/edge_frontier<<"\n";
        }
        if(vertex_frontier == 0) break;
        
        #pragma omp barrier
        
        if(is_top_down) 
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;
            if((edge_remainder / alpha) < edge_frontier)
            {
                is_top_down = false;
                if(VERBOSE)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to bottom up\n";
                }
            }

        }
        else
            if(!is_top_down && !is_top_down_queue && (fq_size*1.0/beta) > vertex_frontier)
            {
                is_top_down_queue = true;
                if(VERBOSE)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to top down queue\n";
                }
            }
        #pragma omp barrier
        level ++;
    }
}

inline void bw_bfs_fq_queue(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *fw_sa,
        vertex_t *bw_sa,
        index_t *vertex_cur,
        index_t *vertex_front,
        vertex_t root,
        index_t tid,
        index_t thread_count,
        const int alpha,
        const int beta,
        const int gamma,
        vertex_t *frontier_queue,
        vertex_t fq_size,
        const double avg_degree,
        vertex_t vertex_visited,
        vertex_t *temp_queue,
        index_t *prefix_sum,
        vertex_t upper_bound,
        vertex_t *thread_queue
        )
{
    depth_t level = 0;
    bw_sa[root] = 0;
    scc_id[root] = 1;
    temp_queue[0] = root;
    vertex_t queue_size = 1;
//    vertex_t upper_bound = vert_count/thread_count * 10;
//    vertex_t *thread_queue = new vertex_t[upper_bound];
    if(VERBOSE)
    {
        if(tid == 0)
        {
            printf("upperbound, %d\n", upper_bound);
        }
    }
//    vertex_t root_out_degree = fw_beg_pos[root+1] - fw_beg_pos[root];
    bool is_top_down = true;
//    bool is_top_down_async = false;
//    if(VERBOSE)
//    {
//        if(tid == 0)
//        {
//            printf("out_degree, %d, limit, %.3lf\n", root_out_degree, alpha * beta * fq_size);
//        }
//    }
//    if(root_out_degree < alpha * beta * fq_size)
//    {
//        is_top_down_async = true;
//    }
    bool is_top_down_queue = false;
//    index_t queue_size = fq_size / thread_count;
    #pragma omp barrier
    while(true)
    {
        double ltm= wtime();
        vertex_t vertex_frontier = 0;
        
        #pragma omp barrier
        if(is_top_down)
        {
            vertex_t step = queue_size / thread_count;
            vertex_t queue_beg = tid * step;
            vertex_t queue_end = (tid == thread_count - 1 ? queue_size: queue_beg + step);
            
//            printf("tid, %d, %d, %d, %d\n", tid, step, queue_beg, queue_end);
            for(vertex_t q_vert_id=queue_beg; q_vert_id<queue_end; q_vert_id++)
            {
                vertex_t vert_id = temp_queue[q_vert_id];
//                printf("vert_id, %d\n", vert_id);
                //in fq, scc_id[vert_id] is always not 0
                if(scc_id[vert_id] == 1)// && bw_sa[vert_id] == level)
                {
                    index_t my_beg = bw_beg_pos[vert_id];
                    index_t my_end = bw_beg_pos[vert_id+1];
//                    printf("my_beg, %d, my_end, %d\n", my_beg, my_end);
                    for(; my_beg<my_end; my_beg++)
                    {
                        vertex_t nebr=bw_csr[my_beg];
                        if(scc_id[nebr] == 0 && bw_sa[nebr] == -1 && fw_sa[nebr] != -1)
                        {
                            bw_sa[nebr] = level+1;
                            thread_queue[vertex_frontier] = nebr;
                            vertex_frontier++;
                            scc_id[nebr] = 1;
                        }
                    }
                }
            }
        }
        else
            if(!is_top_down_queue)
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 0 && fw_sa[vert_id] != -1)
                    {
                        index_t my_beg = fw_beg_pos[vert_id];
                        index_t my_end = fw_beg_pos[vert_id+1];
//                        my_work_curr+=my_end-my_beg;

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=fw_csr[my_beg];
                            if(scc_id[nebr] == 1)
//                            if(scc_id[nebr] == 1 && bw_sa[nebr] == level)
                            {
                                bw_sa[vert_id] = level+1;
                                vertex_frontier++;
                                scc_id[vert_id] = 1;
//                                front_count++;
                                break;
                            }
                        }
                    }
                }
            }
            else
            {
                vertex_t end_queue = upper_bound;
                index_t head = 0;
                index_t tail = 0;
                //std::queue<index_t> q;
                vertex_t step = queue_size / thread_count;
                vertex_t queue_beg = tid * step;
                vertex_t queue_end = (tid == thread_count - 1 ? queue_size: queue_beg + step);

                //Option 1: put current level vertices into fq
                for(vertex_t q_vert_id=queue_beg; q_vert_id<queue_end; q_vert_id++)
                {
                    thread_queue[tail] = temp_queue[q_vert_id];
                    tail ++;
                }
                while(head != tail)
                {
                    vertex_t temp_v = thread_queue[head++];
//                    front_count ++;
                    if(head == end_queue)
                        head = 0;
                    index_t my_beg = bw_beg_pos[temp_v];
                    index_t my_end = bw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = bw_csr[my_beg];
                        
                        if(scc_id[w] == 0 && bw_sa[w] == -1 && fw_sa[w] != -1)
                        {
                            thread_queue[tail++] = w;
                            if(tail == end_queue)
                                tail = 0;
                            scc_id[w] = 1; 
                            bw_sa[w] = level + 1;
                        }
                    }
                }
            }

        
        vertex_front[tid] = vertex_frontier;
        #pragma omp barrier
        vertex_frontier = 0;

        for(index_t i=0; i<thread_count; ++i)
        {
            vertex_frontier += vertex_front[i];
        }
        vertex_visited += vertex_frontier;
//        #pragma omp barrier
        if(VERBOSE)
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
			if(tid==0 && level < 50) 
                std::cout<<"Level-"<<(int)level<<" "
//				<<"-frontier-time-visited:"
				<<vertex_frontier<<" "
                <<fq_size<<" "
                <<(double)(fq_size)/vertex_frontier<<" "
				<<(wtime() - ltm) * 1000<<"ms "
				<<vertex_visited<<" "
                <<edge_frontier<<" "
                <<edge_remaider<<" "
                <<edge_remaider/edge_frontier<<"\n";
        }
        
        if(vertex_frontier == 0) break;
        
        if(is_top_down) 
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;
            if(!is_top_down_queue && (edge_remainder / alpha) < edge_frontier)
            {
                is_top_down = false;
                if(VERBOSE)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to bottom up\n";
                }
            }

        }
        else
            if((!is_top_down && !is_top_down_queue && (fq_size*1.0/beta) > vertex_frontier) || (!is_top_down && !is_top_down_queue && level > gamma))
                    //            if(level > 10)
            {
                vertex_frontier = 0;
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 1 && bw_sa[vert_id] == level + 1)
                    {
                        thread_queue[vertex_frontier] = vert_id;
                        vertex_frontier++;
                    }
                }
                vertex_front[tid] = vertex_frontier;
           

                is_top_down = false;
                is_top_down_queue = true;
                if(VERBOSE)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to top down queue\n";
                }
            }
//                is_top_down = true;
//                is_top_down_queue = true;
//                if(VERBOSE && tid==0)
//                    std::cout<<"--->Switch back to top down\n";
//                vertex_frontier = 0;
//                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
//                {
//                    vertex_t vert_id = frontier_queue[fq_vert_id];
//                    if(scc_id[vert_id] == 1 && bw_sa[vert_id] == level + 1)
//                    {
//                        thread_queue[vertex_frontier] = vert_id;
//                        vertex_frontier++;
//                    }
//                }
//                vertex_front[tid] = vertex_frontier;
//            }
        
//            if(is_top_down && !is_top_down_queue && (fq_size*1.0/BETA) > vertex_frontier)
// switch to async top down queue

//        if(is_top_down && level > gamma)
//        {
//            if(!is_top_down_queue)
//            {
//                vertex_frontier = 0;
//                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
//                {
//                    vertex_t vert_id = frontier_queue[fq_vert_id];
//                    if(scc_id[vert_id] == 1 && bw_sa[vert_id] == level + 1)
//                    {
//                        thread_queue[vertex_frontier] = vert_id;
//                        vertex_frontier++;
//                    }
//                }
//                vertex_front[tid] = vertex_frontier;
//            }
//
//            is_top_down = false;
//            is_top_down_queue = true;
//            if(VERBOSE)
//            {
//                if(tid==0)
//                    std::cout<<"--->Switch to top down queue\n";
//            }
//        }
        
        #pragma omp barrier

        if(is_top_down || is_top_down_queue)
        {
            get_queue(thread_queue,
                    vertex_front,
                    prefix_sum,
                    tid,
                    temp_queue);
            queue_size = prefix_sum[thread_count-1] + vertex_front[thread_count-1];
//            if(VERBOSE && tid == 0)
//            {
//                printf("queue_size, %d\n", queue_size);
//            }
        }

        #pragma omp barrier
        
        level ++;
    }
//    if(tid == 0)
//        printf("bw_level, %d\n", level);
//    delete[] thread_queue;
}

inline void bw_bfs_fq(
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *fw_sa,
        vertex_t *bw_sa,
        index_t *vertex_cur,
        index_t *vertex_front,
        vertex_t root,
        index_t tid,
        index_t thread_count,
        double alpha,
        double beta,
        vertex_t *frontier_queue,
        vertex_t fq_size,
        const double avg_degree,
        vertex_t vertex_visited
        )
{
    bw_sa[root] = 0;
    vertex_t root_in_degree = bw_beg_pos[root+1] - bw_beg_pos[root];
    bool is_top_down = true;
    bool is_top_down_queue = false;
    bool is_top_down_async = false;
    if(DEBUG)
    {
        if(tid == 0)
        {
            printf("in_degree, %d, limit, %.3lf\n", root_in_degree, alpha * beta * fq_size);
        }
    }
    if(root_in_degree < alpha * beta * fq_size)
    {
        is_top_down_async = true;
    }
    index_t level = 0;
    scc_id[root] = 1;
    index_t queue_size = fq_size / thread_count;
    while(true)
    {
        double ltm= wtime();
//        index_t front_count=0;
//        index_t my_work_next=0;
//        index_t my_work_curr=0;
        vertex_t vertex_frontier = 0;

        if(is_top_down)
        {
            if(is_top_down_async)
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 1 && (bw_sa[vert_id] == level || bw_sa[vert_id] == level + 1))
                    {
                        index_t my_beg = bw_beg_pos[vert_id];
                        index_t my_end = bw_beg_pos[vert_id+1];

    //                        printf("my_beg = %d, my_end = %d\n", my_beg, my_end);
                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=bw_csr[my_beg];
                            if(scc_id[nebr] == 0 && bw_sa[nebr] == -1 && fw_sa[nebr] != -1)
                            {
                                bw_sa[nebr] = level+1;
//                                my_work_next+=bw_beg_pos[nebr+1]-bw_beg_pos[nebr];
                                vertex_frontier++;
                                scc_id[nebr] = 1;	
                            }
                        }
                    }
                }
            }
            else
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 1 && bw_sa[vert_id] == level)
                    {
                        index_t my_beg = bw_beg_pos[vert_id];
                        index_t my_end = bw_beg_pos[vert_id+1];

    //                        printf("my_beg = %d, my_end = %d\n", my_beg, my_end);
                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=bw_csr[my_beg];
                            if(scc_id[nebr] == 0 && bw_sa[nebr] == -1 && fw_sa[nebr] != -1)
                            {
                                bw_sa[nebr] = level+1;
//                                my_work_next+=bw_beg_pos[nebr+1]-bw_beg_pos[nebr];
                                vertex_frontier++;
//                                front_count++;
                                scc_id[nebr] = 1;	
                            }
                        }
                    }
                }

            }
//            work_comm[tid]=my_work_next;
        }
        else
            if(!is_top_down_queue)
            {
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
//                    if(scc_id[vert_id] == 0 && bw_sa[vert_id] == -1 && fw_sa[vert_id] != -1)
                    if(scc_id[vert_id] == 0)
                    {
                        index_t my_beg = fw_beg_pos[vert_id];
                        index_t my_end = fw_beg_pos[vert_id+1];
//                        my_work_curr+=my_end-my_beg;

                        for(; my_beg<my_end; my_beg++)
                        {
                            vertex_t nebr=fw_csr[my_beg];
//                            if(scc_id[nebr] == 1 && bw_sa[nebr] != -1 && fw_sa[nebr] != -1)
                            if(scc_id[nebr] == 1)
                            {
                                bw_sa[vert_id] = level+1;
//                                front_count++;
                                vertex_frontier++;
                                scc_id[vert_id] = 1;
                                break;
                            }
                        }
                    }
                }
//                work_comm[tid]=my_work_curr;
            }
            else
            {
//                std::queue<index_t> q;
                index_t *q = new index_t[queue_size];
                index_t head = 0;
                index_t tail = 0;
                
                for(vertex_t fq_vert_id=vert_beg; fq_vert_id<vert_end; fq_vert_id++)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(scc_id[vert_id] == 1 && bw_sa[vert_id] == level)
                    {
                        q[tail++] = vert_id;
//                        front_count ++;
                        //impossible here
//                        if(tail == queue_size)
//                            tail = 0;
                    }
                }

                while(head != tail)
                {
                    vertex_t temp_v = q[head++];
//                    front_count ++;
                    if(head == queue_size)
                        head = 0;
                    index_t my_beg = bw_beg_pos[temp_v];
                    index_t my_end = bw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = bw_csr[my_beg];
                        
                        if(scc_id[w] == 0 && bw_sa[w] == -1 && fw_sa[w] != -1)
                        {
                            q[tail++] = w;
                            if(tail == queue_size)
                                tail = 0;
                            scc_id[w] = 1;
                            bw_sa[w] = level + 1;
                        }
                    }
                }
                delete[] q;
                if(!DEBUG)
                    break;

            }

//        front_comm[tid]=front_count;

        vertex_front[tid] = vertex_frontier;

        #pragma omp barrier
        vertex_frontier = 0;

        for(index_t i=0; i<thread_count; ++i)
        {
            vertex_frontier += vertex_front[i];
        }
        vertex_visited += vertex_frontier;

        if(VERBOSE)
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remaider = (double)(fq_size - vertex_visited) * avg_degree;
			if(tid==0) 
                std::cout<<"Level-"<<(int)level<<" "
//				<<"-frontier-time-visited:"
				<<vertex_frontier<<" "
                <<fq_size<<" "
                <<(double)(fq_size)/vertex_frontier<<" "
				<<(wtime() - ltm) * 1000<<"ms "
				<<vertex_visited<<" "
                <<edge_frontier<<" "
                <<edge_remaider<<" "
                <<edge_remaider/edge_frontier<<"\n";
        }
        if(vertex_frontier == 0) break;
        
        #pragma omp barrier
        
        if(is_top_down) 
        {
            double edge_frontier = (double)vertex_frontier * avg_degree;
            double edge_remainder = (double)(fq_size - vertex_visited) * avg_degree;
            if((edge_remainder / alpha) < edge_frontier)
            {
                is_top_down = false;
                if(VERBOSE)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to bottom up\n";
                }
            }

        }
        else
            if(!is_top_down && !is_top_down_queue && (fq_size*1.0/beta) > vertex_frontier)
            {
                is_top_down_queue = true;
                if(VERBOSE)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to top down queue\n";
                }
            }
        #pragma omp barrier
        level ++;
/*        #pragma omp barrier
        front_count=0;
        my_work_next=0;

        for(index_t i=0;i<thread_count;++i)
        {
            front_count += front_comm[i];
            my_work_next += work_comm[i];
        }

        if(DEBUG)
        {
            if(tid==0) 
                std::cout<<"Level-"<<(int)level
                <<"-frontier-time-futurework(currwork in btup): "
                <<front_count<<" "
                <<(wtime() - ltm) * 1000<<" ms "
                <<my_work_next<<"\n";
        }
            
        if(front_count == 0) 
        {
            break;
        }
        
        if(is_top_down && my_work_next > (alpha*fq_size)) 
        {
            is_top_down=false;
            if(DEBUG)
            {
                if(tid==0)
                    std::cout<<"--->Switch to bottom up"<<my_work_next<<" "<<fq_size<<"<----\n";
            }
        }
        else
            if(!is_top_down && !is_top_down_queue && front_count < (beta * fq_size))
            {
                is_top_down_queue = true;
                if(DEBUG)
                {
                    if(tid==0)
                        std::cout<<"--->Switch to top down queue"<<my_work_next<<" "<<fq_size<<"<----\n";
                }
            }
        #pragma omp barrier

        level ++;
        */
    }
    if(tid == 0)
        printf("bw_level, %d\n", level);
}

inline void init_fw_sa(
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_sa,
        vertex_t *frontier_queue,
        vertex_t fq_size,
        vertex_t *wcc_fq,
        index_t tid,
        vertex_t *color,
        vertex_t &wcc_fq_size
        )
{
    for(index_t i=vert_beg; i<vert_end; ++i)
    {
        fw_sa[i] = -1;
    }
    if(tid == 0)
    {
        std::set<int> s_fq;
        for(vertex_t i=0; i<fq_size; ++i)
        {
//            vertex_t v = frontier_queue[i];
//            vertex_t wcc_id = color[v];
//            if(s_fq.find(v) == s_fq.end()) 
            s_fq.insert(color[frontier_queue[i]]);
        }
        wcc_fq_size = s_fq.size();
        std::set<int>::iterator it;
        int i=0;
        for(it = s_fq.begin(); it != s_fq.end(); ++it, ++i)
        {
//            printf("%d", *it);
            wcc_fq[i] = *it;
        }
    }
}

inline void mice_fw_bw(
        color_t *wcc_color,
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        vertex_t *fw_sa,
        index_t tid,
        index_t thread_count,
        vertex_t *frontier_queue,
        vertex_t fq_size,
        vertex_t *wcc_fq,
        vertex_t wcc_fq_size
        )
{
    index_t step = wcc_fq_size / thread_count;
    index_t wcc_beg = tid * step;
    index_t wcc_end = (tid == thread_count - 1 ? wcc_fq_size : wcc_beg + step);
    
    index_t *q = new index_t[fq_size];
    index_t head = 0;
    index_t tail = 0;

    for(index_t fq_i=0; fq_i<fq_size; ++fq_i)
    {
        vertex_t v = frontier_queue[fq_i];
        if(scc_id[v] == 0)
        {
            vertex_t cur_wcc = wcc_color[v];
            bool in_wcc = false;
            for(index_t i=wcc_beg; i<wcc_end; ++i)
            {
                if(wcc_fq[i] == cur_wcc)
                {
                    in_wcc = true;
                    break;
                }
            }
            if(in_wcc)
            {
//                if(tid == 0)
//                {
//                    printf("v, %d, wcc, %d\n", v, wcc_color[v]);
//                }
                //fw
                fw_sa[v] = v;
                q[tail++] = v;
                if(tail == fq_size)
                    tail = 0;
                while(head != tail)
                {
                    vertex_t temp_v = q[head++];
    //                    front_count ++;
                    if(head == fq_size)
                        head = 0;
                    index_t my_beg = fw_beg_pos[temp_v];
                    index_t my_end = fw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = fw_csr[my_beg];
                        
                        if(scc_id[w] == 0 && fw_sa[w] != v)
                        {
                            q[tail++] = w;
                            if(tail == fq_size)
                                tail = 0;
                            fw_sa[w] = v;
                        }
                    }
                }
                //bw
                scc_id[v] = v;
                q[tail++] = v;
                if(tail == fq_size)
                    tail = 0;

                while(head != tail)
                {
                    vertex_t temp_v = q[head++];
    //                    front_count ++;
                    if(head == fq_size)
                        head = 0;
                    index_t my_beg = bw_beg_pos[temp_v];
                    index_t my_end = bw_beg_pos[temp_v+1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        index_t w = bw_csr[my_beg];
                        
                        if(scc_id[w] == 0 && fw_sa[w] == v)
                        {
                            q[tail++] = w;
                            if(tail == fq_size)
                                tail = 0;
                            scc_id[w] = v;
                        }
                    }
                }
                
            }
            
        }
    }
    delete[] q;
}
#endif
