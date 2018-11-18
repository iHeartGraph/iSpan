#include <iostream>
#include <unistd.h>
#include <omp.h>
#include "graph.h"
#include "wtime.h"
#include <assert.h>
#include <algorithm>
#include "util.h"

template<typename vertex_t,
	typename index_t,
	typename color_t>
void 
openmp_wcc(
		vertex_t *csr,
		index_t *beg_pos,
		color_t *vert_color,
		color_t *color_redirect,
		bool *is_redirect,
		bool *is_change,
		color_t *global_color,
		index_t tid,
		index_t thread_count,
		index_t vert_beg,
		index_t vert_end,
		index_t vert_count,
		index_t edge_count,
        index_t *scc_id,
        index_t *small_queue
){
	index_t level=0;
	color_t global_color_beg=global_color[0];

	while(true)
	{
        if(DEBUG)
        {
		    if(tid==0)std::cout<<"Iteration: "<<level++<<"\n";
        }
		//if(tid==0) global_color=0;
#pragma omp barrier

		bool is_redirect_local=false;
		bool is_change_local=false;

		double tm=wtime();
		for(vertex_t fq_v_id=vert_beg; fq_v_id<vert_end; fq_v_id++)
		{
            vertex_t vert_id = small_queue[fq_v_id];

			if(vert_color[vert_id]!=NEGATIVE && 
					vert_color[vert_id]<global_color_beg) continue;
			index_t my_beg = beg_pos[vert_id];
			index_t my_end = beg_pos[vert_id+1];
			for(; my_beg<my_end; my_beg++)
			{
				vertex_t dest=csr[my_beg];
/// for scc
                if(scc_id[dest] != 0)
                    continue;
//                printf("dest, %d\n", dest);

				if(vert_color[vert_id]==NEGATIVE && vert_color[dest]==NEGATIVE)
				{
					if(!is_change_local) is_change_local=true;

					color_t color=__sync_fetch_and_add(global_color,1);
					vert_color[vert_id]=color;
					vert_color[dest]=color;
				}
				else if(vert_color[vert_id]!=NEGATIVE && vert_color[dest]==NEGATIVE)
				{
					if(!is_change_local) is_change_local=true;
					vert_color[dest]=vert_color[vert_id];	
				}
				else if(vert_color[vert_id]==NEGATIVE && vert_color[dest]!=NEGATIVE)
				{
					if(!is_change_local) is_change_local=true;
					vert_color[vert_id]= vert_color[dest];
				}
				else
				{
					color_t src_color=vert_color[vert_id];
					color_t dest_color=vert_color[dest];

					bool ret=false;
					do{
						while(src_color!=color_redirect[src_color])
							src_color=color_redirect[src_color];

						while(dest_color!=color_redirect[dest_color])
							dest_color=color_redirect[dest_color];

						if(src_color==dest_color)
						{
							break;
						}

						if(!is_change_local) is_change_local=true;
						if(!is_redirect_local) is_redirect_local=true;

						if(src_color>dest_color)
						{
							//possibly, the larger color is just updated.
							//thereby, we have to check again.
							ret=__sync_bool_compare_and_swap(color_redirect+
									src_color,src_color,dest_color);
						}
						else
						{
							ret=__sync_bool_compare_and_swap(color_redirect+
									dest_color,dest_color,src_color);
						}
					}while(ret==false);
				}
			}
		}
		is_redirect[tid]=is_redirect_local;
		is_change[tid]=is_change_local;

#pragma omp barrier
		for(index_t i=0;i<thread_count;i++)
			if(is_redirect[i]) 
			{
				is_redirect_local=true;
				break;
			}

		for(index_t i=0;i<thread_count;i++)
			if(is_change[i]) 
			{
				is_change_local=true;
				break;
			}

        if(DEBUG)
        {
		    if(tid==0) std::cout<<"Global color: "<<global_color[0]
			    <<" "<<wtime()-tm<<" seconds\n";
        }

		//redirect color updates
		if(tid==0){
			index_t color_count=0;
			for(index_t i=global_color_beg;i<global_color[0];i++)
			{
				if(color_redirect[i]==i)
					color_count++;
			}
            if(DEBUG)
            {
			    std::cout<<"WCC groups: "<<color_count<<"\n";
            }
			//global_color=color_count;
		}


		if((is_change_local==false) && (is_redirect_local==false))
        {
            //printf("global color, %d\n", global_color[0]);
            break;
        }
	}
	return ;	
}
