#include <map>
#include "graph.h"
#include "util.h"
#include "wtime.h"

graph *graph_load(
        const char *fw_beg_file, 
        const char *fw_csr_file, 
        const char *bw_beg_file, 
        const char *bw_csr_file,
        double *avg_time
        )
{
    graph *g = new graph(fw_beg_file,
            fw_csr_file,
            bw_beg_file,
            bw_csr_file);

    for(index_t i=0; i<15; ++i)
        avg_time[i] = 0.0;
    return g;
}

void get_scc_result(
        index_t *scc_id,
        index_t vert_count
        )
{
    index_t size_1 = 0;
    index_t size_2 = 0;
    index_t size_3_type_1 = 0;
    index_t size_3_type_2 = 0;
    index_t largest = 0;
    std::map<index_t, index_t> mp;
    for(index_t i=0; i<vert_count+1; ++i)
    {
        if(scc_id[i] == 1)
            largest ++;
        else
            if(scc_id[i] == -1)
                size_1 ++;
            else
                if(scc_id[i] == -2)
                    size_2 ++;
                else
                    if(scc_id[i] == -3)
                    {
                        size_3_type_1 ++;
                    }
                    else
                        if(scc_id[i] == -4)
                            size_3_type_2 ++;
                        else
                        {
                            mp[scc_id[i]] ++;
                        }
    }
    printf("\nResult:\nlargest, %d\ntrimmed size_1, %d\ntrimmed size_2, %d\ntrimmed size_3, %d\nothers, %d\ntotal, %d\n", largest, size_1, size_2/2, size_3_type_1/3 +size_3_type_2/3, mp.size(), (1 + size_1 + size_2/2 + size_3_type_1/3 + size_3_type_2/3 + mp.size()));

}

void print_time_result(
        index_t run_times,
        double *avg_time
        )
{
    //0 trim, 1 largest SCC, 2 small SCC, 3 total time
    //4 trim_size_1, 5 trim_size_2, 6 pivot_selection, 7 fw_bfs, 8 bw_bfs, 9 color propagation, 10 color identify, 11 color_init
    if(run_times > 0)
    {
        for(index_t i=0; i<15; ++i)
            avg_time[i] = (avg_time[i] / run_times) * 1000;
        printf("\nAverage Time Consumption for Running %d Times (ms)\n", run_times);
        printf("Trim, %.3lf\n", avg_time[0]);
        printf("Elephant SCC, %.3lf\n", avg_time[1]);
//        printf("WCC, %.3lf\n", avg_time[12]);
        printf("Mice SCC, %.3lf\n", avg_time[2]);
        printf("Total time, %.3lf\n", avg_time[3]);

        printf("\n------Details------\n");
        printf("Trim size_1, %.3lf\n", avg_time[4]);
        printf("Pivot selection, %.3lf\n", avg_time[6]);
        printf("FW BFS, %.3lf\n", avg_time[7]);
        printf("BW BFS, %.3lf\n", avg_time[8]);
        printf("Trim size_2, %.3lf\n", avg_time[5]);
        printf("Trim size_3, %.3lf\n", avg_time[13]);
//        printf("GFQ, %.3lf\n", avg_time[14]);
        printf("Wcc, %.3lf\n", avg_time[9]);
        printf("Mice fw-bw, %.3lf\n", avg_time[10]);
       // printf("Color init, %.3lf\n", avg_time[11]);
    }
}
