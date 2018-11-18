#include "wtime.h"
#include "graph.h"
//#include "frontier_queue.h"
#include "scc_common.h"
#include <unistd.h>

int main(int args, char **argv)
{
    printf("args = %d\n", args);
	if(args != 11)
    {
	    std::cout<<"Usage: ./scc_cpu <fw_beg_file> <fw_csr_file> <bw_beg_file> <bw_csr_file> <thread_count> <alpha> <beta> <gamma> <theta> <run_times>\n";
        exit(-1);
    }
	
	const char *fw_beg_file = argv[1];
	const char *fw_csr_file = argv[2];
    const char *bw_beg_file = argv[3];
    const char *bw_csr_file = argv[4];
	const index_t thread_count=atoi(argv[5]);
    const int alpha = atof(argv[6]);
    const int beta = atof(argv[7]);
    const int gamma = atof(argv[8]);
    const double theta= atof(argv[9]);
    const index_t run_times = atoi(argv[10]);
	printf("Thread = %d, alpha = %d, beta = %d, gamma = %d, theta = %g, run_times = %d\n", thread_count, alpha, beta, gamma, theta, run_times);
    
    double * avg_time = new double[15]; 
    ///step 1: load the graph
    graph *g = graph_load(fw_beg_file,
                    fw_csr_file, 
                    bw_beg_file, 
                    bw_csr_file,
                    avg_time);
    index_t i=0;

    ///step 2: detect scc
    while(i++ < run_times)
    {
        printf("\nRuntime: %d\n", i);
        scc_detection(g,
                alpha, 
                beta, 
                gamma, 
                theta,
                thread_count,
                avg_time);
//        sleep(2);

    }
//    printf("ready?\n");
    ///step 3: print result
    print_time_result(run_times,
            avg_time);
//    delete[] avg_time;
    return 0;
}
