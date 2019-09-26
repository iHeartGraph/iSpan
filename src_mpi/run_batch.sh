#!/bin/sh
# one hour timelimit:
#SBATCH --time 01:00:00
#SBATCH -o friendster_stat_time_1_2.csv

# default queue, 16 processors (two nodes worth)
# -p: which queue

#SBATCH -p 256gb -N 1 -n 2

###########################################
# -n means number of processes
# -N means number of machines 
# salloc -t 100 -p allgpu-noecc: ask for interactive node @ debug queue
# sinfo: shows the queues status
# scancel jid: cancel the job you do not want to continue
# squeue | grep username: shows my allocated node
#################################

# sleep 3600
# You can also put them in ~/.bashrc
module load openmpi/current
# module load cuda/toolkit/7.5

# num=16

graph_folder=/lustre/groups/huanglab/yuede/graph
# graph=twitter_www
# graph=wikipedia_en
# graph=baidu
# graph=dbpedia
# graph=flickr-growth
# graph=zhishi-hudong-relatedpages
# graph=pokec
# graph=wiki_talk_en
# mpirun -n $num ./scc_cpu $graph_folder/livejournal/fw_begin.bin $graph_folder/livejournal/fw_adjacent.bin $graph_folder/livejournal/bw_begin.bin $graph_folder/livejournal/bw_adjacent.bin 56 30 200 10 0.01 10

# file[1]="livejournal"
# file[2]="flickr-growth"
# file[3]="wikipedia_en"
# file[4]="random"

# file[1]="twitter_www"
# file[1]="twitter_mpi"

 file[1]="livejournal"
# file[2]="flickr-growth"
# file[3]="wikipedia_en"
# file[4]="random"
# file[5]="twitter_www"
# file[6]="twitter_mpi"
# file[7]="wiki_talk_en"
# file[8]="dbpedia"
# file[9]="facebook"
# file[10]="pokec"
# file[11]="baidu"
# file[12]="zhishi-hudong-relatedpages"
# #file[12]="twitter_www"
# #file[12]="twiter_mpi"
# file[13]="RMAT"
# file[14]="wikipedia_link_en"
# file[15]="wiki_communication"
#file[15]="friendster"

# for i in `seq 1 15`;
for i in `seq 1 1`;
do
    graph=${file[$i]}
    echo $graph
    #mpirun -N 2 ./scc_cpu $graph_folder/$graph/fw_begin.bin $graph_folder/$graph/long_int_back/fw_adjacent.bin $graph_folder/$graph/bw_begin.bin $graph_folder/$graph/long_int_back/bw_adjacent.bin 56 30 200 10 0.01 5
    mpirun -N 2 ./scc_cpu $graph_folder/$graph/fw_begin.bin $graph_folder/$graph/fw_adjacent.bin $graph_folder/$graph/bw_begin.bin $graph_folder/$graph/bw_adjacent.bin 56 30 200 10 0.01 5
    wait
done


