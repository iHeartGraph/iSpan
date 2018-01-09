alpha=30
beta=200
gamma=10
theta=0.01
thread_count=56
run_times=20
async_limit=1
thread_num=56

file[1]="livejournal"
file[2]="baidu"
file[3]="flickr-growth"
file[4]="zhishi-hudong-relatedpages"
file[5]="pokec"
file[6]="wikipedia_link_en"
file[7]="wikipedia_en"
file[8]="dbpedia"
file[9]="facebook"
file[10]="twitter_mpi"
file[11]="wiki_talk_en"
file[12]="wiki_communication"
#file[14]="us_patent"
#file[15]="citeseer"
#file[16]="scale25"
file[13]="random"
file[14]="RMAT"
file[15]="twitter_www"

deal () {

#    mkdir /mnt/raid0_huge/yuede/data/scc_result_scalability/${file[$1]}
    ./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $2 $alpha $beta $gamma $theta $run_times > /mnt/raid0_huge/yuede/scc_result/scc_time_result/${file[$1]}.csv
#    echo $2
#    ./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $thread_count $alpha $beta $2 $theta $run_times > test_gamma/result_${file[$1]}_${2}.csv

}
make clean
make
for index in `seq 1 15`;
do
    echo $index
    echo ${file[$index]}
#    for thread_num in 1 2 4 8 16 32 56
#    do
#        echo $thread_num
    deal $index $thread_num
#        sleep 2 
#    done
#    for a in `seq 5 30`;
#    do
#        echo $a
#        deal $index $a
#    done
done
#./scc_cpu /mnt/raid0_huge/yuede/data/livejournal/fw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/fw_adjacent.bin /mnt/raid0_huge/yuede/data/livejournal/bw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/bw_adjacent.bin $alpha $beta $thread_count


