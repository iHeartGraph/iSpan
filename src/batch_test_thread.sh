alpha=0.1
beta=0.0001
thread_count=32
run_times=10
async_limit=1
file[1]="livejournal"
file[2]="baidu"
file[3]="flickr-growth"
file[4]="zhishi-hudong-relatedpages"
file[5]="RMAT"
file[6]="wikipedia_link_en"
file[7]="wikipedia_en"
file[8]="dbpedia"
file[9]="facebook"
file[10]="pokec"
file[11]="random"
file[12]="twitter_mpi"
file[13]="twitter_www"

deal () {
    ./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $alpha $beta $2 $run_times $async_limit > /mnt/raid0_huge/yuede/data/scc_result/${file[$1]}/thread_count_${2}.csv

}
make clean
make
for index in `seq 13 13`;
do
    echo $index
    echo ${file[$index]}
#    mkdir /mnt/raid0_huge/yuede/data/scc_result/${file[$index]}
#    deal $index
    for i in `seq 32 2 128`;
    do
        thread_count=$i
        #$((1<<$i))
        echo $thread_count
        deal $index $thread_count
    done
done
#./scc_cpu /mnt/raid0_huge/yuede/data/livejournal/fw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/fw_adjacent.bin /mnt/raid0_huge/yuede/data/livejournal/bw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/bw_adjacent.bin $alpha $beta $thread_count


