alpha=0.1
beta=0.002
thread_count=56
run_times=10
async_limit=1
#file[1]="livejournal"
#file[2]="baidu"
#file[3]="flickr-growth"
#file[4]="zhishi-hudong-relatedpages"
#file[5]="pokec"
file[1]="wikipedia_link_en"
#file[7]="wikipedia_en"
#file[8]="dbpedia"
#file[9]="facebook"
#file[10]="twitter_mpi"
file[2]="twitter_www"
file[3]="wiki_talk_en"
#file[13]="wiki_communication"
#file[14]="us_patent"
#file[15]="citeseer"
#file[16]="scale25"
#file[15]="random"
file[4]="RMAT"

deal () {
    ./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $alpha $beta $thread_count $run_times $async_limit > result_${file[$1]}.csv

}
make clean
make
for index in `seq 1 4`;
do
    echo $index
    echo ${file[$index]}
    deal $index
done
#./scc_cpu /mnt/raid0_huge/yuede/data/livejournal/fw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/fw_adjacent.bin /mnt/raid0_huge/yuede/data/livejournal/bw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/bw_adjacent.bin $alpha $beta $thread_count


