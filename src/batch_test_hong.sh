#thread_count=32
run_times=10
edge_list="edge_list_hong.bin"
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
file[12]="twiter_mpi"
file[13]="twitter_www"
deal () {
    for i in `seq 1 $run_times`;
    do
#        echo $i
#        echo ${2}_${i}.csv
        ./scc /mnt/raid0_huge/yuede/data/${file[$1]}/$edge_list $2 2 -d >> /mnt/raid0_huge/yuede/data/scc_hong_result/${file[$1]}/thread_count_${2}.csv
    done
#    ./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $alpha $beta $thread_count $run_times $async_limit > result_${file[$1]}.csv

}
make clean
make
for index in `seq 1 13`;
do
#    mkdir result/${file[$index]}
    echo $index
    echo ${file[$index]}
#    mkdir /mnt/raid0_huge/yuede/data/scc_hong_result/${file[$index]}
#    for i in `seq 0 8`;
    for i in `seq 32 2 128`;
    do
        thread_count=$i
        #$((1<<$i))
        echo $thread_count
        deal $index $thread_count
    done
done


