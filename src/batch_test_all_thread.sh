alpha=0.1
beta=0.002
#thread_count=48
run_times=10
async_limit=1
file[1]="livejournal"
file[2]="baidu"
file[3]="flickr-growth"
file[4]="zhishi-hudong-relatedpages"
file[5]="pokec"
file[6]="wikipedia_link_en"
file[7]="wikipedia_en"
file[8]="dbpedia"
file[9]="facebook"
#file[10]="twitter_mpi"
file[10]="wiki_communication"
file[11]="twitter_www"
file[12]="wiki_talk_en"
file[13]="random"
file[14]="RMAT"
#file[14]="us_patent"
#file[15]="citeseer"
#file[16]="scale25"

deal () {
    #./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $alpha $beta $thread_count $run_times $async_limit > result_${file[$1]}.csv
#    echo $1
#    echo $2
#    mkdir /mnt/raid0_huge/yuede/data/scc_result_fisc/${file[$1]}
    ./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $alpha $beta $2 $run_times $async_limit > /mnt/raid0_huge/yuede/data/scc_result_fisc/${file[$1]}/thread_${2}.csv

}
make clean
make
for index in `seq 1 14`;
do
    echo $index
    echo ${file[$index]}
    for thread_num in 1 2 4 8 16 32 56
    do
        echo $thread_num
        deal $index $thread_num
    done

#    for thread_num in `seq 1 64`;
#    do
#        if [ $thread_num -lt 3 ] || [ $thread_num == 4 ]
#        then
#            echo $thread_num
#            deal $index $thread_num
#        else
#            if [ $thread_num > 7 ] && [ $(($thread_num % 8)) == 0 ]
#            then
#                echo $thread_num
#                deal $index $thread_num
#            fi
#        else
#            if [ $thread_num > 48 ] && [ $(($thread_num % 2)) == 0 ]
#            then 
#                echo $thread_num
#                deal $index $thread_num
#            fi
#        fi
#    done
#    deal $index
done
#./scc_cpu /mnt/raid0_huge/yuede/data/livejournal/fw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/fw_adjacent.bin /mnt/raid0_huge/yuede/data/livejournal/bw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/bw_adjacent.bin $alpha $beta $thread_count


