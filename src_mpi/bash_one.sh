make clean
make
mpirun -N 2 ./scc_cpu /mnt/raid0_huge/yuede/data/friendster/fw_begin.bin /mnt/raid0_huge/yuede/data/friendster/fw_adjacent.bin /mnt/raid0_huge/yuede/data/friendster/bw_begin.bin /mnt/raid0_huge/yuede/data/friendster/bw_adjacent.bin 56 30 200 10 0.01 1
