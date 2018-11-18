exe = ispan

cc = "$(shell which g++)" 
#flags = -I. -fopenmp -march=athlon64 -O3
flags = -I. -fopenmp -O3
#flags += -std=c++11

ifeq ($(debug), 1)
	flags+= -DDEBUG 
endif

objs = $(patsubst %.cpp,%.o,$(wildcard ../../lib/*.cpp))\
			$(patsubst %.cpp,%.o,$(wildcard *.cpp))

deps = $(wildcard ../../lib/*.h) \
	   			$(wildcard ./*.hpp) \
				$(wildcard *.h) \
				Makefile

%.o:%.cpp $(deps)
	$(cc) -c $< -o $@ $(flags)

$(exe):$(objs)
	$(cc) $(objs) -o $(exe) $(flags)


test:$(exe)
	./scc_cpu /mnt/raid0_huge/yuede/data/pokec/fw_begin.bin /mnt/raid0_huge/yuede/data/pokec/fw_adjacent.bin /mnt/raid0_huge/yuede/data/pokec/bw_begin.bin /mnt/raid0_huge/yuede/data/pokec/bw_adjacent.bin 

clean:
	rm -rf $(exe) $(objs) 
