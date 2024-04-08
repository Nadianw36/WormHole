#!/usr/bin/bash


#compile code
g++ main.cpp -fopenmp -std=c++17 -O2 -g -o mll

#create dir to hold indexes
if [ ! -d indexes ]; then
  mkdir indexes
fi

cd ..
exe/createL0COO {GRAPH_NAME}_seed_{SIZE}
exe/sample L0-BiBFS {GRAPH_NAME}_seed_{SIZE} 0 10000
exe/queryOnCore {GRAPH_NAME}_seed_{SIZE} 0 10000 #run BiBFS with RandomL0 and RandomL0FromRandomL1 queries
#create index for toyG
# arg0: sub-command
# -d: specify graph name, graph should be put in graphs/
cd MLL

./mll index -d {GRAPH_NAME}_seed_{SIZE}

#answering queries
# arg0: sub-command
# -d: specify graph name, graph should be put in graphs/
# -q: specify path to queries
./mll query -d {GRAPH_NAME}_seed_{SIZE} -q ../results/{GRAPH_NAME}/L0/BiBFS/results/livejournal/L0/BiBFS/0_10000_{GRAPH_NAME}_seed_{SIZE}_L0-BiBFS_core-queries_random-L0.txt
./mll query -d {GRAPH_NAME}_seed_{SIZE} -q ../results/{GRAPH_NAME}/L0/BiBFS/results/livejournal/L0/BiBFS/0_10000_{GRAPH_NAME}_seed_{SIZE}_L0-BiBFS_core-queries_random-L0-RANDOM-L1.txt
