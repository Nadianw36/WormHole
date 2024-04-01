#!/usr/bin/bash

#compile code
g++ main.cpp -fopenmp -std=c++17 -O2 -g -o mll

#create dir to hold indexes
if [ ! -d indexes ]; then
  mkdir indexes
fi

#create index for toyG
# arg0: sub-command
# -d: specify graph name, graph should be put in graphs/
./mll index -d '$0'

#answering queries
# arg0: sub-command
# -d: specify graph name, graph should be put in graphs/
# -q: specify path to queries
./mll query -d '$0' -q '../inputs/$0_input.txt'
