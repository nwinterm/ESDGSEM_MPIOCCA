#!/bin/bash
#checking kernels
nvprof -o Output.csv --force-overwrite --csv --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput nohup ./main CUDA >Results.txt & disown
