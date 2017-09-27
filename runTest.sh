#!/bin/bash
#checking kernels
#nohup nvprof --force-overwrite --log-file ReadableOutput.log --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 1 STDKernelV 1 >Results.txt & disown
#nohup nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 0 STDKernelV 0 >Results1.txt & disown
#nohup nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 1 STDKernelV 1 >Results2.txt & disown
./main CUDA Nelem 3000 N 1 Nepad 32 KernelV 10 STDKernelV 7 >Results.txt
./main CUDA Nelem 2000 N 2 Nepad 16 KernelV 10 STDKernelV 7 >Results.txt