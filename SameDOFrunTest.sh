#!/bin/bash
#checking kernels
#nohup nvprof --force-overwrite --log-file ReadableOutput.log --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 1 STDKernelV 1 >Results.txt & disown
#nohup nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 0 STDKernelV 0 >Results1.txt & disown
#nohup nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 1 STDKernelV 1 >Results2.txt & disown
VerFD=7
VerSD=7
OutFile=Results.txt
./main CUDA Nelem 3000 N 1 Nepad 32 KernelV $VerFD STDKernelV $VerSD &>$OutFile 
./main CUDA Nelem 2000 N 2 Nepad 16 KernelV $VerFD STDKernelV $VerSD &>>$OutFile 
./main CUDA Nelem 1500 N 3 Nepad 4 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
./main CUDA Nelem 1200 N 4 Nepad 16 KernelV $VerFD STDKernelV $VerSD &>>$OutFile 
./main CUDA Nelem 1000 N 5 Nepad 4 KernelV $VerFD STDKernelV $VerSD &>>$OutFile 
./main CUDA Nelem 858 N 6 Nepad 5 KernelV $VerFD STDKernelV $VerSD &>>$OutFile 
./main CUDA Nelem 750 N 7 Nepad 1 KernelV $VerFD STDKernelV $VerSD &>>$OutFile 
./main CUDA Nelem 667 N 8 Nepad 3 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
./main CUDA Nelem 600 N 9 Nepad 4 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
./main CUDA Nelem 546 N 10 Nepad 1 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
./main CUDA Nelem 500 N 11 Nepad 2 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
./main CUDA Nelem 462 N 12 Nepad 2 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
./main CUDA Nelem 429 N 13 Nepad 2 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
./main CUDA Nelem 400 N 14 Nepad 2 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
./main CUDA Nelem 375 N 15 Nepad 2 KernelV $VerFD STDKernelV $VerSD &>>$OutFile
