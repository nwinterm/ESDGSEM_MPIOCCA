#!/bin/bash
#checking kernels
#nohup nvprof --force-overwrite --log-file ReadableOutput.log --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 1 STDKernelV 1 >Results.txt & disown
#nohup nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 0 STDKernelV 0 >Results1.txt & disown
#nohup nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 1 STDKernelV 1 >Results2.txt & disown
StartVerFD=0
EndVerFD=7
VerSD=7
OutFile=Results.txt
rm -f $OutFile
for ((currVersion=$StartVerFD;currVersion<=$EndVerFD;currVersion++)){
#echo $currVersion
	./main CUDA Nelem 3000 N 1 Nepad 32 KernelV $currVersion STDKernelV $VerSD &>>$OutFile 
	./main CUDA Nelem 2000 N 2 Nepad 16 KernelV $currVersion STDKernelV $VerSD &>>$OutFile 
	./main CUDA Nelem 1750 N 3 Nepad 4 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
	./main CUDA Nelem 1300 N 4 Nepad 16 KernelV $currVersion STDKernelV $VerSD &>>$OutFile 
	./main CUDA Nelem 1100 N 5 Nepad 4 KernelV $currVersion STDKernelV $VerSD &>>$OutFile 
	./main CUDA Nelem 1000 N 6 Nepad 5 KernelV $currVersion STDKernelV $VerSD &>>$OutFile 
	./main CUDA Nelem 800 N 7 Nepad 1 KernelV $currVersion STDKernelV $VerSD &>>$OutFile 
	./main CUDA Nelem 800 N 8 Nepad 3 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
	./main CUDA Nelem 700 N 9 Nepad 4 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
	./main CUDA Nelem 650 N 10 Nepad 1 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
	./main CUDA Nelem 600 N 11 Nepad 2 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
	./main CUDA Nelem 550 N 12 Nepad 2 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
	./main CUDA Nelem 550 N 13 Nepad 2 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
	./main CUDA Nelem 500 N 14 Nepad 2 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
	./main CUDA Nelem 475 N 15 Nepad 2 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
}
