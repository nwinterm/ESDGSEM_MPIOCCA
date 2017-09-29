#!/bin/bash
#checking kernels
#nohup nvprof --force-overwrite --log-file ReadableOutput.log --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 1 STDKernelV 1 >Results.txt & disown
#nohup nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 0 STDKernelV 0 >Results1.txt & disown
#nohup nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 100 N 5 Nepad 4 KernelV 1 STDKernelV 1 >Results2.txt & disown
StartVerFD=0
EndVerFD=7
VerSD=7
OutFile=NVprofResults.txt
rm -f $OutFile
for ((currVersion=$StartVerFD;currVersion<=$EndVerFD;currVersion++)){
#echo $currVersion
	nvprof --kernels VolumeKernelFluxDiff0 --metrics flop_sp_efficiency,achieved_occupancy,flop_count_sp,dram_read_throughput,dram_write_throughput ./main CUDA Nelem 400 N 4 Nepad 16 KernelV $currVersion STDKernelV $VerSD &>>$OutFile
}
