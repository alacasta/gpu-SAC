#!/bin/bash
./GPU-SAC_gpu benchmark/ 8 >> sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 16 >> sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 32 >> sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 64 >> sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 128 >> sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 256 >> sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 512 >> sac_gpu_times.out
./GPU-SAC_cpu benchmark/ 8 >> sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 16 >> sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 32 >> sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 64 >> sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 128 >> sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 256 >> sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 512 >> sac_cpu_times.out
