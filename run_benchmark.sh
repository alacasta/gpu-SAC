#!/bin/bash
echo 'Running GPU tests'
./GPU-SAC_gpu benchmark/ 8 >> benchmark/sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 16 >> benchmark/sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 32 >> benchmark/sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 64 >> benchmark/sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 128 >> benchmark/sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 256 >> benchmark/sac_gpu_times.out
./GPU-SAC_gpu benchmark/ 512 >> benchmark/sac_gpu_times.out
echo 'Running CPU tests'
./GPU-SAC_cpu benchmark/ 8 >> benchmark/sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 16 >> benchmark/sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 32 >> benchmark/sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 64 >> benchmark/sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 128 >> benchmark/sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 256 >> benchmark/sac_cpu_times.out
./GPU-SAC_cpu benchmark/ 512 >> benchmark/sac_cpu_times.out
