#!/bin/bash
./sacgpu benchmark/ 8 >> sac_gpu_times.out
./sacgpu benchmark/ 16 >> sac_gpu_times.out
./sacgpu benchmark/ 32 >> sac_gpu_times.out
./sacgpu benchmark/ 64 >> sac_gpu_times.out
./sacgpu benchmark/ 128 >> sac_gpu_times.out
./sacgpu benchmark/ 256 >> sac_gpu_times.out
./sacgpu benchmark/ 512 >> sac_gpu_times.out
./saccpu benchmark/ 8 >> sac_cpu_times.out
./saccpu benchmark/ 16 >> sac_cpu_times.out
./saccpu benchmark/ 32 >> sac_cpu_times.out
./saccpu benchmark/ 64 >> sac_cpu_times.out
./saccpu benchmark/ 128 >> sac_cpu_times.out
./saccpu benchmark/ 256 >> sac_cpu_times.out
./saccpu benchmark/ 512 >> sac_cpu_times.out
