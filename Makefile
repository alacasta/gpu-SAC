#################################################################
## Define wether uning CPU (GPU=0) or GPU (GPU=1)
GPUCOMPUTING=1
#################################################################

# Define C compiler
CC = gcc
NVCC = nvcc
# Define C flags
C_PROFILE_FLAGS = -pg
C_REQ_FLAGS = -D GPU=$(GPUCOMPUTING)
ifeq ($(GPUCOMPUTING),1)
CUDA_REQ_FLAGS =  -O3 -w -I/usr/local/cuda-6.0/include -L/usr/local/cuda-6.0/lib64 -lcuda -lcudart -lcublas -lrt -lm
CULIBFLAGS = -O3 --compiler-bindir=gcc-4.8 -w -m64 -arch sm_35 -Xptxas -dlcm=ca --use_fast_math  -lcuda -lcudart -lcublas -lrt -lm 
endif
ifeq ($(GPUCOMPUTING),0)
CUDA_REQ_FLAGS = -O3 -w -lrt -lm
CULIBFLAGS = -O3 -lrt -lm 
endif
C_DEBUG_FLAGS = --debug-device
C_MUD_FLAGS = -fmudflap -lmudflap
C_OPT_FLAGS =
C_WARN_FLAGS = -Wall
CFLAGS=$(CUDA_REQ_FLAGS) $(C_REQ_FLAGS) $(C_OPT_FLAGS)

ifeq ($(WARN),yes)
      CFLAGS=$(C_REQ_FLAGS) $(C_OPT_FLAGS) $(C_WARN_FLAGS)
endif

ifeq ($(DEBUG),yes)
		CFLAGS=$(C_REQ_FLAGS) $(C_DEBUG_FLAGS)$(D_LEVEL)
endif
LDFLAGS=$(C_FLAGS) -lrt 


OBJGPU = src/sac_utilities.o src/sac_io.o src/sac_g_engine.o GPU-SAC_gpu.o
BINGPU = GPU-SAC_gpu

OBJCPU = src/sac_utilities.o src/sac_io.o src/sac_engine.o GPU-SAC_cpu.o
BINCPU = GPU-SAC_cpu

ifeq ($(GPUCOMPUTING),1)
BINCLEAN = GPU-SAC_gpu
endif
ifeq ($(GPUCOMPUTING),0)
BINCLEAN = GPU-SAC_cpu
endif
## $@ Es el token de la izquierda de los :
## $< Es el token de la derecha de los :
ifeq ($(GPUCOMPUTING),1)
$(BINGPU): $(OBJGPU)
	$(NVCC) -D GPU=1 $(CFLAGS) $(LDFLAGS) -o $@  $(OBJGPU) -lm
endif
ifeq ($(GPUCOMPUTING),0)
$(BINCPU): $(OBJCPU)
	$(CC) -D GPU=0 $(CFLAGS) $(LDFLAGS) -o $@  $(OBJCPU) -lm
endif
#$(BINCPU): $(OBJCPU)
#	$(CC) $(CFLAGS) $(LDFLAGS) -o $@  $(OBJCPU) -lm
	
clean:
	$(RM) $(OBJ) $(BINCLEAN) $(OBJCPU) $(OBJGPU)

cleanEx:
	$(RM) $(OBJ) $(BIN) $(TEST)

	
GPU-SAC_cpu.o: GPU-SAC.c
	$(CC) $(CFLAGS) $< -c -o $@
	
GPU-SAC_gpu.o: GPU-SAC.c
	$(CC) $(CFLAGS) $< -c -o $@
	
src/sac_utilities.o: src/sac_utilities.c
	$(CC) $(CFLAGS) $< -c -o $@

src/sac_io.o: src/sac_io.c
	$(CC) $(CFLAGS) $< -c -o $@
	
src/sac_engine.o: src/sac_engine.c
	$(CC) $(CFLAGS) $< -c -o $@
	
src/sac_g_engine.o: src/sac_g_engine.cu
	$(NVCC) $(CULIBFLAGS) $< -c -o $@

