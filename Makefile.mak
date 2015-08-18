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

ifeq ($(DEBUG),no)
	CFLAGS=$(CUDA_REQ_FLAGS) $(C_REQ_FLAGS) $(C_OPT_FLAGS)
endif

ifeq ($(WARN),yes)
      CFLAGS=$(C_REQ_FLAGS) $(C_OPT_FLAGS) $(C_WARN_FLAGS)
endif

ifeq ($(DEBUG),yes)
		CFLAGS=$(C_REQ_FLAGS) $(C_DEBUG_FLAGS)$(D_LEVEL)

endif
LDFLAGS=$(C_FLAGS) -lrt 


OBJGPU = sac_utilities.o sac_geometry.o sac_io.o sac_g_engine.o sac_gpu.o
BINGPU = sac_gpu

OBJCPU = sac_utilities.o sac_geometry.o sac_io.o sac_engine.o sac_cpu.o
BINCPU = sac_cpu

## $@ Es el token de la izquierda de los :
## $< Es el token de la derecha de los :
ifeq ($(GPUCOMPUTING),1)
$(BINGPU): $(OBJGPU)
	$(NVCC)  $(CFLAGS) $(LDFLAGS) -o $@  $(OBJGPU) -lm
endif
ifeq ($(GPUCOMPUTING),0)
$(BINGPU): $(OBJGPU)
	$(CC) -DGPU=0 $(CFLAGS) $(LDFLAGS) -o $@  $(OBJCPU) -lm
endif
#$(BINCPU): $(OBJCPU)
#	$(CC) $(CFLAGS) $(LDFLAGS) -o $@  $(OBJCPU) -lm
	
clean:
	$(RM) $(OBJ) $(BINCPU) $(BINCPU) $(OBJCPU) $(OBJGPU)

cleanEx:
	$(RM) $(OBJ) $(BIN) $(TEST)

sac_gpu.o: SAC-GPU.c
	$(CC) $(CFLAGS) $< -c -o $@
sac_utilities.o: sac_utilities.c
	$(CC) $(CFLAGS) $< -c -o $@
sac_geometry.o: sac_geometry.c
	$(CC) $(CFLAGS) $< -c -o $@
sac_io.o: sac_io.c
	$(CC) $(CFLAGS) $< -c -o $@
sac_engine.o: sac_engine.c
	$(CC) $(CFLAGS) $< -c -o $@
sac_g_engine.o: sac_g_engine.cu
	$(NVCC) $(CULIBFLAGS) $< -c -o $@

