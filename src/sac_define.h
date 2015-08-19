/*--------------------------------------------------------
GPU-SAC 
GPU Accelerated Sacramento Soil Mosture Accounting Model
----------------------------------------------------------
Definitions for the GPU-SAC model
----------------------------------------------------------
The MIT License (MIT)

Copyright (c) 2015 Asier Lacasta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
---------------------------------------------------------*/

#define _TDATA double
#define RANDOMRAIN 1
#define SILENTMODE 1
#define VTK 0
#define _DUMP_OUTPUT 0
#define nThreads 128
/**< Tolerances */

#define TOL3 1e-3
#define TOL4 1e-4
#define TOL5 1e-5
#define TOL6 1e-6
#define TOL8 1e-8
#define TOL10 1e-10
#define TOL11 1e-11
#define TOL12 1e-12
#define TOL13 1e-13
#define TOL14 1e-14
#define TOL15 1e-15
#ifndef MSGINFO
#define MSGINFO " \033[1;30m[##]\033[0m "
#endif
#ifndef MSGOK
#define MSGOK " \033[1;35m[OK]\033[0m "
#endif
#ifndef MSG1
#define MSG1 " \033[1;36m[_]\033[0m "
#endif
#ifndef MSG2
#define MSG2 " \033[1;37m[_]\033[0m "
#endif
#define CUT_CHECK_ERROR(errorMessage) do {                                 \
        cudaThreadSynchronize();                                           \
         cudaError_t err = cudaGetLastError();                             \
         if( cudaSuccess != err) {                                         \
                     fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                                             errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
                     exit(EXIT_FAILURE);                                                  \
                 } } while (0)
