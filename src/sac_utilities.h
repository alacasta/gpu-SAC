/*--------------------------------------------------------
GPU-SAC 
GPU Accelerated Sacramento Soil Mosture Accounting Model
----------------------------------------------------------
sac_utilities library
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xmlreader.h>
#include <libxml/xpath.h>*/
#include "omp.h"
#include "math.h"
#include "time.h"
#ifdef _WIN32
#include <windows.h>
#endif
#include "sac_define.h"
#include "sac_structs.h"

#ifndef _WIN32
int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p);
#endif
void print_welcome();
double subinterpolate(double x0, double y0, double x1, double y1, double valor);
int buscar_double(double valor,int n, double *x);
double interpolate(double valor,int n, double *x, double *y);
void print_init(sac_stream model, double tmax);
void assign_frozen_ground(sac_stream *model, int state);
int sac_generate_random_prec(sac_grid *grid, int ntimes, double min, double max);
int sac_init_grid(sac_grid *grid, int ntimes, int nx, int ny, double dx);
int sac_generate_sinusoidal_prec(sac_grid *grid, int ntimes, double min, double max);
int update_hours(int *h, int *d, int *m, int *y, int nhours);