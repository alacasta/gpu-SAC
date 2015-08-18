/*--------------------------------------------------------
GPU-SAC 
GPU Accelerated Sacramento Soil Mosture Accounting Model
----------------------------------------------------------
Input/Output library
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

#include "sac_structs.h"
#include "sac_utilities.h"
int sac_read_config(sac_grid *grid, char *dir, char *filename);
int sac_read_read_dimension(char *dir, char *pref, int hour,int day, int month, int year, int *nx, int *ny);
int sac_read_stream_props(sac_stream_properties *props, char *dir,char *filename);
int sac_read_stream_props_v2(sac_stream_properties *props, char *dir,char *filename);
int sac_read_stream_state(sac_stream_state *state, char *dir,char *filename);
int sac_read_stream_state_v2(sac_stream *stream, char *dir,char *filename);
int save_state(sac_stream *stream, char *dir,char *filename, int step, double time);
int save_state_vtk(sac_grid *grid, char *dir,int step, double time);
int sac_read_grid_stream_state(sac_grid *stream, char *dir,char *filename);
int save_state_vtk_surf(sac_grid *grid, char *dir,int step, double time);
int save_state_surf(sac_grid *grid, char *dir,int step, double time);
int save_state_prec(sac_grid *grid, char *dir,int step, double time);
int sac_read_grid_rain(sac_grid *grid, char *dir, char *pref, int hour,int day, int month, int year,int ntimes);
int save_boundary_vtk(sac_grid *grid, char *dir,  char *pref, int hour,int day, int month, int year);
int save_boundary(sac_grid *grid, char *dir,  char *pref, int hour,int day, int month, int year);
