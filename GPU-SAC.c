/*--------------------------------------------------------
GPU-SAC 
GPU Accelerated Sacramento Soil Mosture Accounting Model
----------------------------------------------------------
MAIN function
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
#pragma comment(lib, "cuda.lib")
#pragma comment(lib, "cudart.lib")

#include <stdio.h>

#include "sac_define.h"
#include "sac_utilities.h"
#include "sac_geometry.h"
#include "sac_engine.h"
#include "time.h"
#include "stdint.h"

#include "sac_io.h"
#if GPU==1
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"
extern void sac_g_init_grid(sac_grid *grid, int ntimes, int nx, int ny, double dx);
extern void sac_g_props(sac_grid *gridCPU, sac_stream_properties* propsGPU);
extern void sac_g_cpu_2_gpu(sac_grid *gridCPU, sac_grid* gridGPU);
extern void sac_do_g_grid_timestep(sac_grid *stream, double t, double deltat,
					int index,		 //  Index of the precipitation
					double pedemand,	 //  evapotranspiration-demand
					int ncells,
					sac_stream_properties *props
					);
#endif

#if GPU==1
#include "sac_g_engine.cuh"
#endif

int main(int argc, char *argv[])
{
	char directory[1024];
#if ONESTREAM==1
	sac_stream model;
#else
	sac_grid model;
	sac_grid g_model;
	sac_stream_properties* g_props;
    //cudaDeviceProp properties;
#endif
	double t,tmax;
	double prec,et;
	int currentStep;
	double dt,dthours;
	int ntimes,nx,ny;
	double dx;
	double min,max;
	int month;
	int num_devices,device;
	#ifdef _WIN32
	  DWORD dwStartTime;
	  DWORD dwElapsed;
	#else
	struct timespec start, end;
  	time_t rawtime;
  	struct tm * timeinfo;
	uint64_t timeElapsed;
	#endif

	if(argc>1){
#ifdef _WIN32
		sprintf(directory,"%s\\",argv[1]);
#else
		sprintf(directory,"%s",argv[1]);
#endif
	}

	if(argc>2){
		nx=atoi(argv[2]);
		ny=nx;
	}else{
		nx=512;
		ny=390;
		nx=256;
		ny=256;
		nx=1024;
		ny=1024;

	}
	#if RANDOMRAIN==0
	sac_read_read_dimension(directory,"xmrg",00,1,06,2014,&nx,&ny);

	#endif
#if SILENTMODE==0
	print_welcome();
	printf("The resoultion of the grid is %dx%d\n",nx,ny);
#endif		
	sac_read_stream_props_v2(&model.prop,directory,"SACSMA_BKRM5_BKRM5_UpdateStates.txt");


#if GPU==1
	cudaSetDevice(0);
#endif



#if ONESTREAM
	sac_read_stream_state_v2(&model,directory,"bkrm5_sacsma.csv");
	tmax=model.nrecords*model.prop.fr_out_interval;
	tmax=15*model.prop.fr_out_interval;

	t=0.0;
	
	dt=0.25; // In days!!!!!
	dthours= dt*24.0;

	print_init(model,tmax);

    printf("uz_twc\tuz_fwc\tadimc\tlz_twc\tlz_fsc\tlz_fpc\t\t Surf\t  Grnd\t  TET\n" );

	currentStep=0;
	while(t<tmax){
		prec=model.prec.e[currentStep+1];
		et=model.prop.et_demand_curve[model.records[currentStep].month-1];
		//et/=(24.*60./6.0);
		et/=(24./dthours);
		if(prec<0.0){
			prec=0.0;
		}
		assign_frozen_ground(&model,currentStep);
		sac_timestep(&model,t,dt,prec,et);
		t+=dthours;

		save_state(&model,directory,"bkrm5_sacsma.out",currentStep,t);
		/*printf("t=%lf surf: %.6e ground %.6e ET: %.6e\n",
			t,
			model.state.surf,
			model.state.grnd,
			model.state.tet);*/
		// Move current state as state0 for the next time-step
		model.state0=model.state;
		currentStep++;
	}
#else
	ntimes=256;
	sac_read_config(&model,directory,"config.txt");
	ntimes=model.final_t/model.delta_t;
	tmax=model.final_t;
	//tmax=ntimes*model.prop.fr_out_interval;
	t=0.0;
	dt=model.delta_t/24.0; // In days!!!!!
	dthours= dt*24.0;	

	dx=5.0;
	month=5;

	min=0.0;
	max=24.4;
	
	sac_init_grid(&model,ntimes,nx,ny,dx);
	//sac_generate_random_prec(&model, ntimes,min,max);
	#if RANDOMRAIN==0
		ntimes=sac_read_grid_rain(&model,directory,"xmrg",00,1,06,2014,ntimes);
		model.final_t=ntimes*model.delta_t;
		tmax=model.final_t;
		model.ntimes=ntimes;
		//save_boundary_vtk(&model,directory,"xmrg",00,1,06,2014);
		save_boundary(&model,directory,"xmrg",00,1,06,2014);

	#else
		sac_generate_sinusoidal_prec(&model, ntimes,min,max);
	#endif
	sac_read_grid_stream_state(&model,directory,"initialState.txt");

#if GPU==1
	sac_g_init_grid(&g_model,ntimes,nx,ny,dx);
    //sac_g_props(&model, g_props);
	cudaMalloc((void**)&g_props,sizeof(sac_stream_properties));
	sac_g_props(&model, g_props);
	g_model.g_prop=g_props;
    sac_g_cpu_2_gpu(&model, &g_model);
#endif

	//save_state_vtk(&model,directory,0,0.0);
	//tmax=model.final_t;
	//tmax=ntimes*model.prop.fr_out_interval;
	//t=0.0;
	//dt=0.25; // In days!!!!!
	//dthours= dt*24.0;
	#ifdef _WIN32
	dwStartTime = GetTickCount()	;
	#else
	clock_gettime(CLOCK_MONOTONIC, &start);
  	time (&rawtime);
 	timeinfo = localtime (&rawtime);
#if SILENTMODE==0
	printf ("%s Simulation starts at %s",MSGINFO, asctime(timeinfo));	
#endif
	#endif;
	currentStep=0;
	while(t<tmax){
		et=model.prop.et_demand_curve[month-1];
		et/=(24./dthours);
#if GPU==1
		sac_do_g_grid_timestep(
			&g_model,
			t,
			dt,
			currentStep,
			et,
			model.nCells,
			g_props);
#else
		sac_grid_timestep(&model,t,dt,currentStep,et);
#endif

#if _DUMP_OUTPUT==1
		if((int)t%(int)model.delta_t_dump==0){
			#if GPU==1
				sac_g_gpu_2_cpu(&model, &g_model);
			#endif
			save_state_vtk(&model,directory,currentStep,0.0);
			printf("%s Dumping: %d\n",MSG2,currentStep);

		}
#else
		if((int)t%(int)model.delta_t_dump==0){
			#if GPU==1
				sac_g_gpu_2_cpu(&model, &g_model);
			#endif
#if VTK==1
			save_state_vtk_surf(&model,directory,currentStep,0.0);
#else
			save_state_surf(&model,directory,currentStep,0.0);
			save_state_prec(&model,directory,currentStep,0.0);

#endif
			//printf("%s Dumping: %d\n",MSG2,currentStep);

		}
#endif
		t+=dthours;
#if SILENTMODE==0
		printf("%s t: %lf\n",MSG1,t);
#endif
		currentStep++;
	}
		#ifndef _WIN32
		clock_gettime(CLOCK_MONOTONIC, &end);
		timeElapsed = timespecDiff(&end, &start);

  		time (&rawtime);
 		timeinfo = localtime (&rawtime);
#if SILENTMODE==0
		printf("*************************************\n");
		printf("%s Time elapsed: %lf seconds\n",MSGOK,timeElapsed/1000000000.);
		printf("*************************************\n");
  		printf ("%s Simulation ends at %s",MSGINFO, asctime(timeinfo));
#else
		printf("%d %lf\n",model.nx,timeElapsed/1000000000.);
#endif
		#else
		dwElapsed = GetTickCount() - dwStartTime;
#if SILENTMODE==0
		printf("*************************************\n");
		printf("Time elapsed: %d.%3d seconds\n",dwElapsed/1000, dwElapsed - dwElapsed/1000);
		printf("*************************************\n");
#else
		printf("%d  %d.%d\n",model.nx,dwElapsed/1000, dwElapsed - dwElapsed/1000);
#endif
		#endif
		#if GPU==1
		sac_g_gpu_2_cpu(&model, &g_model);
		#endif
#if VTK==1
			save_state_vtk_surf(&model,directory,currentStep-1,0.0);
#else
			save_state_surf(&model,directory,currentStep-1,0.0);
			save_state_prec(&model,directory,currentStep-1,0.0);

#endif

#endif



  return 0;
}