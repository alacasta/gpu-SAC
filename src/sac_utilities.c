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

#include "sac_utilities.h"
#ifndef _WIN32
int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p)
{
  return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
           ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}
#endif
double subinterpolate(double x0, double y0, double x1, double y1, double valor){
	if((x1-x0)==0) return x1;
	else return y0+(y1-y0)*(valor-x0)/(x1-x0);
}

int buscar_double(double valor,int n, double *x){
//busca entre qué dos valores está un número dado y devuelve el índice del
//menor de esos valores
//
	int i;
	int pos;

	pos=0;

	if(valor>=x[n-1]){
		pos=n-1;
	}else{
		for (i=0;i<n-1;i++){
			if((fabs(valor-x[i])<TOL12) || (valor>x[i])&&(valor<x[i+1])){ //el igual en precision del ordenador lo ponemos con TOL12
				pos=i;
			}
	
		}
		if(fabs(valor-x[n-2])<TOL12){ //falta evaluar con el penúltimno para evitar problemas de precisión
			pos=n-2;
		}
	}
	return pos;

}

double interpolate(double valor,int n, double *x, double *y){
	int index;
	index=buscar_double(valor,n,x);
	return(subinterpolate(x[index],x[index+1],y[index],y[index+1],valor));
}

void print_welcome(){
	printf("-----------------------------------------\n");
	printf("GPU SAC\nGPU Accelerated Sacramento Soil Mosture Accounting Model\n");
	printf("Asier Lacasta. 2015\n");
	printf("-----------------------------------------\n");
}

void print_init(sac_stream model, double tmax){
	printf("Initial Date: %d/%d\n",model.prop.month0,model.prop.day0);
	printf("End Date: %d/%d\n",model.prop.monthN,model.prop.dayN);
	printf("Simulation will be performed for t=%.1f hours with time-step: %d\n",tmax,6);
	printf("------------------------------------------------------------------------------\n");

}
void assign_frozen_ground(sac_stream *model, int state){
	 model->state0.uz_twf=model->records[state].uz_twf;
     model->state0.uz_fwf=model->records[state].uz_fwf;
	 model->state0.lz_twf=model->records[state].lz_twf;
	 model->state0.lz_fsf=model->records[state].lz_fsf;
	 model->state0.lz_fpf=model->records[state].lz_fpf;
}

int sac_init_grid(sac_grid *grid, int ntimes, int nx, int ny, double dx){

	int i,j;
	grid->nx=nx;
	grid->ny=ny;
	grid->nCells=nx*ny;
	grid->deltax=dx;
	grid->ntimes=ntimes;
	grid->prec=(double*) malloc(sizeof(double)*nx*ny*ntimes);
	grid->grnd=(double*) malloc(sizeof(double)*nx*ny);
	grid->surf=(double*) malloc(sizeof(double)*nx*ny);
	grid->tet=(double*) malloc(sizeof(double)*nx*ny);

	grid->uz_twc=(double*) malloc(sizeof(double)*nx*ny);
	grid->uz_fwc=(double*) malloc(sizeof(double)*nx*ny);
	grid->lz_twc=(double*) malloc(sizeof(double)*nx*ny);
	grid->lz_fsc=(double*) malloc(sizeof(double)*nx*ny);
	grid->lz_fpc=(double*) malloc(sizeof(double)*nx*ny);
	grid->adimc=(double*) malloc(sizeof(double)*nx*ny);
	grid->uz_twc0=(double*) malloc(sizeof(double)*nx*ny);
	grid->uz_fwc0=(double*) malloc(sizeof(double)*nx*ny);
	grid->lz_twc0=(double*) malloc(sizeof(double)*nx*ny);
	grid->lz_fsc0=(double*) malloc(sizeof(double)*nx*ny);
	grid->lz_fpc0=(double*) malloc(sizeof(double)*nx*ny);
	grid->adimc0=(double*) malloc(sizeof(double)*nx*ny);

	for(i=0;i<nx*ny;i++){
		grid->surf[i]=0.0;
	}
#if SILENTMODE==0
	printf("There are %d MB allocated\n",sizeof(double)*nx*ny*12/(1024*1024)+sizeof(double)*nx*ny*ntimes/(1024*1024));
#endif

	return(0);

}

int sac_generate_random_prec(sac_grid *grid, int ntimes, double min, double max){

	int i,j;
	double space;

	space=max-min;

	for(i=0;i<ntimes;i++){
		for(j=0;j<grid->nCells;j++){
			grid->prec[i*(grid->nCells)+j]=min+(double)rand()/(double)(RAND_MAX)*(space);
		}
	}

	return(1);
}

int sac_generate_sinusoidal_prec(sac_grid *grid, int ntimes, double min, double max){

	int i,j,k;
	double space;
	double amplitude;

	space=max-min;
	


	for(i=0;i<ntimes;i++){
		for(j=0;j<grid->ny;j++){
			amplitude=(sin(3.141598*((double)i/(double)(ntimes-1))+3.141598*((double)j/(double)(grid->ny-1)))+1)*space;
			for(k=0;k<grid->nx;k++){
				grid->prec[i*(grid->nCells)+j*grid->nx+k]=min+amplitude;
			}
		}
	}

	return(1);
}

int update_hours(int *h, int *d, int *m, int *y, int nhours){
	if(*h+nhours>=24){
		*h=(*h+nhours)%24;
		if(*m==1||*m==3||*m==5||*m==7||*m==8||*m==10||*m==12){ //31
			*d++;
			if(*d==32){
				*d=1;
				if(*m==12){
					*m=1;
					*y++;
				}else{
					*m=*m+1;
				}
			}
		}else{
			if(*m==2){ // 28
				*d=*d+1;
				if(*d==29){
					*d=1;
					if(*m==12){
						*m=1;
						*y=*y+1;
					}else{
						*m=*m+1;
					}
				}
			}else{ //30
				*d=*d+1;
				if(*d==31){
					*d=1;
					if(*m==12){
						*m=1;
						*y=*y+1;
					}else{
						*m=*m+1;
					}
				}
			}
		}
	}else{
		*h+=nhours;
		return(1);
	}
	return(1);
}