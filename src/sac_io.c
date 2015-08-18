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

#include "sac_io.h"

int sac_read_config(sac_grid *grid, char *dir, char *filename){
	FILE *f;
	char ftotal[1024];
	int i;
	sprintf(ftotal,"%s%s",dir,filename);

	f=fopen(ftotal,"r");
	
	if(f==NULL){
		printf("Propierties file %s was not found\n",ftotal);
		return(-1);
	}
	
	fscanf(f,"%*s %lf",&grid->final_t);
	fscanf(f,"%*s %lf",&grid->delta_t);
	fscanf(f,"%*s %lf",&grid->delta_t_dump);
	
	return(0);
}

int sac_read_read_dimension(char *dir, char *pref, int hour,int day, int month, int year, int *nx, int *ny){
	FILE *f;
	char ftotal[1024];
	int i,a,b;


	sprintf(ftotal,"%s%s%02d%02d%4d%02dz.ascii",dir,pref,month,day,year,hour);
	f=fopen(ftotal,"r");
	if(f==NULL){
		printf("sac_read_read_dimension file %s was not found\n",ftotal);
		return(-1);
	}
		
	fscanf(f,"%*s %d",&a);
	fscanf(f,"%*s %d",&b);
	*nx=a;
	*ny=b;

	fclose(f);

	return(0);

}

int sac_read_grid_rain(sac_grid *grid, char *dir, char *pref, int hour,int day, int month, int year,int ntimes){
	FILE *f;
	char ftotal[1024];
	int i,j,k;
	int total;
	int d,h,m,y;
	d=day;
	m=month;
	y=year;
	h=hour;
	total=0;
	free(grid->prec);
	grid->prec =(double*)malloc(sizeof(double)*grid->nCells*ntimes);

	for(k=0;k<ntimes;k++){
		sprintf(ftotal,"%s%s%02d%02d%4d%02dz.ascii",dir,pref,m,d,y,h);
		f=fopen(ftotal,"r");
		if(f==NULL){
			printf("sac_read_grid_rain file %s was not found\n",ftotal);
			return(-1);
		}
		fscanf(f,"%*s %*d");
		fscanf(f,"%*s %*d");
		fscanf(f,"%*s %*lf");
		fscanf(f,"%*s %*lf");
		fscanf(f,"%*s %lf",&grid->deltax);
		fscanf(f,"%*s %*lf");
		for(j=grid->ny-1;j>=0;j--){
			for(i=0;i<grid->nx;i++){
				fscanf(f,"%lf",&grid->prec[k*grid->nCells+j*grid->nx+i]);
				if(grid->prec[k*grid->nCells+j*grid->nx+i]<0.0){
					grid->prec[k*grid->nCells+j*grid->nx+i]=0.0;
				}
			}
		}
		update_hours(&h,&d,&m,&y,1);
		fclose(f);
		total++;
	}
	return(total);

}

int sac_read_stream_props(sac_stream_properties *props, char *dir,char *filename){
	FILE *f;
	char ftotal[1024];
	int i;
	sprintf(ftotal,"%s%s",dir,filename);

	f=fopen(ftotal,"r");
	
	if(f==NULL){
		printf("Propierties file %s was not found\n",ftotal);
		return(-1);
	}

	fscanf(f,"%*s %lf",&props->r_exp);
	fscanf(f,"%*s %d",&props->we_input_option);
	fscanf(f,"%*s %lf",&props->lzpk);
	fscanf(f,"%*s %lf",&props->uzk);
	fscanf(f,"%*s %lf",&props->side);
	fscanf(f,"%*s %lf",&props->lz_fpm);
	fscanf(f,"%*s %d",&props->runoff_interval);
	fscanf(f,"%*s %lf",&props->px_adj);
	fscanf(f,"%*s %lf",&props->pct_free);
	fscanf(f,"%*s %lf",&props->lzsk);
	fscanf(f,"%*s %lf",&props->lz_fsm);
	fscanf(f,"%*s %lf",&props->z_perc);
	fscanf(f,"%*s %d",&props->mape_input);
	fscanf(f,"%*s %lf",&props->riva);
	fscanf(f,"%*s %lf",&props->pe_adj);
	fscanf(f,"%*s %lf",&props->uz_twm);
	fscanf(f,"%*s %lf",&props->lz_twm);
	fscanf(f,"%*s %lf",&props->uz_fwm);
	fscanf(f,"%*s %lf",&props->rserv);
	fscanf(f,"%*s %lf",&props->pct_im);
	fscanf(f,"%*s %lf",&props->adimp);
	fscanf(f,"%*s %lf",&props->efc);
	if(props->mape_input==1){
		fscanf(f,"%*s");
		props->et_demand_curve=(double*)malloc(sizeof(double)*12);
		for(i=0;i<12;i++){
			fscanf(f,"%lf",&props->et_demand_curve[i]);
		}
	}
	fscanf(f,"%*s %d",&props->isfrozen);
	return(1);
}
int sac_read_stream_props_v2(sac_stream_properties *props, char *dir,char *filename){
	FILE *f;
	char ftotal[1024];
	int i;
	sprintf(ftotal,"%s%s",dir,filename);

	f=fopen(ftotal,"r");
	
	if(f==NULL){
		printf("Propierties file %s was not found\n",ftotal);
		return(-1);
	}

	fscanf(f,"%*s %lf",&props->csnow);
	fscanf(f,"%*s %lf",&props->r_exp);
	fscanf(f,"%*s %lf",&props->lzpk);
	fscanf(f,"%*s %lf",&props->fr_temp);
	fscanf(f,"%*s %lf",&props->lz_fpm);
	fscanf(f,"%*s %lf",&props->px_adj);
	fscanf(f,"%*s %lf",&props->pct_free);
	fscanf(f,"%*s %lf",&props->z_perc);
	fscanf(f,"%*s %lf",&props->riva);
	fscanf(f,"%*s %lf",&props->pe_adj);
	fscanf(f,"%*s %lf",&props->lz_twm);
	fscanf(f,"%*s %lf",&props->rserv);
	fscanf(f,"%*s %lf",&props->adimp);
	fscanf(f,"%*s %d",&props->isfrozen);
	fscanf(f,"%*s %lf",&props->fr_exp);
	fscanf(f,"%*s %lf",&props->satr);
	fscanf(f,"%*s %d",&props->sasc_input_option);
	fscanf(f,"%*s %lf",&props->uzk);
	fscanf(f,"%*s %lf",&props->side);
	fscanf(f,"%*s %lf",&props->lz_fsm);
	fscanf(f,"%*s %lf",&props->lzsk);
	fscanf(f,"%*s %lf",&props->csoil);
	fscanf(f,"%*s %lf",&props->uz_twm);
	fscanf(f,"%*s %lf",&props->uz_fwm);
	fscanf(f,"%*s %d",&props->we_input_option);
	fscanf(f,"%*s %lf",&props->pct_im);
	fscanf(f,"%*s %lf",&props->rthaw);
	fscanf(f,"%*s %lf",&props->gch);
	fscanf(f,"%*s %lf",&props->efc);
	props->mape_input=1;
	if(props->mape_input==1){
		fscanf(f,"%*s");
		props->et_demand_curve=(double*)malloc(sizeof(double)*12);
		for(i=0;i<12;i++){
			fscanf(f,"%lf",&props->et_demand_curve[i]);
		}
	}
	fscanf(f,"%*s %d",&props->runoff_interval);
	fscanf(f,"%*s %d",&props->fr_out_interval);
	fscanf(f,"%*s %d",&props->sm_interval);
	return(1);
}
int sac_read_stream_state(sac_stream_state *state, char *dir,char *filename){
	FILE *f;
	char ftotal[1024];
	int i;
	sprintf(ftotal,"%s%s",dir,filename);

	f=fopen(ftotal,"r");
	
	if(f==NULL){
		printf("State file %s was not found\n",ftotal);
		return(-1);
	}

	fscanf(f,"%*s %lf",&state->uz_twc);
	fscanf(f,"%*s %lf",&state->uz_fwc);
	fscanf(f,"%*s %lf",&state->lz_twc);
	fscanf(f,"%*s %lf",&state->lz_fsc);
	fscanf(f,"%*s %lf",&state->lz_fpc);
	fscanf(f,"%*s %lf",&state->adimc);
	fclose(f);
	return(1);
}



int sac_read_stream_state_v2(sac_stream *stream, char *dir,char *filename){
	FILE *f;
	sac_stream_state state;
	int r;
	char ftotal[1024];
	char auxc[1024];
	int month,day;
	int i,k;
	double aux;
	sprintf(ftotal,"%s%s",dir,filename);
	k=0;
	stream->nrecords=0;
	f=fopen(ftotal,"r");
	
	if(f==NULL){
		printf("State file %s was not found\n",ftotal);
		return(-1);
	}

	stream->records=(sac_stream_state*)malloc(sizeof(sac_stream_state));
	stream->prec.t=(double*)malloc(sizeof(double));
	stream->prec.e=(double*)malloc(sizeof(double));

	fscanf(f,"%s,",auxc);

	r=1;
	while(r>0){
		for(i=0;i<27;i++){
				if(i==0){
					fscanf(f,"%d/%d/%*d %*d:%*d,",&month,&day);
				}else{
					r=fscanf(f,"%lf,",&aux);
					// States
					//aux*=24.4;
					aux/=0.0393700787;
					switch (i)
					{
					case 2: state.lz_fpc=aux; break;
					case 6: state.lz_fsc=aux; break;
					case 7: state.uz_fwc=aux; break;
					case 12: state.uz_twc=aux; break;
					case 22: state.adimc=aux; break;
					case 23: state.lz_twc=aux; break;
					case 3: state.lz_fsf=aux; break;
					case 5: state.lz_fpf=aux; break;
					case 9: state.uz_fwf=aux; break;
					case 16: state.adimf=aux; break;
					case 20: state.uz_twf=aux; break;
					case 21: state.lz_twf=aux; break;
					case 4: state.sur_ro=aux; break;
					case 8: state.dir_ro=aux; break;
					case 10: state.sup_ro=aux; break;
					case 17: state.imp_ro=aux; break;
					case 24: state.pri_ro=aux; break;
					case 11: state.infw=aux; break;
					default: 	break;
					}
					
					/*
					// Runoffs
					if(i==4){
						state.sur_ro=aux;
					}
					if(i==8){
						state.dir_ro=aux;
					}
					if(i==10){
						state.sup_ro=aux;
					}
					if(i==17){
						state.imp_ro=aux;
					}
					if(i==24){
						state.pri_ro=aux;
					}



					if(i==2){
						state.lz_fpc=aux;
					}
					if(i==6){
						state.lz_fsc=aux;
					}
					if(i==7){
						state.uz_fwc=aux;
					}
					if(i==12){
						state.uz_twc=aux;
					}
					if(i==22){
						state.adimc=aux;
					}
					if(i==23){
						state.lz_twc=aux;
					}

					// Fractions
					if(i==3){
						state.lz_fsf=aux;
					}
					if(i==5){
						state.lz_fpf=aux;
					}
					if(i==9){
						state.uz_fwf=aux;
					}
					if(i==16){
						state.adimf=aux;
					}
					if(i==20){
						state.uz_twf=aux;
					}
					if(i==21){
						state.lz_twf=aux;
					}
					*/



					if(i==13){
						stream->prec.e[k]=aux;
						stream->prec.t[k]=k*stream->prop.sm_interval;
					}
				}
		}
		state.day=day;
		state.month=month;
		stream->records[k]=state;
		stream->nrecords++;
		stream->prec.ntimes=stream->nrecords;
		stream->records=(sac_stream_state*)realloc(stream->records,sizeof(sac_stream_state)*(stream->nrecords+1));
		stream->prec.t=(double*)realloc(stream->prec.t,sizeof(double)*(stream->nrecords+1));
		stream->prec.e=(double*)realloc(stream->prec.e,sizeof(double)*(stream->nrecords+1));

		if(k==0){
			stream->state0=state;
			stream->prop.month0=month;
			stream->prop.day0=month;
		}else{
			stream->prop.monthN=month;
			stream->prop.dayN=month;
		}
		k++;
	}

	printf("There are %d states\n",stream->nrecords);

	fclose(f);
	return(1);
}

int save_state(sac_stream *stream, char *dir,char *filename, int step, double time){
	FILE *f;
	char ftotal[1024];
	int i;
	sprintf(ftotal,"%s%s",dir,filename);	
	if(step==0){
		f=fopen(ftotal,"w");
	}else{
		f=fopen(ftotal,"a+");
	}
	
	if(f==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}

	fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		time,
		stream->state.uz_fwc,
		stream->state.uz_twc,
		stream->state.adimc,
		stream->state.lz_twc,
		stream->state.lz_fsc,
		stream->state.lz_fpc,
		stream->records[step+1].uz_fwc,
		stream->records[step+1].uz_twc,
		stream->records[step+1].adimc,
		stream->records[step+1].lz_twc,
		stream->records[step+1].lz_fsc,
		stream->records[step+1].lz_fpc,
		stream->state.surf,
		stream->state.grnd,
		stream->records[step+1].lz_fpc,
		stream->records[step+1].lz_fpc
		);

	fclose(f);
	return(0);
}

int save_boundary_vtk(sac_grid *grid, char *dir,  char *pref, int hour,int day, int month, int year){
	FILE *fp,*f;
	char ftotal[1024];
	char ftotal2[1024];
	int i,j;
	double x,y,rx,ry;
	int count;
	int nnodos;
	int ncells;
	int nx,ny;
	double deltay;

	double *info;

	info=(double*)malloc(sizeof(double)*grid->nCells);

	sprintf(ftotal,"%s%s%02d%02d%4d%02dz.ascii",dir,pref,month,day,year,hour);	
	f=fopen(ftotal,"r");

	
	if(f==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}

		fscanf(f,"%*s %*d");
		fscanf(f,"%*s %*d");
		fscanf(f,"%*s %*lf");
		fscanf(f,"%*s %*lf");
		fscanf(f,"%*s %lf",&grid->deltax);
		fscanf(f,"%*s %*lf");
		for(j=grid->ny-1;j>=0;j--){
			for(i=0;i<grid->nx;i++){
				fscanf(f,"%lf",&info[j*grid->nx+i]);
				if(info[j*grid->nx+i]<0.0){
					info[j*grid->nx+i]=1.0;
				}else{
					info[j*grid->nx+i]=0.0;
				}
			}
		}
		fclose(f);

#if GPU==0
	sprintf(ftotal2,"sac_boundary.vtk");
#else
	sprintf(ftotal2,"sac_boundary.vtk");
#endif
	sprintf(ftotal,"%s%s",dir,ftotal2);	

	fp=fopen(ftotal,"w");
	
	if(fp==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}

	fprintf(fp,"# vtk DataFile Version 2.0 \n");
	fprintf(fp,"Titulo \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID \n");

	nx=(grid->nx);
	ny=(grid->ny);
	nnodos=(nx+1)*(ny+1);
	ncells=nx*ny;

	deltay=grid->deltax;

	fprintf(fp,"POINTS %d float \n",nnodos);
	for(i=0;i<ny+1;i++){
	for(j=0;j<nx+1;j++){
	fprintf(fp,"%lf %lf %lf\n",j*1.0*deltay,i*1.0*deltay,0.0);
	}
	}
	fprintf(fp,"CELLS %d %d \n", ncells,ncells * 5);
	count=0;
	for(i=0;i<ny;i++){
	for(j=0;j<nx;j++){
	fprintf(fp,"4 %d %d %d %d\n",
							count+(j),
							count+(nx+1)+(j),
							count+(nx+1)+(j)+1,
							count+(j)+1);	
	}
	count+=(nx+1);
	}
	fprintf(fp,"CELL_TYPES %d\n",ncells);
		 
	for(i=0;i<ncells;i++){
			   fprintf(fp,"9 \n");	
	}
	fprintf (fp, "CELL_DATA %d \n", ncells);
		      
	fprintf (fp, "SCALARS boundary float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",info[i*grid->nx+j]);	
		}
	}

	fclose(fp);
	return(0);
}
int save_boundary(sac_grid *grid, char *dir,  char *pref, int hour,int day, int month, int year){
	FILE *fp,*f;
	char ftotal[1024];
	char ftotal2[1024];
	int i,j;
	double x,y,rx,ry;
	int count;
	int nnodos;
	int ncells;
	int nx,ny;
	double deltay;

	double *info;

	info=(double*)malloc(sizeof(double)*grid->nCells);

	sprintf(ftotal,"%s%s%02d%02d%4d%02dz.ascii",dir,pref,month,day,year,hour);	
	f=fopen(ftotal,"r");

	
	if(f==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}

		fscanf(f,"%*s %*d");
		fscanf(f,"%*s %*d");
		fscanf(f,"%*s %*lf");
		fscanf(f,"%*s %*lf");
		fscanf(f,"%*s %lf",&grid->deltax);
		fscanf(f,"%*s %*lf");
		for(j=grid->ny-1;j>=0;j--){
			for(i=0;i<grid->nx;i++){
				fscanf(f,"%lf",&info[j*grid->nx+i]);
				if(info[j*grid->nx+i]<0.0){
					info[j*grid->nx+i]=1.0;
				}else{
					info[j*grid->nx+i]=0.0;
				}
			}
		}
		fclose(f);

#if GPU==0
	sprintf(ftotal2,"sac_boundary.out");
#else
	sprintf(ftotal2,"sac_boundary.out");
#endif
	sprintf(ftotal,"%s%s",dir,ftotal2);	

	fp=fopen(ftotal,"w");
	
	if(fp==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}


	nx=(grid->nx);
	ny=(grid->ny);
	nnodos=(nx+1)*(ny+1);
	ncells=nx*ny;

	deltay=grid->deltax;

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			if(i==0||j==0||i==ny-1||j==nx-1){
				fprintf(fp,"-10.0 ");
			}else{
				if(info[i*grid->nx+j]==1.0&&(
					info[(i+1)*grid->nx+j]==0.0||
					info[(i-1)*grid->nx+j]==0.0||
					info[(i)*grid->nx+j+1]==0.0||
					info[(i)*grid->nx+j-1]==0.0
					)){
				fprintf(fp,"0.0001 ");	
				}else{
				fprintf(fp,"-10.0 ");	
				}
			}
		}
		fprintf(fp,"\n");	
	}

	fclose(fp);
	return(0);
}


int save_state_vtk(sac_grid *grid, char *dir,int step, double time){
	FILE *fp;
	char ftotal[1024];
	char ftotal2[1024];
	int i,j;
	double x,y,rx,ry;
	int count;
	int nnodos;
	int ncells;
	int nx,ny;
	double deltay;
#if GPU==0
	sprintf(ftotal2,"sac_state_%d.vtk",step);
#else
	sprintf(ftotal2,"sac_gpu_state_%d.vtk",step);
#endif
	sprintf(ftotal,"%s%s",dir,ftotal2);	
	if(step==0){
		fp=fopen(ftotal,"w");
	}else{
		fp=fopen(ftotal,"w");
	}
	
	if(fp==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}

	fprintf(fp,"# vtk DataFile Version 2.0 \n");
	fprintf(fp,"Titulo \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID \n");

	nx=(grid->nx);
	ny=(grid->ny);
	nnodos=(nx+1)*(ny+1);
	ncells=nx*ny;

	deltay=grid->deltax;

	fprintf(fp,"POINTS %d float \n",nnodos);
	for(i=0;i<ny+1;i++){
	for(j=0;j<nx+1;j++){
	fprintf(fp,"%lf %lf %lf\n",j*1.0*deltay,i*1.0*deltay,0.0);
	}
	}
	fprintf(fp,"CELLS %d %d \n", ncells,ncells * 5);
	count=0;
	for(i=0;i<ny;i++){
	for(j=0;j<nx;j++){
	fprintf(fp,"4 %d %d %d %d\n",
							count+(j),
							count+(nx+1)+(j),
							count+(nx+1)+(j)+1,
							count+(j)+1);	
	}
	count+=(nx+1);
	}
	fprintf(fp,"CELL_TYPES %d\n",ncells);
		 
	for(i=0;i<ncells;i++){
			   fprintf(fp,"9 \n");	
	}
	fprintf (fp, "CELL_DATA %d \n", ncells);
		      
	fprintf (fp, "SCALARS prec float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->prec[step*ncells+(i*nx+j)]);	
		}
	}

	fprintf (fp, "SCALARS surf float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->surf[(i*nx+j)]);	
		//	printf("%d %d %lf\n",i,j,grid->surf[(i*nx+j)]);
		}
	//printf("\n");
	}
	fprintf (fp, "SCALARS uz_twc float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->uz_twc[(i*nx+j)]);	
		}
	}
	fprintf (fp, "SCALARS uz_fwc float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->uz_fwc[(i*nx+j)]);	
		}
	}
	fprintf (fp, "SCALARS lz_twc float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->lz_twc[(i*nx+j)]);	
		}
	}
	fprintf (fp, "SCALARS lz_fpc float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->lz_fpc[(i*nx+j)]);	
		}
	}
	fprintf (fp, "SCALARS lz_fsc float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->lz_fsc[(i*nx+j)]);	
		}
	}
	fprintf (fp, "SCALARS lz_fsc float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->lz_fsc[(i*nx+j)]);	
		}
	}
	fprintf (fp, "SCALARS adimc float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->adimc[(i*nx+j)]);	
		}
	}
	fclose(fp);
	return(0);
}

int save_state_vtk_surf(sac_grid *grid, char *dir,int step, double time){
	FILE *fp;
	char ftotal[1024];
	char ftotal2[1024];
	int i,j;
	double x,y,rx,ry;
	int count;
	int nnodos;
	int ncells;
	int nx,ny;
	double deltay;
#if GPU==0
	sprintf(ftotal2,"sac_state_%d.vtk",step);
#else
	sprintf(ftotal2,"sac_gpu_state_%d.vtk",step);
#endif
	sprintf(ftotal,"%s%s",dir,ftotal2);	
	if(step==0){
		fp=fopen(ftotal,"w");
	}else{
		fp=fopen(ftotal,"w");
	}
	
	if(fp==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}

	fprintf(fp,"# vtk DataFile Version 2.0 \n");
	fprintf(fp,"Titulo \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID \n");

	nx=(grid->nx);
	ny=(grid->ny);
	nnodos=(nx+1)*(ny+1);
	ncells=nx*ny;

	deltay=grid->deltax;

	fprintf(fp,"POINTS %d float \n",nnodos);
	for(i=0;i<ny+1;i++){
	for(j=0;j<nx+1;j++){
	fprintf(fp,"%lf %lf %lf\n",j*1.0*deltay,i*1.0*deltay,0.0);
	}
	}
	fprintf(fp,"CELLS %d %d \n", ncells,ncells * 5);
	count=0;
	for(i=0;i<ny;i++){
	for(j=0;j<nx;j++){
	fprintf(fp,"4 %d %d %d %d\n",
							count+(j),
							count+(nx+1)+(j),
							count+(nx+1)+(j)+1,
							count+(j)+1);	
	}
	count+=(nx+1);
	}
	fprintf(fp,"CELL_TYPES %d\n",ncells);
		 
	for(i=0;i<ncells;i++){
			   fprintf(fp,"9 \n");	
	}
	fprintf (fp, "CELL_DATA %d \n", ncells);
		      
	fprintf (fp, "SCALARS prec float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->prec[step*ncells+(i*nx+j)]);	
		}
	}

	fprintf (fp, "SCALARS surf float \n");
	fprintf (fp, "LOOKUP_TABLE default \n");

	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf\n",grid->surf[(i*nx+j)]);	
		//	printf("%d %d %lf\n",i,j,grid->surf[(i*nx+j)]);
		}
	//printf("\n");
	}
	fclose(fp);
	return(0);
}
int save_state_surf(sac_grid *grid, char *dir,int step, double time){
	FILE *fp;
	char ftotal[1024];
	char ftotal2[1024];
	int i,j;
	double x,y,rx,ry;
	int count;
	int nnodos;
	int ncells;
	int nx,ny;
	double deltay;
#if GPU==0
	sprintf(ftotal2,"surf_d.out",step);
#else
	sprintf(ftotal2,"surf_gpu_%d.out",step);
#endif
	sprintf(ftotal,"%s%s",dir,ftotal2);	
	if(step==0){
		fp=fopen(ftotal,"w");
	}else{
		fp=fopen(ftotal,"w");
	}
	
	if(fp==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}


	nx=(grid->nx);
	ny=(grid->ny);
	ncells=nx*ny;

	deltay=grid->deltax;


	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf ",grid->surf[(i*nx+j)]);	
		}
		fprintf(fp,"%\n");	
	}
	fclose(fp);
	return(0);
}
int save_state_prec(sac_grid *grid, char *dir,int step, double time){
	FILE *fp;
	char ftotal[1024];
	char ftotal2[1024];
	int i,j;
	double x,y,rx,ry;
	int count;
	int nnodos;
	int ncells;
	int nx,ny;
	double deltay;
#if GPU==0
	sprintf(ftotal2,"prec_d.out",step);
#else
	sprintf(ftotal2,"prec_gpu_%d.out",step);
#endif
	sprintf(ftotal,"%s%s",dir,ftotal2);	
	if(step==0){
		fp=fopen(ftotal,"w");
	}else{
		fp=fopen(ftotal,"w");
	}
	
	if(fp==NULL){
		printf("Dump file %s was not found\n",ftotal);
		return(-1);
	}


	nx=(grid->nx);
	ny=(grid->ny);
	ncells=nx*ny;

	deltay=grid->deltax;


	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(fp,"%lf ",grid->prec[step*ncells+(i*nx+j)]);	
		}
		fprintf(fp,"%\n");	
	}
	fclose(fp);
	return(0);
}
int sac_read_grid_stream_state(sac_grid *stream, char *dir,char *filename){
	FILE *f;
	sac_stream_state state;
	int r;
	char ftotal[1024];
	char auxc[1024];
	int month,day;
	int i,k;
	double aux;
	sprintf(ftotal,"%s%s",dir,filename);
	k=0;
	f=fopen(ftotal,"r");
	
	if(f==NULL){
		printf("State file %s was not found\n",ftotal);
		return(-1);
	}

	r=1;
	r=fscanf(f,"%*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf",
		&stream->uz_twc[0],
		&stream->uz_fwc[0],
		&stream->adimc[0],
		&stream->lz_twc[0],
		&stream->lz_fpc[0],
		&stream->lz_fsc[0]
	);		

	for(i=1;i<stream->nCells;i++){
		stream->uz_twc[i]=stream->uz_twc[0];
		stream->uz_fwc[i]=stream->uz_fwc[0];
		stream->adimc[i]=stream->adimc[0];
		stream->lz_twc[i]=stream->lz_twc[0];
		stream->lz_fpc[i]=stream->lz_fpc[0];
		stream->lz_fsc[i]=stream->lz_fsc[0];
	}

	for(i=0;i<stream->nCells;i++){
		stream->uz_twc0[i]=stream->uz_twc[0];
		stream->uz_fwc0[i]=stream->uz_fwc[0];
		stream->adimc0[i]=stream->adimc[0];
		stream->lz_twc0[i]=stream->lz_twc[0];
		stream->lz_fpc0[i]=stream->lz_fpc[0];
		stream->lz_fsc0[i]=stream->lz_fsc[0];
	}
	fclose(f);
	return(1);
}
