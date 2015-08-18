/*--------------------------------------------------------
GPU-SAC 
GPU Accelerated Sacramento Soil Mosture Accounting Model
----------------------------------------------------------
Sacramento GPU model numerical engine
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
#include "sac_g_engine.cuh"

extern "C"  void sac_do_g_grid_timestep(sac_grid *stream, double t, double deltat,
					int index,		 //  Index of the precipitation
					double pedemand,	 //  evapotranspiration-demand
					int ncells,
					sac_stream_properties *props
					){
		sac_g_grid_timestep<<<ncells/nThreads+1,nThreads>>>(
			*stream,
			t,
			deltat,
			index,
			pedemand,
			ncells,
			props
			);
}


extern "C" void sac_g_cpu_2_gpu(sac_grid *gridCPU, sac_grid* gridGPU){

		cudaMemcpy(gridGPU->prec,gridCPU->prec, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny*gridCPU->ntimes, cudaMemcpyHostToDevice );

		cudaMemcpy(gridGPU->grnd,gridCPU->grnd, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->surf,gridCPU->surf, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->tet,gridCPU->tet, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );

		cudaMemcpy(gridGPU->uz_twc,gridCPU->uz_twc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->uz_fwc,gridCPU->uz_fwc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->lz_twc,gridCPU->lz_twc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->lz_fsc,gridCPU->lz_fsc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->lz_fpc,gridCPU->lz_fpc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->adimc,gridCPU->adimc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );

		cudaMemcpy(gridGPU->uz_twc0,gridCPU->uz_twc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->uz_fwc0,gridCPU->uz_fwc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->lz_twc0,gridCPU->lz_twc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->lz_fsc0,gridCPU->lz_fsc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->lz_fpc0,gridCPU->lz_fpc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
		cudaMemcpy(gridGPU->adimc0,gridCPU->adimc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyHostToDevice );
}
extern "C" void sac_g_gpu_2_cpu(sac_grid *gridCPU, sac_grid* gridGPU){

		cudaMemcpy(gridCPU->prec,gridCPU->prec, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny*gridCPU->ntimes, cudaMemcpyHostToDevice );

		cudaMemcpy(gridCPU->grnd,gridGPU->grnd, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->surf,gridGPU->surf, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->tet,gridGPU->tet, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );

		cudaMemcpy(gridCPU->uz_twc,gridGPU->uz_twc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->uz_fwc,gridGPU->uz_fwc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->lz_twc,gridGPU->lz_twc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->lz_fsc,gridGPU->lz_fsc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->lz_fpc,gridGPU->lz_fpc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->adimc,gridGPU->adimc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );

		/*cudaMemcpy(gridCPU->uz_twc0,gridGPU->uz_twc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->uz_fwc0,gridGPU->uz_fwc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->lz_twc0,gridGPU->lz_twc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->lz_fsc0,gridGPU->lz_fsc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->lz_fpc0,gridGPU->lz_fpc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );
		cudaMemcpy(gridCPU->adimc0,gridGPU->adimc, sizeof(_TDATA)*gridCPU->nx*gridCPU->ny, cudaMemcpyDeviceToHost );*/
}
__global__ void sac_g_assign_fields(sac_stream_properties* propsGPU,sac_stream_properties propsCPU){
	//*propsGPU=propsCPU;
	*propsGPU=propsCPU;
}

extern "C" void sac_g_props(sac_grid *gridCPU, sac_stream_properties* propsGPU){
		//cudaMalloc((void**)&propsGPU,sizeof(sac_stream_properties));
		//cudaMemcpy(propsGPU,&gridCPU->prop,sizeof(sac_stream_properties),cudaMemcpyHostToDevice);
		sac_g_assign_fields<<<1,1>>>(propsGPU,gridCPU->prop);
}




extern "C" void sac_g_init_grid(sac_grid *grid, int ntimes, int nx, int ny, double dx){


	grid->nx=nx;
	grid->ny=ny;
	grid->nCells=nx*ny;
	grid->deltax=dx;
	grid->ntimes=ntimes;

	cudaMalloc((void**)&grid->prec,sizeof(_TDATA)*nx*ny*ntimes);
	cudaMalloc((void**)&grid->grnd,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->surf,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->tet,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->uz_twc,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->uz_fwc,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->lz_twc,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->lz_fsc,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->lz_fpc,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->adimc,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->uz_twc0,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->uz_fwc0,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->lz_twc0,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->lz_fsc0,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->lz_fpc0,sizeof(_TDATA)*nx*ny);
	cudaMalloc((void**)&grid->adimc0,sizeof(_TDATA)*nx*ny);
#if SILENTMODE==0
	printf("There are %d MB allocated\n",sizeof(double)*nx*ny*12/(1024*1024)+sizeof(double)*nx*ny*ntimes/(1024*1024));
#endif
}

__global__ void sac_g_grid_timestep(sac_grid stream, double t, double deltat,
					int index,		 //  Index of the precipitation
					double pedemand,	 //  evapotranspiration-demand
					int ncells,
					sac_stream_properties *props
					){	
	/**
	 *  @brief Definition of the variables
	 */

      double red,uzrat,duztwc,ratlz,ratlzt,
		twx,duz,dlzs,simpvt,roimp,dlzp,
		uztwm,uzfwm,lztwm,lzfsm,lzfpm,adimp,	/**< are maximum water storages capacities*/ 
		uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc,	/**< are total water storages*/
		uztwh,uzfwh,lztwh,lzfsh,lzfph,			/**< are unfrozen water storages*/
	  	eused,parea,tbf,bfcc,bfp,bfs,tci,
		ratio,adsur,addro,bfncc,
		bf,xx1,percm,perc,defr,check,perct,percf,
		hpl,ratlp,ratls,fracp,excess,percs,percp,
		surf,									/**< Surface flow */
		grnd,									/**< Subsurface flow */
		tet									/**< Total Evapotranspiration*/
		;

	  double aesc=0.0;

	  double aux;

	  double sbf,ssur,sur,sif,sperc,sdro,spbf;		/**<time interval sums.*/

	  double e1,e2,e3,e4,e5,xx,sfh,del;

	  int i,k;

	  double dinc; /**<dinc=length of each increment in days.*/
	  int  ninc; /**<ninc=number of time increments that the time interval
					is divided into for further
					soil-moisture accounting.  no one increment
					will exceed 5.0 millimeters of uzfwc+pav
				*/
	double  pinc; /**< pinc=amount of available moisture for each increment.*/
	double precip;
	/**
	*  @brief Maximum capacities
	*/
	uztwm=props->uz_twm;
	uzfwm=props->uz_fwm;
	adimp=props->adimp;
	lztwm=props->lz_twm;
	lzfsm=props->lz_fsm;
	lzfpm=props->lz_fpm;
	parea=1.0-props->pct_im-adimp;

	i = threadIdx.x+(blockIdx.x*blockDim.x);

	if(i<ncells){
		//printf("Time: %d Thread: %d\n",i,index);


		precip=stream.prec[index*stream.nCells+i];
		/**
		*  @brief Total water storages (state)
		*/
		 uztwc=stream.uz_twc0[i];
		 uzfwc=stream.uz_fwc0[i];
		 adimc=stream.adimc0[i];
		 lztwc=stream.lz_twc0[i];
		 lzfsc=stream.lz_fsc0[i];
		 lzfpc=stream.lz_fpc0[i];
		//printf("uztwc %.2e uzfwc %.2e adimc %.2e lztwc %.2e lzfsc %.2e lzfpc %.2e\n",uztwc, uzfwc, adimc, lztwc, lzfsc, lzfpc);

		/**
		*  @brief Unfrozen water storages. In this verion
		* it is assumed equal to the total water storages
		* keeping unfrozen water=total. It'd be desireble to
		* include free water movement:
		* ruzice & rlzice are reduction of free water movement 
		* based on kulik's theory: kfrz = kunfrz/(1+fgpm(4)*ice)**2
		*/    
		 //ruzice=uzk
		 //rlzice=lzsk
		 //ruzperc=1.0

		 uztwh=uztwc;
		 uzfwh=uzfwc;
		 lztwh=lztwc;
		 lzfsh=lzfsc;
		 lzfph=lzfpc;
	 
	 
		/* uztwh=uztwc-s0->uz_twf;
		 uzfwh=uzfwc-s0->uz_fwf;
		 lztwh=lztwc-s0->lz_twf;
		 lzfsh=lzfsc-s0->lz_fsf;
		 lzfph=lzfpc-s0->lz_fpf;
		 */
		/**
		*  @brief 
		compute evapotranspiration loss for the time interval.
		edmnd is the et-demand for the time interval
			edmnd=ep*epdist(kint)
		adjust edmnd for effect of snow & forest cover.
			edmnd=(1.-(1.0-sacpar(17))*aesc)*edmnd
		Otherwise, compute et from upper zone
		*/
		//

		pedemand=(1.-(1.0-props->efc)*aesc)*pedemand;
		e1=pedemand*(uztwh/uztwm);
		// Residual evap depand
		red=pedemand-e1;
		uztwh=uztwh-e1;
    
		e2=0.0;
		if(uztwh>=0.0){
			uztwc=uztwc-e1;
			//  go to 220
		}
		if((uztwh<0.0||(uztwh>=0.0&&(uztwc/uztwm)<(uzfwc/uzfwm)))){ //this is jump to 220 and then to 225
			if(uztwh<0){
				//  e1 can not exceed uztwc
				e1=e1+uztwh;
				uztwh=0.0;
	
				// reduce total tension water by actual e1
				uztwc=uztwc-e1;
				if(uztwc<0.0){
					uztwc=0.0;
				}
			}
			red=pedemand-e1;
			if(uztwh>=0.0||uzfwh>=red){
					/*221*/
					if(uztwh<0.0){
					e2=red;
					uzfwc=uzfwc-e2;
					uzfwh=uzfwh-e2;
					red=0.0;
					/*220*/
					}
					if(uzfwh>=0&&(uztwc/uztwm)<=(uzfwc/uzfwm)){
						uzrat=(uztwc+uzfwc)/(uztwm+uzfwm);
						duztwc=uztwm*uzrat-uztwc;
						if(duztwc>uzfwh) duztwc=uzfwh; 
						//transfered water can not exceed unfrozen free water
						uztwc=uztwc+duztwc;
						uztwh=uztwh+duztwc;
						uzfwc=uzfwc-duztwc;
						uzfwh=uzfwh-duztwc;
					}
			}else{
				e2=uzfwh;
				uzfwh=0.0;
				uzfwc=uzfwc-e2;
				if(uzfwc<0.0){
					uzfwc=0.0;
				}
				red=red-e2;
			}
		}

		/*225*/
		if (uztwc<TOL5) {
		   uztwc=0.0;
		   uztwh=0.0;
		} 
		if (uzfwc<TOL5) {
		   uzfwc=0.0;
		   uzfwh=0.0;
		} 

		  e3=red*(lztwh/(uztwm+lztwm));
		  lztwh=lztwh-e3;
		  if(lztwh>=0.0){
			lztwc=lztwc-e3;
		  }else{
			e3=e3+lztwh;
			lztwh=0.0;
			// reduce total tension water by e3
			lztwc=lztwc-e3;
		  }

		  /*226*/
		  // sacpar(16)=props->rserv;
		  ratlzt=lztwc/lztwm;
		  ratlz=(lztwc+lzfpc+lzfsc-props->rserv)/
				(lztwm+lzfpm+lzfsm-props->rserv);
		  if(ratlzt>=ratlz){
			  //NULL
		  }else{
			//resupply lower zone tension water from lower
			//zone free water if more water available there.
			del=(ratlz-ratlzt)*lztwm;
			//only unfrozen water can be transfered
			sfh=lzfsh+lzfph;
			if(del>sfh){
				del=sfh;
			}
			lzfsh=lzfsh-del;
			if(lzfsh>=0.0){
				// transfer from lzfsc to lztwc.      
			   lzfsc=lzfsc-del;
			}else{
				//if transfer exceeds lzfsc then remainder comes from lzfpc
			   lzfpc=lzfpc+lzfsh;
			   lzfph=lzfph+lzfsh;
			   xx=lzfsh+del;
			   lzfsc=lzfsc-xx;
			   lzfsh=0.0;
			}
			lztwc=lztwc+del;
			lztwh=lztwh+del;
		  }
		  if (lztwc<TOL5){
		   lztwc=0.0;
		   lztwh=0.0; 
		  }
		  e5=e1+(red+e2)*((adimc-e1-uztwc)/(uztwm+lztwm));
		//adjust adimc,additional impervious area storage, for evaporation.
		  adimc=adimc-e5;
		  if(adimc<0.0){
			//e5 can not exceed adimc.
			e5=e5+adimc;
			adimc=0.0;
		  }
		  e5=e5*adimp;
		  //compute percolation and runoff amounts.
		  twx=precip+uztwc-uztwm ;
		  if(twx<0.0){
			  uztwc=uztwc+precip;
			//adjust unfrozen tension water
			  uztwh=uztwh+precip;     
			  twx=0.0;
		  }else{
			uztwh=uztwh+(uztwm-uztwc);
			uztwc=uztwm;
		  }
		  adimc=adimc+precip-twx;
		  //compute impervious area runoff.
		  roimp=precip*props->pct_im;
		  //roimp is runoff from the minimum impervious area.

			/******************************
			NOTE: In the original code, this was not
			initialized. Required for the
			following intruction. Dit it implicitly init?
			*******************************/
	  		  simpvt=0.0;
			/*******************************/
			simpvt=simpvt+roimp;

		  //initialize time interval sums.
		  sbf=0.0;
		  ssur=0.0;
		  sif=0.0;
		  sperc=0.0;
		  sdro=0.0;
		  spbf=0.0;
		  //determine computational time increments for the basic time
		  //interval
		  //ninc=1.0+0.2*(uzfwc+twx)
		  //percolate unfrozen water only
		  ninc=1+(int)(0.2*(uzfwh+twx));
		  /**
		  *  @brief 
		  ninc=number of time increments that the time interval
		  is divided into for further
		  soil-moisture accounting.  no one increment
		  will exceed 5.0 millimeters of uzfwc+pav
		  */
		  dinc=1.0/(double)ninc*deltat;
		  //dinc=length of each increment in days.
		  pinc=twx/ninc;
		  pinc=(twx/(1.0+(0.2*(uzfwh+twx))));

		  //introduced reduction (ruzice & rlzice) due frozen ground
		  //however, primary runoff is unchanged 
		  duz =1.0-pow((1.0-props->uzk),dinc);
		  dlzs=1.0-pow((1.0-props->lzsk),dinc);
		  dlzp=1.0-pow((1.0-props->lzpk),dinc);

		//****************************************************************
		//----- end incremental do loop for the time interval. -------
		//****************************************************************
		  for(k=0;k<ninc;k++){
			adsur=0.0;
			// compute direct runoff (from adimp area).
			ratio=(adimc-uztwc)/lztwm;
			if (ratio<0.0){
				ratio=0.0;
			}
			 addro=pinc*(ratio*ratio);
			//baseflow from unfrozen water only   
			bf=lzfph*dlzp;
			lzfph=lzfph-bf;
			if (lzfph>TOL4){
			 lzfpc=lzfpc-bf;
			}else{
			  bf=bf+lzfph;
			  lzfph=0.0;
			  lzfpc=lzfpc-bf;
			  if(lzfpc<=TOL4){
				  lzfpc=0.0;
			  }
			}
			sbf=sbf+bf;
			spbf=spbf+bf;
			// supplamental flow from unfrozen water only (note, dlzs
			// note, dlzs is reduced due frozen ground
			bf=lzfsh*dlzs;
			lzfsh=lzfsh-bf;
			if(lzfsh>TOL4){
				lzfsc=lzfsc-bf;
			}else{
				bf=bf+lzfsh;
				lzfsh=0.0;
				lzfsc=lzfsc-bf;
				if(lzfsc<=TOL4){
					lzfsc=0.0;
				}
			}
			sbf=sbf+bf;
			//compute percolation-if no water available then skip
			xx1=pinc+uzfwh;
			if(xx1<=0.01){
				uzfwc=uzfwc+pinc;
				//add to unfrozen water also
				uzfwh=uzfwh+pinc;
			}else{
				 percm=lzfpm*dlzp+lzfsm*dlzs;
				 perc=percm*(uzfwh/uzfwm);
				 //if(ivers .ne. 0) perc=perc*ruzperc
				 defr=1.0-((lztwc+lzfpc+lzfsc)/(lztwm+lzfpm+lzfsm));
				 aux=pow(defr,props->r_exp);
				 perc=perc*(1.0+props->z_perc*pow(defr,props->r_exp));
				  if(perc>=uzfwh){
					  perc=uzfwh;
				  }
					uzfwc=uzfwc-perc;
					uzfwh=uzfwh-perc; 

				// check to see if percolation exceeds lower zone deficiency.
				check=lztwc+lzfpc+lzfsc+perc-lztwm-lzfpm-lzfsm;
				if(check>0.0){
					  perc=perc-check;
					  uzfwc=uzfwc+check;
					//adjust unfrozen starage also
					  uzfwh=uzfwh+check;  
				}
				//sperc is the time interval summation of perc
					sperc=sperc+perc;
				/*
					 compute interflow and keep track of time interval sum.
					 note...pinc has not yet been added
					 del=uzfwc*duz*fi
				 interflow also reduced due frofen ground (duz reduced by ruzice)
				 additional reduction due impervious frozen areas (fi) is optional
				 in the new version. basic option is fi=1
				*/
				del=uzfwh*duz;
				sif=sif+del;
				uzfwc=uzfwc-del;
				//adjust unfrozen storage also
				uzfwh=uzfwh-del;
	
				/*distribe percolated water into the lower zones
				 tension water must be filled first except for the pfree area.
				 perct is percolation to tension water and percf is percolation
				 going to free water.*/
				 perct=perc*(1.0-props->pct_free);
				 xx1=perct+lztwc;
				 if (xx1<=lztwm){
					 lztwc=lztwc+perct;
					 lztwh=lztwh+perct;      
					 percf=0.0;
				 }else{
					percf=perct+lztwc-lztwm;
					// change unfrozen water storage
					lztwh=lztwh+lztwm-lztwc;  
					lztwc=lztwm;
				 }
				 percf=percf+perc*props->pct_free;
				 if(percf!=0.0){ //In other case, goto 245
					// hpl is the relative size of the primary storage
					// as compared with total lower zone free water storage.
					 hpl=lzfpm/(lzfpm+lzfsm);
					  if(lzfpm != 0.){
					   ratlp=lzfpc/lzfpm;
					  }else{
					   ratlp = 1.0;
					  }
					  if(lzfsm != 0.){
					   ratls=lzfsc/lzfsm;
					  }else{
					   ratls = 1.0;
					  }
        
					// ratlp and ratls are content to capacity ratios, or
					// in other words, the relative fullness of each storage
					fracp=(hpl*2.0*(1.0-ratlp))/((1.0-ratlp)+(1.0-ratls));
					// fracp is the fraction going to primary.
					  if (fracp>1.0){
						  fracp=1.0;
					  }
					  percp=percf*fracp;
					  percs=percf-percp;
				// percp and percs are the amount of the excess
				// percolation going to primary and supplemental
				// storges,respectively.
					  lzfsc=lzfsc+percs;
					  if(lzfsc<=lzfsm){
						 lzfsh=lzfsh+percs;
					  }else{
						  percs=percs-lzfsc+lzfsm;
						//adjust unfrozen storage also
						  lzfsh=lzfsh+percs;
						  lzfsc=lzfsm;
					  }
					  lzfpc=lzfpc+(percf-percs);
					  if (lzfpc<=lzfpm){
						  lzfph=lzfph+(percf-percs);
					  }else{
						  excess=lzfpc-lzfpm;
						  lztwc=lztwc+excess;
						//  adjust unfrozen storages also
						  lztwh=lztwh+excess;
						  lzfph=lzfph+(percf-percs)-excess;
						  lzfpc=lzfpm;
					  }
				 }
				 /*245*/
				 if(pinc!=0.0){
					 xx1=pinc+uzfwc;
					 if(xx1<=uzfwm){
						 uzfwc=uzfwc+pinc;
						 uzfwh=uzfwh+pinc;
					 }else{
						sur=pinc+uzfwc-uzfwm;
						uzfwc=uzfwm;
						//adjust unfrozen storage also
						uzfwh=uzfwh+pinc-sur;
						ssur=ssur+sur*parea;
						adsur=sur*(1.0-addro/pinc);
						//     adsur is the amount of surface runoff which comes
						//     from that portion of adimp which is not
						//     currently generating direct runoff.  addro/pinc
						//     is the fraction of adimp currently generating
						//     direct runoff.
						ssur=ssur+adsur*adimp;
					 }
				 }
			}
			/*249*/
			adimc=adimc+pinc-addro-adsur; 
			xx1=uztwm+lztwm;
			 if (adimc>xx1){
				addro=addro+adimc-xx1;
				adimc=xx1;
			 }
			sdro=sdro+addro*adimp;
			if (adimc<TOL5){
				  adimc=0.0;
			}
		  } // END
		//****************************************************************
		// ----- end incremental do loop for the time interval. -------
		//****************************************************************
		eused=e1+e2+e3;
		// eused is the et from parea which is 1.0-adimp-pctim
		sif=sif*parea;
		//separate channel component of baseflow
		// from the non-channel component
		tbf=sbf*parea;
		//tbf is total baseflow
		bfcc=tbf*(1.0/(1.0+props->side ));
		//bfcc is baseflow, channel component
		bfp=spbf*parea/(1.0+props->side);
		bfs=bfcc-bfp;
		if(bfs<0.0){
			bfs=0.0;
		}
		bfncc=tbf-bfcc;
		//bfncc is baseflow,non-channel component

		//add to monthly sums.
		//
		// compute total channel inflow for the time interval.
		tci=roimp+sdro+ssur+sif+bfcc;
		grnd = sif + bfcc;   // interflow is part of ground flow
		//grnd = bfcc           // interflow is part of surface flow
		surf = tci - grnd;
		//compute e4-et from riparian vegetation.
		e4=(pedemand-eused)*props->riva;
		//subtract e4 from channel inflow
		tci=tci-e4;
		if(tci<0.0){
		e4=e4+tci;
		tci=0.0	;	
		}
		grnd = grnd - e4;
		if (grnd<0.0){
		   surf = surf + grnd;
		   grnd = 0.0;
			if (surf<= 0.0){
				 surf = 0.0;
			}
		}
		//compute total evapotranspiration-tet
		eused=eused*parea;
		tet=eused+e5+e4;
		if (adimc<uztwc){
			adimc=uztwc;
		}
		/**
		*  @brief Total water storages (state)
		*/
		 stream.uz_twc[i]=uztwc;
		 stream.uz_fwc[i]=uzfwc;
		 stream.adimc[i]=adimc;
		 stream.lz_twc[i]=lztwc;
		 stream.lz_fsc[i]=lzfsc;
		 stream.lz_fpc[i]=lzfpc;	

		 stream.uz_twc0[i]=uztwc;
		 stream.uz_fwc0[i]=uzfwc;
		 stream.adimc0[i]=adimc;
		 stream.lz_twc0[i]=lztwc;
		 stream.lz_fsc0[i]=lzfsc;
		 stream.lz_fpc0[i]=lzfpc;		

		 stream.surf[i]+=surf;
		// printf("Surface: %lf\n",stream.surf[i]);

		 stream.grnd[i]=grnd;
		 stream.tet[i]=tet;
	 }
/*	 printf("%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2e %6.2e %6.2e\n",
		 s->uz_twc,
		  s->uz_fwc,
		  s->adimc,
		  s->lz_twc,
		  s->lz_fsc,
		  s->lz_fpc,
		  s->surf,
		  s->grnd,
		  s->tet
		  );
		  */
	 // Check negative values


}

