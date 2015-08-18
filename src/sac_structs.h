/*--------------------------------------------------------
GPU-SAC 
GPU Accelerated Sacramento Soil Mosture Accounting Model
----------------------------------------------------------
Datatypes definition
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

#ifndef DEF_SAC_STRUCTS_
#define DEF_SAC_STRUCTS_

typedef struct sac_grid_ sac_grid;
typedef struct sac_conf_ sac_conf;
typedef struct sac_stream_ sac_stream;
typedef struct sac_stream_properties_ sac_stream_properties;
typedef struct sac_stream_state_ sac_stream_state;

typedef struct sm_stream_ sm_stream;

typedef struct sac_uz_ sac_uz;
typedef struct sac_lz_ sac_lz;
typedef struct sac_prec_ sac_prec;
typedef struct sac_evap_ sac_evap;


struct sac_prec_{
	int ntimes;				/**< Number of time series of precipitation*/
	double *t;				/**< Time*/
	double *e;				/**< Precipitation*/
};

struct sac_uz_{
	double iflow;			/**< Interflow*/
	double percolation;		/**< Percolation*/
};

struct sac_lz_{
	double sup_bflow;		/**< Supplemental Baseflow*/
	double prim_bflow;		/**< Primary Baseflow*/
	double groundloss;		/**< loss to groundwater*/
};

struct sac_stream_properties_{
	// Propierties
	/**
	 *  @brief Standard parameters
	 */
	int sasc_input_option;/**<
							 If present with value of 1, use
							input SASC time series, which is
							required to be present in input xml
							file. It is used in regular SACSMA
							calculation, and in the Frozen
							Ground calculation.
							 */
	int runoff_interval;/**<
						   Runoff component time series
							(ROCL) interval. Must be multiple
							of the input RAIM time series
							interval

						   */
	int sm_interval;/**<
					   Soil moisture storages time series
						(SMZC) interval. Must be multiple
						of the input RAIM time series
						interval
						   */
	int mape_input;/**< (NOT USED!)
						If absent or present with value of
						No, input MAPE time series is not
						used, regardless whether it is present
						or not. Then etDemand 12 values
						represent ET-demand 16th of each
						month, daily values are computed by
						linear interpolation
						   */
	double *et_demand_curve;/**<
						   if mape_input==1, 
						   ET-demand or PE-adjustment factor
							(the table contains 12 values).
							o If PE data used (i.e. MAPE_INPUT
							is 1), then values are PEadjustments;
							if PE data not used,
							then values represent ET-demand. 
							o Both 16th of each month (Jan.
							through Dec.; units of MM/day);
							daily values are computed by linear
							interpolation.
						   */
	double px_adj;/**<
						  Precipitation adjustment factor
						   */
	double pe_adj;/**<
						  ET-demand adjustment factor
						   */
	double uz_twm;/**< (mm)
						 Upper zone tension water capacity;
						   */
	double uz_fwm;/**<(mm)
						 Upper zone free water capacity;
						   */
	double lz_twm;/**<(mm)
						 Lower zone tension water capacity;
						   */
	double lz_fsm;/**< (mm)
						 Lower zone supplemental free water
						capacity;
					LZFSM are input as
					total values and not just as the
					visible (channel component) portion.
						   */
	double lz_fpm;/**< (mm)
					Lower zone primary free water
					capacity;
					LZFPM are input as
					total values and not just as the
					visible (channel component) portion.
						   */
	double uzk;/**<
						Fractional daily upper zone free
						water withdrawal rate
						   */	
	double lzsk;/**<
						Fractional daily supplemental
						withdrawal rate
						   */	
	double lzpk;/**<
						Fractional daily primary withdrawal
						rate
						   */	
	double pct_im;/**< (percent)
						Minimum impervious area;
						   */	
	double adimp;/**< (percent)
						 Additional impervious area;
						   */	
	double riva;/**< (percent)
						 Riparian vegetation area;
						   */	
	double efc;/**< (percent)
						  Effective forest cover;
						   */	
	double z_perc;/**<
						Maximum percolation rate
						   */		
	double r_exp;/**<
						Exponent for the percolation
						equation
						   */		
	double pct_free;/**< (percent)
						Percent/100 of percolated water
						which always goes directly to lower
						zone free water storages
						   */		
	double rserv;/**< (percent)
				 Percent/100 of lower zone free water
				which cannot be transferred to lower
				zone tension water
						   */		
	double side;/**<
				Ratio of non-channel baseflow (deep
				recharge) to channel (visible)
				baseflow
						   */			
	/**
	 *  @brief Frozen ground parameters
	 */
	int isfrozen;/**<
				Enables frozen ground calculation
					 */	
	int we_input_option;/**<
				Enables frozen ground calculation
					 */	
	int fr_out_interval;/**<
				 Output time series interval. Default
				value is RAIM interval;
					 */	
	double csoil;/**< (Cdeg^-1 h^-1)
				Bare ground frost coefficient for a
				given time interval;

					 */	
	double csnow;/**<
				Reduction in CSOIL per MM of
				snow water equivalent
					 */	
	double gch;/**< (Cdeg day^-1)
				Daily thaw rate from ground head;
					 */	
	double rthaw;/**< (Cdeg mm^-1)
				Thaw coefficient for water entering
				the soil;
					 */	
	double fr_temp;/**<(Cdeg)
				FI value above which there is no
				reduction in percolation or interflow
				withdrawal;
					 */	
	double satr;/**< (Cdeg^-1 h^-1)
				Reduction in percolation and
				interflow withdrawal per DEGC of
				FI below FIl under saturated soil
				conditions;
					 */	
	double fr_exp;/**<
				Exponent
					 */	
	int month0,monthN;
	int day0,dayN;

};
struct sac_stream_state_{

	int day, month;/**< 
					Day and month of the current state
							 */

	double uz_twc;/**< (mm)
					Upper zone tension water contents
							 */
	double uz_fwc;/**< (mm)
					Upper zone free water contents;
							 */	
	double lz_twc;/**< (mm)
					Lower zone tension water contents; 
							 */	
	double lz_fsc;/**< (mm)
					Lower zone free supplemental contents; 
							 */	
	double lz_fpc;/**< (mm)
					Lower zone free primary contents;
							 */	
	double adimc;/**< (mm)
					Tension water contents of the ADIMP area;
					If not known then use UZTWC+LZTWC
							 */	
	double uz_twf;/**< (mm)
					Upper zone tension water contents
							 */
	double uz_fwf;/**< (mm)
					Upper zone free water contents;
							 */	
	double lz_twf;/**< (mm)
					Lower zone tension water contents; 
							 */	
	double lz_fsf;/**< (mm)
					Lower zone free supplemental contents; 
							 */	
	double lz_fpf;/**< (mm)
					Lower zone free primary contents;
							 */	
	double adimf;/**< (mm)
					Tension water contents of the ADIMP area;
					If not known then use UZTWC+LZTWC
							 */	
	double fgix;/**< (cdeg)
					Initial value of the frost index; units of DEGC. If
					starting with no frozen ground the value is zero.
							 */	

	double sur_ro;/**< (mm)
					Surface runoff
							 */		

	double dir_ro;/**< (mm)
					direct runoff
							 */		

	double imp_ro;/**< (mm)
					impervious runoff
							 */		

	double sup_ro;/**< (mm)
					supplemental runoff
							 */		

	double pri_ro;/**< (mm)
					primary runoff
							 */	
	double infw;/**< (mm)
					primary runoff
							 */	
	int unit;/**<
					0 metric, 1 english
							 */		
	double surf;/**< ( )
					surface flow
							 */		
	double grnd;/**< ( )
					subsurface flow
							 */		
	double tet;/**< ( )
					total evapotranspiration
							 */	



};

struct sac_grid_{

	sac_stream_properties prop;/**<
					Stream properties for the sacramento model
					*/
	sac_stream_properties *g_prop;/**<
					Stream properties for the sacramento model
							 */	
	int nCells;
	int ntimes;
	int nx; 
	int ny;
	double deltax;

	double *position; /*
					  The matrix position is mapped as a vector of nx*ny components where
						it is accesed to the position (i,j) as
						position[i*nx+j]
					  */
	double *grnd;
	double *surf;
	double *tet;

	double *uz_twc; /**< (mm)
					Upper zone tension water contents
							 */
	double *uz_fwc;/**< (mm)
					Upper zone free water contents;
							 */	
	double *lz_twc;/**< (mm)
					Lower zone tension water contents; 
							 */	
	double *lz_fsc;/**< (mm)
					Lower zone free supplemental contents; 
							 */	
	double *lz_fpc;/**< (mm)
					Lower zone free primary contents;
							 */	
	double *adimc;/**< (mm)
					Tension water contents of the ADIMP area;
					If not known then use UZTWC+LZTWC
					 */
	double *uz_twc0; /**< (mm)
					Upper zone tension water contents
							 */
	double *uz_fwc0;/**< (mm)
					Upper zone free water contents;
							 */	
	double *lz_twc0;/**< (mm)
					Lower zone tension water contents; 
							 */	
	double *lz_fsc0;/**< (mm)
					Lower zone free supplemental contents; 
							 */	
	double *lz_fpc0;/**< (mm)
					Lower zone free primary contents;
							 */	
	double *adimc0;/**< (mm)
					Tension water contents of the ADIMP area;
					If not known then use UZTWC+LZTWC
					 */
	double *prec;	/**<
					  The matrix position is mapped as a vector of ntimes*nx*ny components where
						it is accesed to the position (i,j) as
						position[i*nx+j]
					*/
	sac_stream_state *records;/**<
					records for each variable 
							 */	
	double final_t;
	double delta_t;
	double delta_t_dump;							 
};



struct sac_stream_{

	sac_stream_properties prop;/**<
					Stream properties for the sacramento model
							 */	
	sac_stream_state state0;/**<
					Stream state at previous time t-Delta t
							 */	
	sac_stream_state state;/**<
					Stream state at present time t
							 */	
	sac_stream_state *grid_state0;/**<
					  The matrix grid_state0 is mapped as a vector of nx*ny components where
						it is accesed to the position (i,j) as
						position[i*nx+j]
					  */
	sac_stream_state *grid_state;/**<
					  The matrix grid_state is mapped as a vector of nx*ny components where
						it is accesed to the position (i,j) as
						position[i*nx+j]
					  */
	sac_grid *grid;/**<
					Grid information
							 */	

	sac_prec prec;/**<
					Precipitation series
							 */	
	int nrecords;

	sac_stream_state *records;/**<
					records for each variable 
							 */	


};

struct sm_stream_{
	// Propierties
	double wmax;			/**< Capacity to hold water (mm)*/
	double alpha;			/**< Is de inverse of the response time in the baseflow*/
	double mu;              /**< Dimensionless parameter determines the portion of the
							subsurface flow that becomes baseflow in the channels draining 
							out from the area of interes*/
	double Ep;				/**< Potential evaportranspiration rate*/
	double Z0;				/**< Soil water capacity */
	double I;				/**< Fraction of baseflow lost to deep aquifers*/
	double m;
	// Main Members
	sac_uz uz;
	sac_lz lz;
	sac_prec prec;
	double sr;				/**< Surface Runoff*/
	double w0;				/**< Soil water content at time t-1*/
	double w;				/**< Soil water content at time t*/

	double Q;				/**< Flow*/
};

#endif