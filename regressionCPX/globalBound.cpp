#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include "regression.hpp"
#include <fstream>

// time averaged observations: 4320 x 2160
#define FILE_NAME_OBS "climatology_time_averaged.nc"
// time averaged and interpolated model output: 4320 x 2160
#define FILE_NAME_MOD "chl_8d_fmcd_interpolated.nc"
#define NX 4320
#define NY 2160

// netCDF: error message and exiting with a non-zero status
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main()
{
   /* netCDF ID for the file and data variable */
   int ncid, varid, retv;

   // READ OBSERVATION GRID AND TIME AVERAGED OBSERVATIONS:
   double *lon, *lat, *chlaObs;
   lon = ( double* )calloc( NX, sizeof( double ) );
   lat = ( double* )calloc( NY, sizeof( double ) );
   chlaObs = ( double* )calloc( NX * NY, sizeof( double ) );

   if ( ( retv = nc_open( FILE_NAME_OBS, NC_NOWRITE, &ncid ) ) ) ERR( retv );
   if ( ( retv = nc_inq_varid(ncid, "LON", &varid ) ) ) ERR( retv );
   if ( ( retv = nc_get_var_double( ncid, varid, lon ) ) ) ERR( retv );
   printf( "read observation grid longitudes from %s!\n", FILE_NAME_OBS );

   if ( ( retv = nc_inq_varid(ncid, "LAT", &varid ) ) ) ERR( retv );
   if ( ( retv = nc_get_var_double( ncid, varid, lat ) ) ) ERR( retv );
   printf( "read observation grid latitudes from %s!\n", FILE_NAME_OBS );

   if ( ( retv = nc_inq_varid(ncid, "TAVEFEIN", &varid ) ) ) ERR( retv );
   if ( ( retv = nc_get_var_double( ncid, varid, chlaObs ) ) ) ERR( retv );
   if ( ( retv = nc_close( ncid ) ) ) ERR( retv );
   printf( "read observed time averaged CHLA from %s!\n", FILE_NAME_OBS );

   // READ TIME AVERAGED AND SPACIALLY INTERPOLATED MODEL OUTPUT:
   double *chlaMod;
   chlaMod = ( double* )calloc( NX * NY, sizeof( double ) );

   if ( ( retv = nc_open( FILE_NAME_MOD, NC_NOWRITE, &ncid ) ) ) ERR( retv );
   if ( ( retv = nc_inq_varid(ncid, "CHLA", &varid ) ) ) ERR( retv );
   if ( ( retv = nc_get_var_double( ncid, varid, chlaMod ) ) ) ERR( retv );
   if ( ( retv = nc_close( ncid ) ) ) ERR( retv );
   printf( "read time averaged model CHLA output from %s!\n", FILE_NAME_MOD );

   // MASK TRUSTED BOXES AND CALCULATE RMSE(OBS-SIM)
   // We will only trust chlorophyl values between 0 and 1 mg/m^3:
   int *mask;
   mask = ( int* )calloc( NX * NY, sizeof( int ) );
   for ( int i = 0; i < NX; i++ )
   {
      for ( int j = 0; j < NY; j++ )
      {
         mask[ NX * j + i ] = 0;
         if (    0 <= chlaObs[ NX * j  +  i ] && chlaObs[ NX *  j +  i ] <= 1 
              && 0 <= chlaMod[ NX * j + i ] && chlaMod[ NX * j + i ] <= 1 )
         {
            mask[ NX * j + i ] = 1;
         }
      }
   }

   // OPTIONALLY MASK ONLY OCEAN REGION (E.G., SOUTH ATLANTIC)
   // I.E., WE UNMASK ALL OTHER REGIONS
/*   for ( int i = 0; i < NX; i++ ) for ( int j = 0; j < NY; j++ )
   {
      if ( lat[ j ] > -60 ) // for Southern Ocean
//      if ( lon[ i ] < -69 || lon[ i ] > 21 || lat[ j ] < -60 || lat[ j ] > 0 ) // for South Atlantic
//      if ( lon[ i ] < -69 || lon[ i ] > 21 || lat[ j ] > -60 ) // for Southern Ocean adjacent to South Atlantic
      {
         mask[ NX * j + i ] = 0;
      }
   }
*/

   // CALCULATE RMSE(OBS-SIM)
   int nObs = 0;
   double sse = 0;
   for ( int i = 0; i < NX; i++ )
      for ( int j = 0; j < NY; j++ )
         if ( mask[ NX * j + i ] == 1 )
         {
            nObs++;
            double chlaDiff = chlaObs[ NX * j + i ] - chlaMod[ NX * j + i ];
            sse = sse + chlaDiff * chlaDiff;
         }
   printf( "%i boxes in desired region provide average observed value below 1\n", nObs );

   // CALCULATE LOWER BOUNDS ON RMSE(OBS-SIM)
   // SUBJECT TO ASSUMPTIONS FOR THE MODEL SIMULATIONS
   regression reg;
   double sseBound = 0;
   int n = 0;
   int N = 200;
   for ( int i = 0; i < NX; i++ )
   {
      // calculate local lower bound

      // generate local data vector:
      double *td; // latitudes
      double *xd; // observed chlorophyl values
      double *xr; // regression chlorophyl values
      td = ( double* )calloc( N, sizeof( double ) );
      xd = ( double* )calloc( N, sizeof( double ) );
      xr = ( double* )calloc( N, sizeof( double ) );

      int kSeg = 0;
      int nLocal = 0;
      double der = 0;
      int lastJ = 0;
      for ( int j = 0; j < NY; j++ ) if ( mask[ NX * j + i ] == 1 ) lastJ = j;
      for ( int j = 0; j < NY; j++ ) if ( mask[ NX * j + i ] == 1 )
      {
         td[ nLocal ] = lat[ j ];
         xd[ nLocal ] = chlaObs[ NX * j + i ];
         nLocal++;
         n++;

         // update der for current longitude segment
         if ( j < NY - 1 && mask[ NX * ( j + 1 ) + i ] )
         {
            double d2 = chlaMod[ NX * ( j + 1 ) + i ] - chlaMod[ NX * j + i ];
            double d = fabs( d2 / ( lat[ j + 1 ] - lat[ j ] ) );
            if ( d > der ) der = d;

            if ( j > 0 && mask[ NX * ( j - 1 ) + i ] )
            {
               double d1 = chlaMod[ NX * j + i ] - chlaMod[ NX * ( j-1 ) + i ];
            }
         }

         if ( nLocal == N - 1 || j == lastJ )
         {
            // longitude segment data and properties scanned,
            // apply lower bound method to local data and add value
            double sseBoundLocal;
            printf( "Consider chunk %i of %i-th longitude (%f degree): N = %i, dMax = %f\n", kSeg, i, lon[ i ], nLocal + 1, 1.4 * der );
            sseBoundLocal = reg.slopeReg( nLocal, -1.4 * der, 1.4 * der, td, xd, xr );
            sseBound = sseBound + sseBoundLocal;
            printf( "Local rmse bound = %f\n", sqrt( sseBoundLocal / nLocal ) );
            printf( "Current rmse bound = %f\n\n", sqrt( sseBound / n ) );
            kSeg++;
            nLocal = 0;
            der = 0;
         }
      }
      free( td );
      free( xd );
      free( xr );
   }
   if ( n != nObs )
   {
		printf( "Wrong number of observations used in lower bound calculation!" );
        exit( 1 );
   }
   double rmse = sqrt( sse / nObs );
   double rmseBound = sqrt( sseBound / nObs );
   double ratio = double( floor( 1000 * ( rmseBound / rmse ) ) ) / 10;
   printf( "RMSE(Obs-Sim) = %f\n", rmse );
   printf( "Lower bound = %f\n", rmseBound );
   printf( "Ratio between bound and actual misfit: %.1f\%\n", ratio );

   free( lon );
   free( lat );
   free( chlaObs );
   free( chlaMod );
   free( mask );

   return 0;
}
