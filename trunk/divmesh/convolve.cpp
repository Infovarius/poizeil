#include "mesh.h"
# define min_d min(dx,dy)

# define M 2         // number of edge structural points
# define koef 1      //koefficient in Fourier transform
double sinkx[M][Nx],sinky[M][Ny];
extern double    vx[2][Nx+2][Ny+2][Nz+2],
		 vy[2][Nx+2][Ny+2][Nz+2],
		 vz[2][Nx+2][Ny+2][Nz+2];
extern double tm, dx,dy,dz, lx,ly,lz, tmmax,divmax ;

double fvx[M],fvy[M],fvz[M];
double mul[M] = {0.01,0.01};     //dissipation parameters

void fourier(int nl, int nn)      //fourier transform with dissipation
{
double klim=2*M_PI/min_d;
int i,j,k;
//initial filling of sine tables
for(i=0;i<=M-1;i++)
   {
	for(j=0;j<=Nx-1;j++) sinkx[i][j] = sin(klim /pow(2.,i)*(j+0.5)*dx);
	for(j=0;j<=Ny-1;j++) sinky[i][j] = sin(klim /pow(2.,i)*(j+0.5)*dy);
   };
//direct transform
for(i=0;i<=M-1;i++)
   {
   fvx[i] = fvy[i] = fvz[i] = 0;
   for(j=1;j<=Nx;j++)
   	for(k=1;k<=Ny;k++)
      	{
      	fvx[i] += vx[nn][j][k][nl]*sinkx[i][j-1]*sinky[i][k-1]*dx*dy;
      	fvy[i] += vy[nn][j][k][nl]*sinkx[i][j-1]*sinky[i][k-1]*dx*dy;
      	fvz[i] += vz[nn][j][k][nl]*sinkx[i][j-1]*sinky[i][k-1]*dx*dy;
         };
   };
//inverse transform
for(j=1;j<=Nx;j++)
  	for(k=1;k<=Ny;k++)
       for(i=0;i<=M-1;i++)
           {
           vx[nn][j][k][nl] += (mul[i]-1)*koef*fvx[i]*sinkx[i][j-1]*sinky[i][k-1];
           vy[nn][j][k][nl] += (mul[i]-1)*koef*fvy[i]*sinkx[i][j-1]*sinky[i][k-1];
           vz[nn][j][k][nl] += (mul[i]-1)*koef*fvz[i]*sinkx[i][j-1]*sinky[i][k-1];
           };
}

void Gauss(int nl, int nn)      //fourier transform with dissipation
{
double klim=2*M_PI/min_d;
int i,j,k,jj,kk;
double all, mass, aa=1.;

for(j=1;j<=Nx;j++)
  	for(k=1;k<=Ny;k++)
        {
        all=0;mass=0;
       for(jj=1;jj<=Nx;jj++)
  	for(kk=1;kk<=Ny;kk++)
           {
           all += exp(-(pow((j-jj)*dx,2)+pow((k-kk)*dy,2))/aa)*vx[nn][jj][kk][nl];
           mass += exp(-(pow((j-jj)*dx,2)+pow((k-kk)*dy,2))/aa);
           }
           vx[nn][j][k][nl] = all/mass;
        all=0;mass=0;
       for(jj=1;jj<=Nx;jj++)
  	for(kk=1;kk<=Ny;kk++)
           {
           all += exp(-(pow((j-jj)*dx,2)+pow((k-kk)*dy,2))/aa)*vy[nn][jj][kk][nl];
           mass += exp(-(pow((j-jj)*dx,2)+pow((k-kk)*dy,2))/aa);
           }
           vy[nn][j][k][nl] = all/mass;
        all=0;mass=0;
       for(jj=1;jj<=Nx;jj++)
  	for(kk=1;kk<=Ny;kk++)
           {
           all += exp(-(pow((j-jj)*dx,2)+pow((k-kk)*dy,2))/aa)*vz[nn][jj][kk][nl];
           mass += exp(-(pow((j-jj)*dx,2)+pow((k-kk)*dy,2))/aa);
           }
           vz[nn][j][k][nl] = all/mass;
}
}
