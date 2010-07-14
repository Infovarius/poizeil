//undone
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <conio.h>

# define  Nx   5
# define  Ny   5
# define  Nz   5

double dt;
FILE  *ff;

/*void prepout(char *x)  //opening of file to ff
{
  strcat(x,OUTNAME);
  if ((ff = fopen (x,"w+"))==NULL)
    {
      printf ("Can't open file %s !\n",x);
      exit(-1);
    }
} */

/*=========================================================================*/

void main()
{
int  i,j,k,n0,n1 ;

long ns;

double vx[2][Nx+2][Ny+2][Nz+2],
       vy[2][Nx+2][Ny+2][Nz+2],
       vz[2][Nx+2][Ny+2][Nz+2],
       p[2][Nx+2][Ny+2][Nz+2],
       nut[Nx+2][Ny+2][Nz+2],
       div[Nx+2][Ny+2][Nz+2];

 double p1, p2;//pressure on the ends

  double Re, gamma /*, varrand, sp[NTT]*/;

//  double epsp, epsvx, epsvy, gamma, disv, enk, enki;

  double tm, dx,dy,dz, lx,ly,lz, tmmax;//vxmax, vymax, dopmin, dopmax;

/*============ Initial condition ==============*/

  Re = 100.0;

  tmmax=50.0;
  tm = 0;
  dt = 1e-03;

  ly = 1; lx = ly; dx = lx/Nx; dy = ly/Ny;

  p1 = 8*lx/(ly*Re); p2 = 0;

  gamma = 0.001; //  ns = 0;

  for(i=0;i<Nx+2;i++)
	for(j=0;j<Ny+2;j++)
		for(k=0;k<Nz+2;j++)
		       {
		       p[0][i][j][k] = p1+(i-0.5)*(p2-p1)/Nx;
		       nut[i][j][k]= 0.01; //constant for a while
		       }

/*   prepout("vv.dat");
   for (j=1; j<Ny+1; j++)
     fprintf (ff,"%e %e %e\n",dy*(j-0.5),vx[0][Nx%2][j],vy[0][5][j]);
   fclose (ff);*/

/*------------ Step I Boundary condition --------------------*/

 //прилипание на горизонтальных площадках
 for(i=1;i<Nx+1;i++)
     for(j=1;j<Ny+1;j++)
	{
	 vx[0][i][j][0] = -vx[0][i][j][1];
	 vy[0][i][j][0] = -vy[0][i][j][1];
	 vz[0][i][j][0] = -vz[0][i][j][1];

	 vx[0][i][j][Nz+1] = -vx[0][i][j][Nz];
	 vy[0][i][j][Nz+1] = -vy[0][i][j][Nz];
	 vz[0][i][j][Nz+1] = -vz[0][i][j][Nz];
	}
 //периодические условия на вертикальных площадках
  for(i=1;i<Nx+1;i++)
     for(k=1;k<Nz+1;k++)
	 {
	 vx[0][i][0][k] = vx[0][i][Ny][k];
	 vy[0][i][0][k] = vy[0][i][Ny][k];
	 vz[0][i][0][k] = vz[0][i][Ny][k];

	 vx[0][i][Ny+1][k] = vx[0][i][1][k];
	 vy[0][i][Ny+1][k] = vy[0][i][1][k];
	 vz[0][i][Ny+1][k] = vz[0][i][1][k];
	 }

/*============ Main block ====================================*/

   while(tm<tmmax)
   {
    /* printf("\n %9f",tm);*/

     n0 = ns%2;
     n1 = (ns+1)%2;

/*------------ Step I (по скоростям без давления) ------------------*/

/*------------ Step II Boundary condition --------------------------*/

 //прилипание на горизонтальных площадках
 for(i=1;i<Nx+1;i++)
     for(j=1;j<Ny+1;j++)
	{
	 vx[n1][i][j][0] = -vx[n1][i][j][1];
	 vy[n1][i][j][0] = -vy[n1][i][j][1];
	 vz[n1][i][j][0] = -vz[n1][i][j][1];

	 vx[n1][i][j][Nz+1] = -vx[n1][i][j][Nz];
	 vy[n1][i][j][Nz+1] = -vy[n1][i][j][Nz];
	 vz[n1][i][j][Nz+1] = -vz[n1][i][j][Nz];
	}
 //периодические условия на вертикальных площадках
  for(i=1;i<Nx+1;i++)
     for(k=1;k<Nz+1;k++)
	 {
	 vx[n1][i][0][k] = vx[n1][i][Ny][k];
	 vy[n1][i][0][k] = vy[n1][i][Ny][k];
	 vz[n1][i][0][k] = vz[n1][i][Ny][k];

	 vx[n1][i][Ny+1][k] = vx[n1][i][1][k];
	 vy[n1][i][Ny+1][k] = vy[n1][i][1][k];
	 vz[n1][i][Ny+1][k] = vz[n1][i][1][k];
	 }

/*------------ Step II Divergention---------------------------*/

     for(i=1;i<Nx+1;i++)
	 for(j=1;j<Ny+1;j++)
	    for(k=1;k<Nz+1;k++)
	    {
	    div[i][j][k] = (vx[n1][i+1][j][k]-vx[n1][i-1][j][k])/(2*dx)+
		       (vy[n1][i][j+1][k]-vy[n1][i][j-1][k])/(2*dy)+
		       (vz[n1][i][j][k+1]-vz[n1][i][j][k-1])/(2*dz);
	    }

/*------------ Step III (по давлениям)---------------------------*/

      for(i=1;i<Nx+1;i++)
	  for(j=1;j<Ny+1;j++)
	      for(k=1;k<Nz+1;k++)
		  p[n1][i][j][k] = p[n0][i][j][k] - dt*div[i][j][k]/gamma;

/*------------ Step IV (давление - в скорости)-----------------------*/

      for(i=1;i<Nx+1;i++)
	  for(j=1;j<Ny+1;j++)
		{
	  p[n1][i][j][0] = p[n1][i][j][1];
	  p[n1][i][j][Nz+1] = p[n1][i][j][Nz];
		}

      for(j=0;j<Ny+2;j++)
	  for(k=0;k<Nz+2;k++)
	 {
	  p[n1][0][j][k]= 2*p1-p[n1][1][j][k];
	  p[n1][Nx+1][j][k]= 2*p2-p[n1][Nx][j][k];
	 }

      for(i=1;i<Nx+1;i++)
	  for(j=1;j<Ny+1;j++)
	      for(k=1;k<Nz+1;k++)
		     {
	    vx[n1][i][j][k] -= dt*(p[n1][i+1][j][k]-p[n1][i-1][j][k])/(2*dx);
	    vy[n1][i][j][k] -= dt*(p[n1][i][j+1][k]-p[n1][i][j-1][k])/(2*dy);
	    vz[n1][i][j][k] -= dt*(p[n1][i][j][k+1]-p[n1][i][j][k-1])/(2*dz);
		     }

/*--------------------Printing of results -----------------------------*/
    if(!(ns%100))
    {
    //printing
    }


     ns++;
     tm+=dt;
   }
   printf("\n The End \n");

}//main