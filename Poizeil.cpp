//working version:
//gradient periodic condition+max v=1+print to file
/*rending in Mathematics by string:
ListPlot[#,PlotJoined->True,PlotRange->{0,1}]&
	  /@(l=ReadList["e:\\bc31\\mine\\vv.dat"]);*/

//turbulence constant nu
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <conio.h>

# define  Nx   10
# define  Ny   10
# define  Nz   10
# define  EPS  1e-15
# define min(x,y) ((x>y)?y:x)
# define abs(x)  ((x>0)?x:-x)

double dt,
		 vx[2][Nx+2][Ny+2][Nz+2],
		 vy[2][Nx+2][Ny+2][Nz+2],
		 vz[2][Nx+2][Ny+2][Nz+2];
FILE  *ff;

void prepout(char *x)  //opening of file to ff
{
  if ((ff = fopen (x,"w+"))==NULL)
	 {
		printf ("Can't open file %s !\n",x);
		exit(-1);
	 }
}

void velocitybounder(int ind)  //boundary conditions on velocities
{
int i,j,k;
			 //sticking conditions on horizontal surfaces
 for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
	{
	 vx[ind][i][j][0] = -vx[ind][i][j][1];
	 vy[ind][i][j][0] = -vy[ind][i][j][1];
	 vz[ind][i][j][0] = -vz[ind][i][j][1];

	 vx[ind][i][j][Nz+1] = -vx[ind][i][j][Nz];
	 vy[ind][i][j][Nz+1] = -vy[ind][i][j][Nz];
	 vz[ind][i][j][Nz+1] = -vz[ind][i][j][Nz];
	}
			 //periodic conditions on vertical surfaces
  for(i=1;i<=Nx;i++)
	  for(k=1;k<=Nz;k++)
	 {
	 vx[ind][i][0][k] = vx[ind][i][Ny][k];
	 vy[ind][i][0][k] = vy[ind][i][Ny][k];
	 vz[ind][i][0][k] = vz[ind][i][Ny][k];

	 vx[ind][i][Ny+1][k] = vx[ind][i][1][k];
	 vy[ind][i][Ny+1][k] = vy[ind][i][1][k];
	 vz[ind][i][Ny+1][k] = vz[ind][i][1][k];
	 }

			 //periodic conditions on stream surfaces
	for(j=1;j<=Ny;j++)
		for(k=1;k<=Nz;k++)
	 {
	 vx[ind][0][j][k] = vx[ind][Ny][j][k];
	 vy[ind][0][j][k] = vy[ind][Ny][j][k];
	 vz[ind][0][j][k] = vz[ind][Ny][j][k];

	 vx[ind][Ny+1][j][k] = vx[ind][1][j][k];
	 vy[ind][Ny+1][j][k] = vy[ind][1][j][k];
	 vz[ind][Ny+1][j][k] = vz[ind][1][j][k];
	 }
} //velocitybounder


/*=========================================================================*/

void main()
{
int  i,j,k,n0,n1 ;

long ns;

double p[2][Nx+2][Ny+2][Nz+2],
		 nut[Nx+2][Ny+2][Nz+2],
		 div[Nx+1][Ny+1][Nz+1];

 double p1, p2;//pressure on the ends

 double Re, gamma ,ksi, mn;

 double epsp, epsvx,epsvy,epsvz ;

 double tm, dx,dy,dz, lx,ly,lz, tmmax, vxmax,vymax,vzmax,divmax;

/*============ Initial condition ==============*/

  Re = 10.0;

  tmmax=50.0;
  tm = 0;
  dt = 1e-05;

  lx = 2; lz = ly= 1; dx = lx/Nx; dy = ly/Ny; dz=lz/Nz;

  p1 = 8*lx/(lz*Re); p2 = 0;

  gamma = 0.001; ksi = 0;
//nut->(nut+1)/Re
  mn = Re*(p2-p1)/lx*(1-coshl(ksi))/ksi;
//mn=0.01;
  ns = 0;

  for(i=0;i<=Nx+1;i++)
	for(j=0;j<=Ny+1;j++)
		for(k=0;k<=Nz+1;k++)
				 {
				 vx[0][i][j][k]=vy[0][i][j][k]=vz[0][i][j][k]=0;
				 p[0][i][j][k] = p1+(i-0.5)*(p2-p1)/Nx;
				 nut[i][j][k]= mn*(2*(k-0.5)*dz-1)/sinh(ksi*(2*(k-0.5)*dz-1));
				 }

/*  for(i=0;i<=Nx+1;i++)
	for(j=0;j<=Ny+1;j++)
		for(k=0;k<=Nz+1;k++)
				 nut[i][j][k]=mn*min(k,Nz-k)*min(k,Nz-k)*dz*dz*
					abs((vx[0][i][j][k+1]-vx[0][i][j][k-1])/2/dz);*/
	prepout("vv.dat");
	fprintf(ff,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(ff,"%g,",vx[0][Nx/2][Ny/2][k]);
	fprintf(ff,"%g}\n",vx[0][Nx/2][Ny/2][Nz]);

/*=================== Main block ====================================*/

	while(tm<tmmax)
	{

	  n0 = ns%2;
	  n1 = (ns+1)%2;

/*------------ Step 0 Boundary condition --------------------*/

velocitybounder(n0);

/*------------ Step I (velocities without pressure) ------------------*/
for(i=1;i<=Nx;i++)
	for(j=1;j<=Ny;j++)
		for(k=1;k<=Nz;k++)
			{
	 vx[n1][i][j][k] = vx[n0][i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-vx[n0][i][j][k] )*
								(vx[n0][i+1][j][k]-vx[n0][i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-vy[n0][i][j][k] )*
								(vx[n0][i][j+1][k]-vx[n0][i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-vz[n0][i][j][k] )*
								(vx[n0][i][j][k+1]-vx[n0][i][j][k-1])/(2*dz)+
			  ((vx[n0][i+1][j][k]-2*vx[n0][i][j][k]+vx[n0][i-1][j][k])/(dx*dx)+
				(vx[n0][i][j+1][k]-2*vx[n0][i][j][k]+vx[n0][i][j-1][k])/(dy*dy)+
				(vx[n0][i][j][k+1]-2*vx[n0][i][j][k]+vx[n0][i][j][k-1])/(dz*dz))*
						  (nut[i][j][k]) );
	 vy[n1][i][j][k] = vy[n0][i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-vx[n0][i][j][k] )*
								(vy[n0][i+1][j][k]-vy[n0][i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-vy[n0][i][j][k] )*
								(vy[n0][i][j+1][k]-vy[n0][i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-vz[n0][i][j][k] )*
								(vy[n0][i][j][k+1]-vy[n0][i][j][k-1])/(2*dz)+
			  ((vy[n0][i+1][j][k]-2*vy[n0][i][j][k]+vy[n0][i-1][j][k])/(dx*dx)+
				(vy[n0][i][j+1][k]-2*vy[n0][i][j][k]+vy[n0][i][j-1][k])/(dy*dy)+
				(vy[n0][i][j][k+1]-2*vy[n0][i][j][k]+vy[n0][i][j][k-1])/(dz*dz))*
						  (nut[i][j][k]) );
	 vz[n1][i][j][k] = vz[n0][i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-vx[n0][i][j][k] )*
								(vz[n0][i+1][j][k]-vz[n0][i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-vy[n0][i][j][k] )*
								(vz[n0][i][j+1][k]-vz[n0][i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-vz[n0][i][j][k] )*
								(vz[n0][i][j][k+1]-vz[n0][i][j][k-1])/(2*dz)+
			  ((vz[n0][i+1][j][k]-2*vz[n0][i][j][k]+vz[n0][i-1][j][k])/(dx*dx)+
				(vz[n0][i][j+1][k]-2*vz[n0][i][j][k]+vz[n0][i][j-1][k])/(dy*dy)+
				(vz[n0][i][j][k+1]-2*vz[n0][i][j][k]+vz[n0][i][j][k-1])/(dz*dz))*
						  (nut[i][j][k]) );
			}

/*------------ Step I Boundary condition --------------------*/

velocitybounder(n1);

/*------------ Step II Divergention---------------------------*/
	divmax=0;
	  for(i=1;i<=Nx;i++)
	 for(j=1;j<=Ny;j++)
		 for(k=1;k<=Nz;k++)
		 {
		 div[i][j][k] = (vx[n1][i+1][j][k]-vx[n1][i-1][j][k])/(2*dx)+
				 (vy[n1][i][j+1][k]-vy[n1][i][j-1][k])/(2*dy)+
				 (vz[n1][i][j][k+1]-vz[n1][i][j][k-1])/(2*dz);
		 divmax=max(divmax,abs(div[i][j][k]));
		 }

/*------------ Step III (counting pressure)---------------------------*/

		for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
			for(k=1;k<=Nz;k++)
		  p[n1][i][j][k] = p[n0][i][j][k] - dt*div[i][j][k]/gamma;

/*------------ Step IV (pressure to velocities)-----------------------*/

			  //free conditions on horizontal surfaces
		for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
		{
	  p[n1][i][j][0] = p[n1][i][j][1];
	  p[n1][i][j][Nz+1] = p[n1][i][j][Nz];
		}

				//periodic conditions on vertical surfaces
		for(i=1;i<=Nx;i++)
	  for(k=0;k<=Nz+1;k++)
		{
	  p[n1][i][0][k] = p[n1][i][Ny][k];
	  p[n1][i][Ny+1][k] = p[n1][i][1][k];
		}

				//gradient-periodic conditions on stream surfaces
		for(j=0;j<=Ny+1;j++)
	  for(k=0;k<=Nz+1;k++)
	 {
	  p[n1][0][j][k]= p[n1][Nx][j][k]-(p2-p1);
	  p[n1][Nx+1][j][k]= p[n1][1][j][k]+(p2-p1);
	 }

//correction of velocities by pressure
		for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
			for(k=1;k<=Nz;k++)
			  {
		 vx[n1][i][j][k] -= dt*(p[n1][i+1][j][k]-p[n1][i-1][j][k])/(2*dx);
		 vy[n1][i][j][k] -= dt*(p[n1][i][j+1][k]-p[n1][i][j-1][k])/(2*dy);
		 vz[n1][i][j][k] -= dt*(p[n1][i][j][k+1]-p[n1][i][j][k-1])/(2*dz);
			  }
//correction nut in Prandtl model
/*  for(i=0;i<=Nx+1;i++)
	for(j=0;j<=Ny+1;j++)
		for(k=0;k<=Nz+1;k++)
				 nut[i][j][k]=mn*min(k,Nz-k)*min(k,Nz-k)*dz*dz*
					abs((vx[n1][i][j][k+1]-vx[n1][i][j][k-1])/2/dz);*/

	  for(i=1;i<=Nx;i++)
	 for(j=1;j<=Ny;j++)
		 for(k=1;k<=Nz;k++)
		 {
		 div[i][j][k] = (vx[n1][i+1][j][k]-vx[n1][i-1][j][k])/(2*dx)+
				 (vy[n1][i][j+1][k]-vy[n1][i][j-1][k])/(2*dy)+
				 (vz[n1][i][j][k+1]-vz[n1][i][j][k-1])/(2*dz);
		 }

/*--------------------Printing of results -----------------------------*/
	 if(!(ns%100))
	 {
	 //printing
		epsvx = 0;
		epsvy = 0;
		epsvz = 0;
		vxmax = 0;
		vymax = 0;
		vzmax = 0;
		epsp  = 0;

		for(i=1;i<=Nx;i++)
		  for(j=1;j<=Ny;j++)
			  for(k=1;k<=Nz;k++)
			 {
				 if(vxmax < fabs(vx[n1][i][j][k])) vxmax = fabs(vx[n1][i][j][k]);
				 if(vymax < fabs(vy[n1][i][j][k])) vymax = fabs(vy[n1][i][j][k]);
				 if(vzmax < fabs(vz[n1][i][j][k])) vzmax = fabs(vz[n1][i][j][k]);
				 if(epsvx < fabs(vx[n1][i][j][k]-vx[n0][i][j][k]))
					 epsvx = fabs(vx[n1][i][j][k]-vx[n0][i][j][k]);
				 if(epsvy < fabs(vy[n1][i][j][k]-vy[n0][i][j][k]))
					 epsvy = fabs(vy[n1][i][j][k]-vy[n0][i][j][k]);
				 if(epsvz < fabs(vz[n1][i][j][k]-vz[n0][i][j][k]))
					 epsvz = fabs(vz[n1][i][j][k]-vz[n0][i][j][k]);
				 if(epsp  < fabs(p[n1][i][j][k]-p[n0][i][j][k]))
					 epsp  = fabs(p[n1][i][j][k]-p[n0][i][j][k]);
			 }
		 clrscr();
		 printf("\r %f    %e %e %e %e\n",tm,epsp,epsvx,epsvy,epsvz);
		 printf("             %e %e %e %e\n",vxmax,vymax,vzmax,divmax);
	fprintf(ff,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(ff,"%g,",vx[n0][Nx/2][Ny/2][k]);
	fprintf(ff,"%g}\n",vx[n0][Nx/2][Ny/2][Nz]);

	 if(kbhit()&&getch()=='q')break;
	 }

	 if(epsvx<EPS && epsvy<EPS && epsvz<EPS ||
			vxmax>1.0/EPS || vymax>1.0/EPS ||vzmax>1.0/EPS) break;

	  ns++;
	  tm+=dt;
	}
	printf("\n The End \n");
	putch(getch());
	printf("\n vx[%u][%u][k]:\n",Nx/2,Ny/2);
	for(k=1;k<=Nz;k++)
		printf("%g\n",vx[n0][Nx/2][Ny/2][k]);
	putch(getch());
	fclose(ff);
}//main
