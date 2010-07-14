//reduced poizeil 1.10 (without Prandtl)

//divergence under small Reynolds

//structural functions(along rectangle with periodic condtions)
# include "mesh.h"

# define  EPS  1e-9
# define STEP (0.01)
# define uplevel 100000.
# define Noise 0.01
# include "macros.h"

# include "vars.h"

 double dt=5e-3,
		 vx[2][Nx+2][Ny+2][Nz+2],
		 vy[2][Nx+2][Ny+2][Nz+2],
		 vz[2][Nx+2][Ny+2][Nz+2];

 double p[2][Nx+2][Ny+2][Nz+2],
		 nut[Nx+2][Ny+2][Nz+2],
		 diver[Nx+1][Ny+1][Nz+1];
 double s_func[Nz+2][Nx/2+Ny/2+Nz];
 long num_points[Nx/2+Ny/2+Nz];
 double gamma ;

 double tm, lx,ly,lz, tmmax,divmax ;

 double p1, p2;//pressure on the ends

 bool konez;
 FILE  *fv, *fnu, *fsf;

FILE *prepout(char *x)  //opening of file to ff
{
FILE *ff;
  if ((ff = fopen (x,"w+"))==NULL)
	 {
		printf ("Can't open file %s !\n",x);
		exit(-1);
	 }
return(ff);
}

# define maschtab 10
# define modul(a,b,c) (sqrt(norm(a,b,c)))

void nut_by_flux(int ind) //calculating nu_turbulent by velocity fluctuations
{
double flux, maxflux = 0;
int i,j,k,kol = 0;
for(i=0;i<=Nx+1;i++)
    for(j=0;j<=Ny+1;j++)
        for(k=0;k<=Nz+1;k++)
            {
            flux = 0;
            if(i>0)
            	{flux += modul(vx[ind][i][j][k]-vx[ind][i-1][j][k],vy[ind][i][j][k]-vy[ind][i-1][j][k],vz[ind][i][j][k]-vz[ind][i-1][j][k]);
            	kol++;
               }
            if(i<Nx+1)
            	{flux += modul(vx[ind][i][j][k]-vx[ind][i+1][j][k],vy[ind][i][j][k]-vy[ind][i+1][j][k],vz[ind][i][j][k]-vz[ind][i+1][j][k]);
            	kol++;
               }
            if(j>0)
            	{flux += modul(vx[ind][i][j][k]-vx[ind][i][j-1][k],vy[ind][i][j][k]-vy[ind][i][j-1][k],vz[ind][i][j][k]-vz[ind][i][j-1][k]);
            	kol++;
               }
            if(j<Ny+1)
	            {flux += modul(vx[ind][i][j][k]-vx[ind][i][j+1][k],vy[ind][i][j][k]-vy[ind][i][j+1][k],vz[ind][i][j][k]-vz[ind][i][j+1][k]);
            	kol++;
               }
            if(k>0)
            	{flux += modul(vx[ind][i][j][k]-vx[ind][i][j][k-1],vy[ind][i][j][k]-vy[ind][i][j][k-1],vz[ind][i][j][k]-vz[ind][i][j][k-1]);
            	kol++;
               }
            if(k<Nz+1)
            	{flux += modul(vx[ind][i][j][k]-vx[ind][i][j][k+1],vy[ind][i][j][k]-vy[ind][i][j][k+1],vz[ind][i][j][k]-vz[ind][i][j][k+1]);
            	kol++;
               }
            flux /= kol;
            if(flux>maxflux) maxflux = flux;
            nut[i][j][k] = (1. + maschtab * flux)/Re;
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

void printing(int ind)
{
 int i,j,k;
 double epsp=0, epsvx=0, epsvy=0, epsvz=0;
 double vxmax=0, vymax=0, vzmax= 0;
// double avervx[Nz+2], avernu[Nz+2];

		for(i=1;i<=Nx;i++)
		  for(j=1;j<=Ny;j++)
			  for(k=1;k<=Nz;k++)
			 {
				 if(vxmax < fabs(vx[ind][i][j][k])) vxmax = fabs(vx[ind][i][j][k]);
				 if(vymax < fabs(vy[ind][i][j][k])) vymax = fabs(vy[ind][i][j][k]);
				 if(vzmax < fabs(vz[ind][i][j][k])) vzmax = fabs(vz[ind][i][j][k]);
				 if(epsvx < fabs(vx[ind][i][j][k]-vx[(ind+1)%2][i][j][k]))
					 epsvx = fabs(vx[ind][i][j][k]-vx[(ind+1)%2][i][j][k]);
				 if(epsvy < fabs(vy[ind][i][j][k]-vy[(ind+1)%2][i][j][k]))
					 epsvy = fabs(vy[ind][i][j][k]-vy[(ind+1)%2][i][j][k]);
				 if(epsvz < fabs(vz[ind][i][j][k]-vz[(ind+1)%2][i][j][k]))
					 epsvz = fabs(vz[ind][i][j][k]-vz[(ind+1)%2][i][j][k]);
				 if(epsp  < fabs(p[ind][i][j][k]-p[(ind+1)%2][i][j][k]))
					 epsp  = fabs(p[ind][i][j][k]-p[(ind+1)%2][i][j][k]);
			 }

/*	for(k=1;k<=Nz;k++)
		{
		avervx[k] = avernu[k] = 0;
		for(i=1;i<=Nx;i++)
			for(j=1;j<=Ny;j++)
				{
				avervx[k] += vx[(ind+1)%2][i][j][k];
				avernu[k] += nut[i][j][k];
				}
		avervx[k] /= Nx*Ny;
		avernu[k] /= Nx*Ny;
		}*/

		 clrscr();
		 printf("\r %9f    %e %e %e %e\n",tm,epsp,epsvx,epsvy,epsvz);
		 printf("              %e %e %e %e\n",vxmax,vymax,vzmax,divmax);
	//putting velocities to file
	fprintf(fv,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(fv,"%0.5f,",vx[(ind+1)%2][Nx/2][Ny/2][k]);
	fprintf(fv,"%0.5f}\n",vx[(ind+1)%2][Nx/2][Ny/2][Nz]);
	//putting viscosities to file
/*	fprintf(fnu,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(fnu,"%0.5f,",avernu[k]);
	fprintf(fnu,"%0.5f}\n",avernu[Nz]);*/
	//putting structural functions to file
//	struct_func(2,ind);
	fprintf(fsf,"{%0.6f",s_func[1][0]);
	for(k=1;k<=Nx/2+Ny/2+Nz-1&&num_points[k];k++)
		fprintf(fsf,",%0.6f",s_func[1][k]);
	fprintf(fsf,"}\n");

	 if(kbhit()&&getch()=='q') konez=true;
	 if(epsvx<EPS && epsvy<EPS && epsvz<EPS ||
			vxmax>uplevel || vymax>uplevel ||vzmax>uplevel)
			konez=true;

} //printing

void step_of_time(int ns)      //one of the time steps of calculation

{
 int n0 = ns%2, n1 = (ns+1)%2;
 int  i,j,k ;

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
						  ( nut[i][j][k]) );
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

	  for(i=1;i<=Nx;i++)
	 for(j=1;j<=Ny;j++)
		 for(k=1;k<=Nz;k++)
		 {
		 diver[i][j][k] = (vx[n1][i+1][j][k]-vx[n1][i-1][j][k])/(2*dx)+
				 (vy[n1][i][j+1][k]-vy[n1][i][j-1][k])/(2*dy)+
				 (vz[n1][i][j][k+1]-vz[n1][i][j][k-1])/(2*dz);
		 }

/*------------ Step III (counting pressure)---------------------------*/

		for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
			for(k=1;k<=Nz;k++)
		  p[n1][i][j][k] = p[n0][i][j][k] - dt*diver[i][j][k]/gamma;

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

}  //step_of_time

/*=========================================================================*/

void main()
{
 long ns;
 int i,j,k;

/*============ Initial condition ==============*/

  tmmax=500.0;
  tm = 0;

  lx = 2; lz = ly= 1; dx = lx/Nx; dy = ly/Ny; dz=lz/Nz;

  gamma = 0.01;

  p1 = 8*lx/(lz*Re) ; p2 = 0;

  ns = 0;

  //initial filling of arrays (parabole profil)
  for(i=0;i<=Nx+1;i++)
	for(j=0;j<=Ny+1;j++)
		for(k=0;k<=Nz+1;k++)
				 {
				 vx[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX
								+ 1 - 4./Nz/Nz*(k-1-Nz/2)*(k-1-Nz/2);
				 vy[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
				 vz[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
				 p[0][i][j][k] = p1+(i-0.5)*(p2-p1)/Nx;
				 nut[i][j][k]= 1./Re;
				 }
//	fill_num_points();

	fv = prepout("vv.dat");
	fprintf(fv,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(fv,"%0.5g,",vx[0][Nx/2][Ny/2][k]);
	fprintf(fv,"%0.5g}\n",vx[0][Nx/2][Ny/2][Nz]);

	fnu = prepout("nut.dat");
	fprintf(fnu,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(fnu,"%0.5g,",nut[Nx/2][Ny/2][k]);
	fprintf(fnu,"%0.5g}\n",nut[Nx/2][Ny/2][Nz]);

	fsf = prepout("strfunc.dat");

/*=================== Main block ====================================*/

 konez=false;

	while(tm<tmmax && ! konez)
	{
	 step_of_time(ns);
    nut_by_flux(ns);
	 divmax = 0;
	  for(i=1;i<=Nx;i++)
	 for(j=1;j<=Ny;j++)
		 for(k=1;k<=Nz;k++)
		 {
		 diver[i][j][k] = (vx[ns%2][i+1][j][k]-vx[ns%2][i-1][j][k])/(2*dx)+
				 (vy[ns%2][i][j+1][k]-vy[ns%2][i][j-1][k])/(2*dy)+
				 (vz[ns%2][i][j][k+1]-vz[ns%2][i][j][k-1])/(2*dz);
		 divmax=max(divmax,fabs(diver[i][j][k]));
		 }
//	 if((int)((tm+dt/2)/STEP)-(int)((tm-dt/2)/STEP))
	 {
	 printing((ns+1)%2);
	 }

	  ns++;
	  tm+=dt;
	}

	printf("\n The End \n");
	putch(getch());
	printf("\n vx[%u][%u][k]:\n",Nx/2,Ny/2);
	for(k=1;k<=Nz;k++)
		printf("%g\n",vx[ns%2][Nx/2][Ny/2][k]);
	putch(getch());
	fclose(fv);
	fclose(fsf);
}//main
