//reduced poizeil 1.10 (without Prandtl)

//divergence under small Reynolds

//structural functions(along rectangle with periodic condtions)
# include "mesh.h"
typedef double mesh_arr[Nx+2][Ny+2][Nz+2];

# define  EPS  1e-9
# define STEP (0.1)
# define uplevel 10000.
# define Noise 0.0

# include "macros.h"

# define maschtab 00
# define modul(a,b,c) (sqrt(norm(a,b,c)))

//# include "integrating.cpp"

 double Re = 100000000;

 double dt=5e-4,
		 vx[2][Nx+2][Ny+2][Nz+2],
		 vy[2][Nx+2][Ny+2][Nz+2],
		 vz[2][Nx+2][Ny+2][Nz+2];

 double p[2][Nx+2][Ny+2][Nz+2],
		 nut[Nx+2][Ny+2][Nz+2];

 double s_func[Nz+2][Nx/2+Ny/2+Nz];
 long num_points[Nx/2+Ny/2+Nz];
 double gamma ;

 double tm, dx,dy,dz, lx,ly,lz, tmmax=500,divmax ;

 double p1, p2;//pressure on the ends

 bool konez;
 FILE  *fv, *fnu, *fsf;

FILE *prepout(char *x, int mode)  //opening of file to ff
{
FILE *ff;
if(mode) ff=fopen(x,"a");
else if ((ff = fopen (x,"w+"))==NULL)
	 {
		printf ("Can't open file %s !\n",x);
		exit(-1);
	 }
return(ff);
}


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
dt = 0.2*dx * dx*Re/(1.+ maschtab*maxflux);
if(dt>0.1) dt = 0.1;
}


void velocitybounder(mesh_arr Vx, mesh_arr Vy, mesh_arr Vz)  //boundary conditions on velocities
{
int i,j,k;
			 //sticking conditions on horizontal surfaces
 for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
	{
	 Vx[i][j][0] = -Vx[i][j][1];
	 Vy[i][j][0] = -Vy[i][j][1];
	 Vz[i][j][0] = -Vz[i][j][1];

	 Vx[i][j][Nz+1] = -Vx[i][j][Nz];
	 Vy[i][j][Nz+1] = -Vy[i][j][Nz];
	 Vz[i][j][Nz+1] = -Vz[i][j][Nz];
	}
			 //periodic conditions on vertical surfaces
  for(i=1;i<=Nx;i++)
	  for(k=1;k<=Nz;k++)
	 {
	 Vx[i][0][k] = Vx[i][Ny][k];
	 Vy[i][0][k] = Vy[i][Ny][k];
	 Vz[i][0][k] = Vz[i][Ny][k];

	 Vx[i][Ny+1][k] = Vx[i][1][k];
	 Vy[i][Ny+1][k] = Vy[i][1][k];
	 Vz[i][Ny+1][k] = Vz[i][1][k];
	 }

			 //periodic conditions on stream surfaces
	for(j=1;j<=Ny;j++)
		for(k=1;k<=Nz;k++)
	 {
	 Vx[0][j][k] = Vx[Ny][j][k];
	 Vy[0][j][k] = Vy[Ny][j][k];
	 Vz[0][j][k] = Vz[Ny][j][k];

	 Vx[Ny+1][j][k] = Vx[1][j][k];
	 Vy[Ny+1][j][k] = Vy[1][j][k];
	 Vz[Ny+1][j][k] = Vz[1][j][k];
	 }
} //velocitybounder

void printing(int ns)
{
 int i,j,k;
 int ind = (ns+1)%2;
 double epsp=0, epsvx=0, epsvy=0, epsvz=0;
 double vxmax=0, vymax=0, vzmax= 0;
 double avervx[Nz+2], avernu[Nz+2];

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

for(k=1;k<=Nz;k++)
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
		}

clrscr();
printf("\r %9f    %e %e %e %e\n",tm,epsvx,epsvy,epsvz,epsp);
printf("              %e %e %e %e\n",vxmax,vymax,vzmax,dt);
	//putting velocities to file
        fv = prepout("vv.dat",ns);
	fprintf(fv,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(fv,"%0.5f,",vx[(ind+1)%2][Nx/2][Ny/2][k]);
	fprintf(fv,"%0.5f}\n",vx[(ind+1)%2][Nx/2][Ny/2][Nz]);
        fclose(fv);
	//putting viscosities to file
        fnu = prepout("nut.dat",ns);
	fprintf(fnu,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(fnu,"%0.10f,",avernu[k]);
	fprintf(fnu,"%0.10f}\n",avernu[Nz]);
        fclose(fnu);
	//putting structural functions to file
/*	struct_func(2,ind);
	fprintf(fsf,"{%0.6f",s_func[1][0]);
	for(k=1;k<=Nx/2+Ny/2+Nz-1&&num_points[k];k++)
		fprintf(fsf,",%0.6f",s_func[1][k]);
	fprintf(fsf,"}\n");*/

	 if(kbhit()&&getch()=='q') konez=true;
	 if(epsvx<EPS && epsvy<EPS && epsvz<EPS ||
			vxmax>uplevel || vymax>uplevel ||vzmax>uplevel)
			konez=true;
         if(dt>0.5*min_d/modul(vxmax,vymax,vzmax))
               dt = 0.5*min_d/modul(vxmax,vymax,vzmax);
         if(dt<1e-20) konez=true;
         if(dt>0.1) dt = 0.1;
} //printing

void next_layer(  mesh_arr VxOld, mesh_arr VxNew,
                  mesh_arr VyOld, mesh_arr VyNew,
                  mesh_arr VzOld, mesh_arr VzNew,
                  mesh_arr pOld,  mesh_arr pNew     )
     //calculating next velocities
{
 int  i,j,k ;
 double diver[Nx+1][Ny+1][Nz+1];

/*------------ Step 0 Boundary condition --------------------*/

velocitybounder(VxOld, VyOld, VzOld);

/*------------ Step I (velocities without pressure) ------------------*/
for(i=1;i<=Nx;i++)
	for(j=1;j<=Ny;j++)
		for(k=1;k<=Nz;k++)
			{
	 VxNew[i][j][k] = VxOld[i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-VxOld[i][j][k] )*
								(VxOld[i+1][j][k]-VxOld[i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-VyOld[i][j][k] )*
								(VxOld[i][j+1][k]-VxOld[i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-VzOld[i][j][k] )*
								(VxOld[i][j][k+1]-VxOld[i][j][k-1])/(2*dz)+
			  ((VxOld[i+1][j][k]-2*VxOld[i][j][k]+VxOld[i-1][j][k])/(dx*dx)+
				(VxOld[i][j+1][k]-2*VxOld[i][j][k]+VxOld[i][j-1][k])/(dy*dy)+
				(VxOld[i][j][k+1]-2*VxOld[i][j][k]+VxOld[i][j][k-1])/(dz*dz))*
						  (nut[i][j][k]) );
	 VyNew[i][j][k] = VyOld[i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-VxOld[i][j][k] )*
								(VyOld[i+1][j][k]-VyOld[i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-VyOld[i][j][k] )*
								(VyOld[i][j+1][k]-VyOld[i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-VzOld[i][j][k] )*
								(VyOld[i][j][k+1]-VyOld[i][j][k-1])/(2*dz)+
			  ((VyOld[i+1][j][k]-2*VyOld[i][j][k]+VyOld[i-1][j][k])/(dx*dx)+
				(VyOld[i][j+1][k]-2*VyOld[i][j][k]+VyOld[i][j-1][k])/(dy*dy)+
				(VyOld[i][j][k+1]-2*VyOld[i][j][k]+VyOld[i][j][k-1])/(dz*dz))*
						  ( nut[i][j][k]) );
	 VzNew[i][j][k] = VzOld[i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-VxOld[i][j][k] )*
								(VzOld[i+1][j][k]-VzOld[i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-VyOld[i][j][k] )*
								(VzOld[i][j+1][k]-VzOld[i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-VzOld[i][j][k] )*
								(VzOld[i][j][k+1]-VzOld[i][j][k-1])/(2*dz)+
			  ((VzOld[i+1][j][k]-2*VzOld[i][j][k]+VzOld[i-1][j][k])/(dx*dx)+
				(VzOld[i][j+1][k]-2*VzOld[i][j][k]+VzOld[i][j-1][k])/(dy*dy)+
				(VzOld[i][j][k+1]-2*VzOld[i][j][k]+VzOld[i][j][k-1])/(dz*dz))*
						  (nut[i][j][k]) );
			}

/*------------ Step I Boundary condition --------------------*/

velocitybounder(VxNew, VyNew, VzNew);

/*------------ Step II Divergention---------------------------*/

  for(i=1;i<=Nx;i++)
	 for(j=1;j<=Ny;j++)
		 for(k=1;k<=Nz;k++)
		 {
		 diver[i][j][k] = (VxNew[i+1][j][k]-VxNew[i-1][j][k])/(2*dx)+
				 (VyNew[i][j+1][k]-VyNew[i][j-1][k])/(2*dy)+
				 (VzNew[i][j][k+1]-VzNew[i][j][k-1])/(2*dz);
		 }

/*------------ Step III (counting pressure)---------------------------*/

	for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
			for(k=1;k<=Nz;k++)
		  pNew[i][j][k] = pOld[i][j][k] - dt*diver[i][j][k]/gamma;

/*------------ Step IV (pressure to velocities)-----------------------*/

			  //free conditions on horizontal surfaces
	for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
		{
	  pNew[i][j][0] = pNew[i][j][1];
	  pNew[i][j][Nz+1] = pNew[i][j][Nz];
		}

				//periodic conditions on vertical surfaces
	for(i=1;i<=Nx;i++)
	  for(k=0;k<=Nz+1;k++)
		{
	  pNew[i][0][k] = pNew[i][Ny][k];
	  pNew[i][Ny+1][k] = pNew[i][1][k];
		}

				//gradient-periodic conditions on stream surfaces
	for(j=0;j<=Ny+1;j++)
	  for(k=0;k<=Nz+1;k++)
	 {
	  pNew[0][j][k]= pNew[Nx][j][k]-(p2-p1);
	  pNew[Nx+1][j][k]= pNew[1][j][k]+(p2-p1);
	 }

//correction of velocities by pressure
	for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
			for(k=1;k<=Nz;k++)
			  {
		 VxNew[i][j][k] -= dt*(pNew[i+1][j][k]-pNew[i-1][j][k])/(2*dx);
		 VyNew[i][j][k] -= dt*(pNew[i][j+1][k]-pNew[i][j-1][k])/(2*dy);
		 VzNew[i][j][k] -= dt*(pNew[i][j][k+1]-pNew[i][j][k-1])/(2*dz);
			  }

}  //next_layer

/*=========================================================================*/

void main()
{
 long ns;
 int i,j,k;

/*============ Initial condition ==============*/

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
								//+ 1 - 4./Nz/Nz*(k-1-Nz/2)*(k-1-Nz/2)
                                                                ;
				 vy[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
				 vz[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
				 p[0][i][j][k] = p1+(i-0.5)*(p2-p1)/Nx;
				 nut[i][j][k]= 1./Re;
				 }
//	fill_num_points();

	fv = prepout("vv.dat",0);

	fnu = prepout("nut.dat",0);

//	fsf = prepout("strfunc.dat");

/*=================== Main block ====================================*/

 konez=false;
 int n0 = ns%2, n1 = (ns+1)%2;
	while(tm<tmmax && ! konez)
	{
        n0 = ns%2;  n1 = (ns+1)%2;
	 next_layer(vx[n0],vx[n1],vy[n0],vy[n1],vz[n0],vz[n1],p[n0],p[n1]);
    nut_by_flux(n0);
	 if((int)((tm+dt/2)/STEP)-(int)((tm-dt/2)/STEP))
	 {
	 printing(ns);
	 }

	  ns++;
	  tm+=dt;
	}

	printf("\n The End \n");
	putch(getch());
	printf("\n vx[%u][%u][k]:\n",Nx/2,Ny/2);
	for(k=1;k<=Nz;k++)
		printf("%g\n",vx[n0][Nx/2][Ny/2][k]);
	putch(getch());
	fclose(fv);
	fclose(fsf);
}//main
