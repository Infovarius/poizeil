//approximate solving of transfer equation;
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <conio.h>

# define pi 3.141592653579893
# define n 100
# define tmmax 20          //number of meshing points
# define STEP 0.1
# define EACH 10
# define uplevel 50
# define Re 10
# define Noise 0.01
double xb=0,xe=1,         //line boundaries
       dt=0.0000001;        //time increasement
double pgrad=.3;          //pressure gradient
double u[2][n+1],         //function meanings
       nu[n+1];
FILE *fil;
double alpha,h,tm;
long i,k,i1,i2;
bool konez;

FILE *openfile(const char *st)
{
FILE *f;
if((f = fopen(st,"w+"))==NULL)
		{
      printf("File %s can't be open",st);
      exit(-1);
      }
  return(f);
}

//initial conditions
double f(double x)
{
return(5-4*fabs(x-(xb+xe)/2.));
}

//boundary conditions
double psi1(double t)
{
return(0);
}

double psi2(double t)
{
return(0);
}

double check_stop()
{
double eps = 0;
for(int i=0;i<=n;i++)
    eps = max(fabs(u[0][i]-u[1][i]), eps);
if(eps>uplevel) konez = true;
return(eps);
}

void vivod(int ind, double t) //printing of meanings on each layer
{
long i;
clrscr();
printf("at t=%0.5f:\n",t);
fprintf(fil,"{");
for(i=0;i<=n;i++)
    {if(!(i%EACH)) printf("u[%0.3f]=%8.6g\n",i*h+xb,u[ind][i]);
     fprintf(fil,"%0.6f",u[ind][i]);
     if(i==n) fprintf(fil,"}\n");
        else fprintf(fil,",");
    }
printf("EPS = %f",check_stop());
if(kbhit()&&getch()=='q') konez = true;
}

void periodic(double *u)              //periodic boundary conditions
{
u[0]=u[n-1];
u[n]=u[1];
}

void free(double *u)                //free gradient on border
{
u[0]=u[1];
u[n]=u[n-1];
}

void harsh(double *u, double left, double right) //hard conditions
{
u[0] = 2*left - u[1];
u[n] = 2*right - u[n-1];
}

void stream(double *u,double left)         //some condition
{
u[0] = 2*left - u[1];
u[n] = u[n-1];
}

void main()
{
tm = 0;
k = 2;
h=(xe-xb)/n;  //steps
alpha=dt/h/h;              	  //convergence parameter
for(i=0;i<=n+1;i++) nu[i] = 1./Re;
for(i=0;i<=n;i++) u[0][i]=f(xb+h*i) + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
                                  //first layer(from I init.condition)
fil = openfile("res");
fprintf(fil,"%f\n",alpha);
vivod(0,0);
harsh(u[0],3,3);
konez=false;

while(tm<tmmax && ! konez)
    {
    i1 = k % 2;
    i2 = (k+1) % 2;
    for(i=1;i<=n-1;i++)
    	u[i2][i] = u[i1][i]+dt*(   pgrad-u[i1][i]*(u[i1][i+1]-u[i1][i-1])/2/h +
                                (nu[i+1]-nu[i-1])*(u[i1][i+1]-u[i1][i-1])/4/h/h+
                                nu[i]*(u[i1][i+1]-2*u[i1][i]+u[i1][i-1])/h/h
                                );
    harsh(u[i2],3,3);
    check_stop();
if((int)((tm+dt/2)/STEP)-(int)((tm-dt/2)/STEP))
    vivod(i2,tm);
    tm+=dt;
    k++;
    };
vivod(i2,tm);
fclose(fil);
printf("\nwork is done");
while(kbhit())getch();
getch();
}
