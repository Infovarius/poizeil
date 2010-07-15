# define  Nx   30
# define  Ny   30
# define  Nz   30
# define Noise 0.0
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
double dx,dy,dz,lx,ly,lz,
		 vx[2][Nx+2][Ny+2][Nz+2],
		 vy[2][Nx+2][Ny+2][Nz+2],
		 vz[2][Nx+2][Ny+2][Nz+2];

FILE *fsf;


# define norm(a,b,c) ((a)*(a)+(b)*(b)+(c)*(c))
# define min_d min(min(dx,dy),dz)

 double s_func[Nz+2][Nx/2+Ny/2+Nz];
 long num_points[Nx/2+Ny/2+Nz];

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


void struct_func(int q,int ind)//structural function of order q
{
long nn=ind;
int i,j,k,l,m,n,dist;
double d,rx,ry,rz;
for(k=1;k<=Nz;k+=1)
	{
for(i=0;i<=Nx/2+Ny/2+Nz-1;i++) s_func[k][i] = num_points[i] = 0;
	for(i=1;i<=Nx;i++)
	for(j=1;j<=Ny;j++)
			for(l=1;l<=Nx;l++)
				for(m=1;m<=Ny;m++)
					for(n=1;n<=Nz;n++)
						 {
						 if ((i==l)&&(j==m)&&(k==n)) continue;
						 rx = min( abs(l-i),Nx-abs(l-i) )*dx;
						 ry = min( abs(m-j),Ny-abs(m-j) )*dy;
						 rz = abs(n-k)*dz;
						 d = sqrt(rx*rx+ry*ry+rz*rz);
						 d = 2*log(d/min_d)/log(2);
						 //(per wave number add "minus")
						 dist=floor(d);
						 s_func[k][dist]+=pow(
									norm(vx[nn][l][m][n]-vx[nn][i][j][k],
										 vy[nn][l][m][n]-vy[nn][i][j][k],
										 vz[nn][l][m][n]-vz[nn][i][j][k]),
										 q/2.);
						 num_points[dist]++;
						 }
for(i=0;i<=Nx/2+Ny/2+Nz-1&&num_points[i];i++)
	s_func[k][i]=(s_func[k][i]/num_points[i]);
for(i=0;i<=Nx/2+Ny/2+Nz-1&&num_points[i];i++)
	printf("%0.5f;",s_func[k][i]);
printf("\n");
fprintf(fsf,"{%0.6f",s_func[k][0]);
for(i=1;i<=Nx/2+Ny/2+Nz-1&&num_points[i];i++)
	fprintf(fsf,",%0.6f",s_func[k][i]);
fprintf(fsf,"}\n");
	 }//"for" per layers
}//struct_func

void main()
{

fsf = prepout("strfunc.dat");
//fprintf(fsf,"{");

long i,j,k;

  lx = 2; lz = ly= 1; dx = lx/Nx; dy = ly/Ny; dz=lz/Nz;
  //initial filling of arrays
  for(i=0;i<=Nx+1;i++)
	for(j=0;j<=Ny+1;j++)
		for(k=0;k<=Nz+1;k++)
				 {
				 vx[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX
								+ 1 - 4./Nz/Nz*(k-(1+Nz)/2.)*(k-(1+Nz)/2.);
				 vy[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
				 vz[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
				 }
	 struct_func(2,0);
//fprintf(fsf,"}");
fclose(fsf);
 }
