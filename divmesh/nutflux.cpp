# include <math.h>
double vx[2][Nx+2][Ny+2][Nz+2],
       vy[2][Nx+2][Ny+2][Nz+2],
       vz[2][Nx+2][Ny+2][Nz+2];

extern double nut[Nx+2][Ny+2][Nz+2];
# define maschtab 10
# define modul(a,b,c) (sqrt((a)*(a)+(b)*(b)+(c)*(c)))

void nut_by_flux() //calculating nu_turbulent by velocity fluctuations
{
double flux;
int i,j,k;
for(i=1;i<=Nx;i++)
    for(j=1;j<=Ny;j++)
        for(k=1;k<=Nz;k++)
            {
            flux = 0;
            flux += modul(vx[i][j][k]-vx[i-1][j][k],vy[i][j][k]-vy[i-1][j][k],vz[i][j][k]-vz[i-1][j][k];
            flux += modul(vx[i][j][k]-vx[i+1][j][k],vy[i][j][k]-vy[i+1][j][k],vz[i][j][k]-vz[i+1][j][k]);
            flux += modul(vx[i][j][k]-vx[i][j-1][k],vy[i][j][k]-vy[i][j-1][k],vz[i][j][k]-vz[i][j-1][k];
            flux += modul(vx[i][j][k]-vx[i][j+1][k],vy[i][j][k]-vy[i][j+1][k],vz[i][j][k]-vz[i][j+1][k]);
            flux += modul(vx[i][j][k]-vx[i][j][k-1],vy[i][j][k]-vy[i][j][k-1],vz[i][j][k]-vz[i][j][k-1];
            flux += modul(vx[i][j][k]-vx[i][j][k+1],vy[i][j][k]-vy[i][j][k+1],vz[i][j][k]-vz[i][j][k+1]);
            nut[i][j][k] = maschtab * flux;
            }
}

