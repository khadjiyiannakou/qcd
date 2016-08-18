#include <stdio.h>

typedef struct{
  double r;
  double i;
}complex;

int main(){
  int dummy;
  complex values_part[21][4];
  FILE *pt_in,*pt_out;
  pt_in = fopen("uu3.dat","r");
  pt_out=fopen("uur3.dat","w");
  for(int itime=0;itime<21;itime++)
    for(int itype=0;itype<4;itype++)
      fscanf(pt_in,"%d %d %d %d %lf %lf %d",&dummy,&dummy,&dummy,&dummy,&values_part[itime][itype].r,&values_part[itime][itype].i,&dummy);
  int itype=3;
    for(int itime=0;itime<21;itime++)
      fprintf(pt_out,"%+d %+d %+e %+e\n",itime,itime,values_part[itime][itype].r,values_part[itime][itype].i);
}
