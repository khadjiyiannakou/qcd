#include <stdio.h>
#include <stdlib.h>
void qcd_swap_8(char *R, int N)
{
   register char *i,*j,*k;
   char swap;
   char *max;

   max = R+(N<<3);
   for(i=R;i<max;i+=8)
   {
      j=i; k=j+7;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
   }
}

int main () {
  FILE * pFile;
  long lSize;
  char * buffer;
  size_t result;
int i;
  double* d_buffer;
  d_buffer=(double*)buffer;
  pFile = fopen ( "kale.dat" , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);
  printf("%d\n",lSize);
  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
  
  /* the whole file is now loaded in the memory buffer. */
qcd_swap_8(buffer,2*4*3*16*8*8*8);
  // terminate
//printf("%c \n",buffer[1]);
for(i==0;i<2*4*3*16*8*8*8;i++){
 printf("%c \n",buffer[i]);
fflush(stdout);
}
  fclose (pFile);
  free (buffer);
  return 0;
}
