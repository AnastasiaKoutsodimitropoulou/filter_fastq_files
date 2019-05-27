#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#define MAXLINE 450
#define WinLen 20
#define WinThres 20
#define min(X, Y) (((X) < (Y)) ? (X) : (Y))

int countlines(char *filename)
{
  // count the number of lines in the file
  FILE *fp = fopen(filename,"r");
  int ch=0;
  int count=0;

  do{
    ch = fgetc(fp);
    if(ch == '\n')
      count++;
  }while(ch!=EOF);
  fclose(fp);
  return count;
}

int main(int argc,char **argv)
{

  // Open the file given from CLI for input
  FILE * Fin= fopen(argv[1], "r");
  // Open the file given from CLI for output
  FILE * Fout= fopen(argv[2], "w");

  int i,j;
  int m=0;
  int Line, lines, dest;
  int rank, size, dummy;
  int tag=55;
  lines = countlines(argv[1])-1;
  printf("%d\n", lines);

  //initialize mpi
  MPI_Init(&argc, &argv);
  //get current process id
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //get number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;
  char ** buffer;
  buffer=(char**)malloc(sizeof(char*)*lines);
  for(i=0;i<lines;i++)
    buffer[i]=(char*)malloc(sizeof(char)*MAXLINE);
  int * buf;
  buf=(int*)malloc(sizeof(int)*(lines/4));
  size_t len = 0;

  for(Line=0;Line<lines;Line++)
    getline(&buffer[Line], &len, Fin);
 
  //master process sends data to all slave processes
  if(rank == 0){
    for(i=1; i<size; i++){
       MPI_Send(&m, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
    }
    printf("Sent\n");
  }else if(rank!=0){
  //all slaves receive message from master
  for(i=1; i<size; i++){
   MPI_Recv(&m, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
   //printf("Received m:%d\n", m);
 
  
  for(j=0;j<((lines/4)/size);j++){
  // The number of nucleotides in the second line
  // or equally in the last line
  int MaxLen=strlen(buffer[(j*4)+1])-1;
  printf("Number of Nucelotides %d:\n", MaxLen);

  // length of line[1] and line[3] MUST be equally
  if (strlen(buffer[(j*4)+3])!=strlen(buffer[(j*4)+1]))
  {
    if(strlen(buffer[(j*4)+1])<strlen(buffer[(j*4)+3])){
        buffer[(j*4)+3][strlen(buffer[(j*4)+1])]="\0";
    }else{
        buffer[(j*4)+1][strlen(buffer[(j*4)+3])]="\0";
    }

  }


  float Qual=0;

  //for (i=0;i<MaxLen;i++)
  //printf("%c %d\n",line[3][i],line[3][i]-33);

  // start and end position of the sliding window
  int start=0;
  int end=start+WinLen;

  Qual=WinThres+1;
  // slide the window while:
  // the end position has not reached the end of the line
  // and the mean quality score is above the minimum threshold
  while ((end<=MaxLen)&&Qual>WinThres)
  {
    //  printf("%c %d\n",line[3][start],line[3][start]-33);
    // printf("Window at position (%d,%d)\n",start,end-1);

    // calculate the mean quality score
    Qual=0;
    for (int k=start;k<end;k++)
      Qual+=buffer[(j*4)+3][k]-33;
    //printf("Accumulated Qual=%f \t Mean Qual=%f \n",Qual,Qual/WinLen);
    Qual/=WinLen;

    //slide the window by one position to the right
    start++;
    end=start+WinLen;
  }

  start--;
  printf("Nucelotides after position %d have mean window quality under %d\n",start,WinThres);

  // trim out the filter positions from
  // the second and the last lines up to
  strncpy(buffer[(j*4)+1],buffer[(j*4)+1],start);
  buffer[(j*4)+1][start]='\0';
  buffer[(j*4)+3][start]='\0';

  //write the filtered fastq to the output file
  fprintf(Fout,"%s",buffer[(j*4)] );
  fprintf(Fout,"%s\n",buffer[((j*4))+1] );
  fprintf(Fout,"%s",buffer[((j*4))+2] );
  fprintf(Fout,"%s\n",buffer[((j*4))+3] );
 
  printf("%d\n",m);
  }
  }
  }
  MPI_Finalize();
  //close the files opened
  fclose(Fin);
  fclose(Fout);
  printf("Files closed\n");

  // free the allocated memory
  for (i=0;i<lines;i++) {
    free(buffer[i]);
  }
  printf("Freed\n");
  free(buffer);
  printf("All good\n");
  
  exit(0);

}
