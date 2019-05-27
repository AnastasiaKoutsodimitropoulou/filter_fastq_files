#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

#define MAXLINE 450
#define WinLen 20
#define WinThres 19
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

  int i,j,k;
  int Line, lines;
  int tid;
  int START,STOP,STEP;
  //take number of theads from CLI
  int THREADS = atoi(argv[1]);

  lines = countlines(argv[1]);
  //set number of threads to be used
  omp_set_num_threads(THREADS);
  //number of records that every thread examines
  STEP=(lines/4)/THREADS;
  
  char ** buffer;
  buffer=(char**)malloc(sizeof(char*)*lines);
  for(i=0;i<lines;i++)
    buffer[i]=(char*)malloc(sizeof(char)*MAXLINE);

  size_t len = 0;
  int *buf;
  buf=(int*)malloc(sizeof(int)*(lines/4));
  //fill buffer serially
  //otherwise threads read the lines out of order 
  for(Line=0;Line<lines;Line++)
    getline(&buffer[Line], &len, Fin);


  //run parallel 
  //private(tid, j, START, STOP)-->every thread must have these variable private to know which records they should filter 
  #pragma omp parallel private(tid, j, START, STOP)
  {
  //get id(class) of each thread 
  tid=omp_get_thread_num();

  //each thread has to know from which part of the buffer they should start filtering and where they should stop
  START=STEP*tid;
  STOP=START+STEP;
  //for each thread that runs parallelly
  for(j=START;j<STOP;j++){
  // The number of nucleotides in the second line
  // or equally in the last line
  //j*4-->every record has 4 lines, first line of record
  //(j*4)+1-->second line of record
  //(j*4)+2->third line of record
  //(j*4)+3-->third line of record
  int MaxLen=strlen(buffer[(j*4)+1])-1;
  printf("Number of Nucelotides %d:\n", MaxLen);

  // length of line[1] and line[3] MUST be equally
  if (strlen(buffer[(j*4)+3])!=strlen(buffer[(j+4)*1]))
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
 // printf("end:%d\n", end);
  Qual=WinThres+1;
  // slide the window while:
  // the end position has not reached the end of the line
  // and the mean quality score is above the minimum threshold
  while ((end<=MaxLen)&&Qual>WinThres)
  {
     // printf("%c %d\n",line[3][start],line[3][start]-33);
     // printf("Window at position (%d,%d)\n",start,end-1);

    // calculate the mean quality score
    Qual=0;
    for (k=start;k<end;k++)
      Qual+=buffer[(j*4)+3][k]-33;
    //printf("Accumulated Qual=%f \t Mean Qual=%f \n",Qual,Qual/WinLen);
    Qual/=WinLen;

    //slide the window by one position to the right
    start++;
    end=start+WinLen;
  }
  start--;
  buf[j]=start;
  printf("Nucelotides after position %d have mean window quality under %d\n",start,WinThres);

  // trim out the filter positions from
  // the second and the last lines up to
  strncpy(buffer[(j*4)+1],buffer[(j*4)+1],start);
  buffer[(j*4)+1][start]='\0';
  strncpy(buffer[(j*4)+3],buffer[(j*4)+3],start);
  buffer[(j*4)+3][start]='\0';

  }
  }
  printf("Ended parallel\n");
  //write in output file seriallly
  //for the number of records check if the quality of the record is above
  //the threshold and then write in output file
  for(j=0; j<lines/4; j++){
     if(buf[j]>0){
 	    fprintf(Fout,"%s",buffer[j*4]);
 	    fprintf(Fout,"%s\n",buffer[(j*4)+1] );
 	    fprintf(Fout,"%s",buffer[(j*4)+2] );
 	    fprintf(Fout,"%s\n",buffer[(j*4)+3] );
     }
  }
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
  free(buf);
  printf("All good\n");

  exit(0);

}
