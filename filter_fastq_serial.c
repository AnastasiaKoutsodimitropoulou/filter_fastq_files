#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//change size of sliding window to 20 and mean quality to 19
//MAXLINE has been increased due to problems that occurred while freeing the buffer
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
  int m=0;
  int Line, lines;
  lines = countlines(argv[1]);
  
  // Malloc for a 2-dimensional array of strings with
  // as many lines as the input file and MAXLINE of characters per line
  char ** buffer;
  buffer=(char**)malloc(sizeof(char*)*lines);
  for(i=0;i<lines;i++)
    buffer[i]=(char*)malloc(sizeof(char)*MAXLINE);
  //create one dimemsional array with size = number of records of the input file
  int * buf;
  buf=(int*)malloc(sizeof(int)*(lines/4));
  size_t len = 0;
  
  // read line-by-line all the lines of the file
  // and store each in the array named buffer
  for(Line=0;Line<lines;Line++){
    getline(&buffer[Line], &len, Fin);
  }

  //run for the entirety of the records of the file
  for(j=0;j<lines/4;j++){
  // The number of nucleotides in the second line
  // or equally in the last line
  //m --> first line of the record we are examining
  //m+1-->second line of the record we are examining
  //m+2-->third line of the record we are examining
  //m+3-->fourth line of the record we are examining
  int MaxLen=strlen(buffer[m+1])-1;
  printf("Number of Nucelotides %d:\n", MaxLen);

  // length of line[1] and line[3] MUST be equally
  if (strlen(buffer[m+3])!=strlen(buffer[m+1]))
  {
    //if line[1] and line[3] are not equal make them have the same size
    //the longest line takes the length of the shortest one
    if(strlen(buffer[m+1])<strlen(buffer[m+3])){
        buffer[m+3][strlen(buffer[m+1])]="\0";
    }else{
        buffer[m+1][strlen(buffer[m+3])]="\0";
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
    for (k=start;k<end;k++)
      Qual+=buffer[m+3][k]-33;
    //printf("Accumulated Qual=%f \t Mean Qual=%f \n",Qual,Qual/WinLen);
    Qual/=WinLen;

    //slide the window by one position to the right
    start++;
    end=start+WinLen;
  }
  start--;
  //store start in buf-->used in determining if a record has mean quality below 19
  buf[j]=start;
  printf("Nucelotides after position %d have mean window quality under %d\n",start,WinThres);

  // trim out the filter positions from
  // the second and the last lines up to
  strncpy(buffer[m+1],buffer[m+1],start);
  buffer[m+1][start]='\0';
  strncpy(buffer[m+3],buffer[m+3],start);
  buffer[m+3][start]='\0';

  //write the filtered fastq to the output file
  //check if the quality of the record is below the mean quality
  //if buf[j] <= 0, the filtered record is empty
  //then don't write record j inside the file
  if(buf[j]>0){
  fprintf(Fout,"%s",buffer[m] );
  fprintf(Fout,"%s\n",buffer[m+1] );
  fprintf(Fout,"%s",buffer[m+2] );
  fprintf(Fout,"%s\n",buffer[m+3] );
  }
  //increase m to filter the next record(next four lines)
  m+=4;
  //printf("%d\n", m);

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
