#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <sys/resource.h>

double myclock();
double myclock() {
  static time_t t_start = 0;
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  if( t_start == 0) t_start = ts.tv_sec;
  return (double) (ts.tv_sec - t_start) + ts.tv_nsec * 1.0e-9;
}

int main(int argc, char **argv) {
  int nwords, maxwords = 50;
  int nlines, maxlines = 100;
  int Start, End, Size;
  int ID, nPro;
  int nCores, nNodes;
  int i, k, n, err, *count;

  double nchars = 0;
  double tstart, ttotal;
  FILE *fd;
  char **word, **line;

  MPI_Status status;

  maxwords = atoi(argv[1]);
  maxlines = atoi(argv[2]);
  nCores = atoi(argv[3]);
  nNodes = atoi(argv[4]);
  struct rusage r_usage;

  //Initializes MPI environment.
  err = MPI_Init(&argc, &argv);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &ID);
  err = MPI_Comm_size(MPI_COMM_WORLD, &nPro);

  //Allocates space for ALL the keywords for each thread.
  word = (char**) malloc(maxwords * sizeof(char *));
  count = (int *) malloc(maxwords * sizeof(int));
  for(i = 0; i < maxwords; i++) {
    word[i] = malloc(10);
    count[i] = 0;
  }

  //Allocates space for the wiki data for a given thread.
  line = (char **) malloc(maxlines * sizeof(char*));
  for(i = 0; i < maxlines; i++) {
    line[i] = malloc(2001);
  }

  //Loads keyword data into buffer.
  fd = fopen("/homes/dan/625/keywords.txt", "r");
  nwords = -1;
  do {
    err = fscanf(fd, "%[^\n]\n", word[++nwords]);
  } while(err != EOF && nwords < maxwords);
  fclose(fd);

  //Loads wiki data into buffer.
  fd = fopen("/bulk/daveturner/wiki/wiki.1Mshort", "r");
  if(fd != NULL) {
    nlines = -1;
    do {
      err = fscanf(fd, "%[^\n]\n", line[++nlines]);
      if(line[nlines] != NULL) nchars += (double) strlen(line[nlines]);
    } while(err != EOF && nlines < maxlines);
    fclose(fd);
  }

  //Starts timers.
  if(ID == 0) {
    tstart = myclock(); //Set zero time method from find_keys
    tstart = myclock(); //Start the clock method from find_keys
  }

  //Gets the start and end bounds of the set of keywords that this
  //thread searches through.
  Size = nwords / nPro;
  Start = ID * Size;
  End = Start + Size;
  if(ID == nPro - 1) {
    End = nwords;
  }

  //Compares each keyword to each line. If the keyword is present in the
  //line, the count increases by one.
  for (i = Start; i < End; i++) {
    for (k = 0; k < nlines; k++) {
      if (strstr(line[k], word[i]) != NULL) count[i]++;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int done = 0;
  if(ID == 0) {
    //Master thread prints out the keyword and the number of times it appears
    //in the text to a given output file. Master thread sends data to thread
    //with ID = 1.
    fd = fopen("woutput.out", "w");
    for(k = 0; k < End; k++) {
      fprintf(fd, "%d %s %d\n", k, word[k], count[k]);
    }
    fclose(fd);
    done = 1;
    if(nPro != 1) err = MPI_Send(&done, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
  } else {
    //All other threads append their data to the end of a given file.
    err = MPI_Recv(&done, 1, MPI_INT, ID - 1, 0, MPI_COMM_WORLD, &status);
    fd = fopen("wiki.out", "a");
    for(k = Start; k < End; k++) {
      fprintf(fd, "%d %s %d\n", k, word[k], count[k]);
    }
    fclose(fd);
    if(ID != nPro - 1) {
      err = MPI_Send(&done ,1 , MPI_INT, ID + 1, 0, MPI_COMM_WORLD);
    }
  }
  err = MPI_Barrier(MPI_COMM_WORLD);

  //Prints runtime data to standard out.
  if (ID == 0) {
    ttotal = myclock() - tstart;
    fd = fopen("wiki.out", "a");
    fclose(fd);
    getrusage(RUSAGE_SELF, &r_usage);
    printf("%f, %ld, %d, %d, %d\n", ttotal, r_usage.ru_maxrss, nCores, nNodes, nPro);
  }

  MPI_Finalize();
}
