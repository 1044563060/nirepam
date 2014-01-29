#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// mmap
//#include <fcntl.h>
//#include <sys/mman.h>

//
// in R
// library(cluster)
// data = read.table("data0.txt", sep=",")
// m = data[,2:ncol(data)]
// pam(m, k=2, metric="euclidean")
// pam(m, k=2, keep.diss=FALSE, trace.lev=2, metric="euclidean")
// plot(m, col=c(3,1,3,3,1))
//
// pam(distmatrix, diss=TRUE, k=10, keep.diss=FALSE, trace.lev=2)
//


// constant variables
float  floatMAX = (float)2147483647.0; //2^31
double doubleMAX = 2147483647.0; //2^31
int num_OF_SEQ;
int seq_LEN;
int clusterK;


// distances
double euclidean(int x[seq_LEN], int lx, int y[seq_LEN], int ly);
double levenshtein(int *word1, int len1, int *word2, int len2);
// other functions
int ismedoidf(int *medoids, int caseindex);
int nearestmedoid(int caseindex, int *medoids, FILE *fp, int bannedcase, float *dist);
void findnearests(int *medoids, float *nearest1, float *nearest2, FILE *fp);
void findnearests1(int *medoids, float *nearest1, FILE *fp);
double objectivefunc(int *medoids, FILE *fp);
float getDM(int i, int j, FILE *fp);
void getrowDM(int i, FILE *fp, float row[num_OF_SEQ]);


int main(int argc, char *argv[]){
  // PARAMETERS
  // Paremeter control
  if(argc<7){
    printf("Parametro kopuru okerra: ./nirepam.out 1241 50 10 train_000_1.txt dm.CData 1\n");
    exit(1);
  }
  // input parameters
  // pogram name: argv[0]
  num_OF_SEQ = atoi(argv[1]);
  seq_LEN = atoi(argv[2]);
  clusterK = atoi(argv[3]);
  char *dbname = argv[4];
  char *distma = argv[5];
  int computeDM = atoi(argv[6]);

  // other parameters
  clock_t time1, time2;
  double elapsed;
  int i, j, k, h;



  // pogram parameters
      // malloc matrix
      // int **mat = (int **)malloc(rows * sizeof(int*));
      // for(int i = 0; i < rows; i++) mat[i] = (int *)malloc(cols * sizeof(int));
      // int medoids[clusterK];
  int *medoids = (int *)malloc(clusterK * sizeof(int));


 if(computeDM==1){
  // int sequences[num_OF_SEQ][seq_LEN];
  int (*sequences)[seq_LEN] = (int (*)[seq_LEN])malloc(num_OF_SEQ * sizeof(int*) * seq_LEN);
  // int sequencesLen[num_OF_SEQ];
  int *sequencesLen = (int *)malloc(num_OF_SEQ * sizeof(int));

  // READ FILE: Read the sequences file &
  // initialize the sequences matrix
  printf("READ FILE\n"); fflush(stdout);
  time1 = time(0);
  FILE *fp;
  int x, ci;
  char numStr[100];
  int isID = 1;
  fp = fopen(dbname, "r");
  ci = 0;
  i = 0;
  j = 0;
  while( (x=fgetc(fp)) != EOF ){
    if(x=='0' || x=='1' || x=='2' || x=='3' || 
       x=='4' || x=='5' || x=='6' || x=='7' ||
       x=='8' || x=='9' || x=='.'){
      numStr[ci] = x;
      ci++;
    } else {
      if(x==',' || x=='\n'){
        numStr[ci] = '\0';
        ci = 0;
        if(isID==0){
          sequences[i][j] = atoi(numStr);
          j++;
        } else{ // the 1st element is the sequence ID
          isID=0;
        }
      }
      if(x=='\n'){
        // save the sequences length
        sequencesLen[i] = j;
        // fill the remain positions with -1
        for(;j<seq_LEN;j++){sequences[i][j] = -1;}
        // update the indexes
        isID = 1;
        j = 0;
        i++;
      }
    }
  }
  fclose(fp);
  time2 = time(0);
  elapsed = difftime(time2, time1);
  printf("  (%d sec.) reading the file.\n", (int)elapsed); fflush(stdout);


  // create the DISTANCE MATRIX
  // create the distance matrix (in the hard drive)
  printf("DISTANCE MATRIX\n"); fflush(stdout);
  //float (*dm)[num_OF_SEQ] = (float (*)[num_OF_SEQ])malloc(num_OF_SEQ * sizeof(float*) * num_OF_SEQ);
  time1 = time(0);
  // empty the file
  FILE *fp2 = fopen(distma, "wb");
  fclose(fp2);
  // open the file in append mode
  fp2 = fopen(distma, "ab");
  float dmW1[num_OF_SEQ];
  for(i=0;i<num_OF_SEQ;i++){
    for(j=0;j<num_OF_SEQ;j++){
      //dm[i][j] = (float)levenshtein(sequences[i],sequencesLen[i], sequences[j],sequencesLen[j]);
      dmW1[j] = (float)levenshtein(sequences[i],sequencesLen[i], sequences[j],sequencesLen[j]);
      //dmW1[j] = (float)euclidean(sequences[i],sequencesLen[i], sequences[j],sequencesLen[j]);
      //printf("%f,", dmW1[j]);
    }
    fwrite(dmW1, sizeof(float), num_OF_SEQ, fp2);
    //printf("\n");
  }
  fclose(fp2);
  free(sequences);
  free(sequencesLen);
  time2 = time(0);
  elapsed = difftime(time2, time1);
  printf("  (%d sec.) creating the distance matrix.\n", (int)elapsed); fflush(stdout);
 }


  // open the file descriptor
  FILE *fp3 = fopen(distma, "rb");


  //// BUILD phase ////
  printf("BUILD phase\n"); fflush(stdout);
  time1 = time(0);
  // initialize medoid array with -1: medoids collection is empty
  for(i=0;i<clusterK;i++){medoids[i]=-1;}

  // given medoids
  int nmed = 0;
  for(i=7;i<argc;i++){
    medoids[nmed] = atoi(argv[i]);
    nmed++;
  }

 if(nmed==0){
  // select the medoids of the database: medoids[0]
  float rowi[num_OF_SEQ];
  double minsumdist = floatMAX;
  for(i=0;i<num_OF_SEQ;i++){
    double rowsum = 0.0;
    getrowDM(i, fp3, rowi);
    for(j=0;j<num_OF_SEQ;j++){
      //rowsum = rowsum + (double)dm[i][j];
      //rowsum = rowsum + getDM(i,j,fp3);
      rowsum = rowsum + rowi[j];
    }
    if(minsumdist>=rowsum){
      minsumdist = rowsum;
      medoids[0] = i;
    }
  }
  time2 = time(0);
  elapsed = difftime(time2, time1);
  printf("  (%d sec.) new repr. %d\n", (int)elapsed, medoids[0]); fflush(stdout);
  nmed++;
 }

  // select the other medoids: from medoids[1] to medoids[clusterK-1]
  float nearest1[num_OF_SEQ];
  float rowj[num_OF_SEQ];
  for(i=nmed;i<clusterK;i++){
    time1 = time(0);
    double contributionMAX = -1.0;
    int representant = -1;
    findnearests1(medoids, nearest1, fp3);

    for(j=0;j<num_OF_SEQ;j++){     // i
      // j is a no-medoid case
      int ismedoid = ismedoidf(medoids,j);
      if(ismedoid==1){ continue; }
      // compute the contribution j does
      getrowDM(j, fp3, rowj);
      double contribution = 0.0;
      for(h=0;h<num_OF_SEQ;h++){      // j
        // h is a no-medoid case
        ismedoid = ismedoidf(medoids,h);
        if(ismedoid==1){ continue; }
        // compute the contribution hj does
        //float distmedh;
        //int medoid = nearestmedoid(h, medoids, fp3, -2, &distmedh);
        float distmedh = nearest1[h];
        //float rest = getDM(medoid, h, fp3) - getDM(j, h, fp3);
        float rest = distmedh - rowj[h];
        float chj = rest>0 ? rest : 0.0;
        contribution = contribution + chj;
      }
      if(contributionMAX<=contribution){
        contributionMAX = contribution;
        representant = j;
      }
    }
    medoids[i] = representant;
    time2 = time(0);
    elapsed = difftime(time2, time1);
    printf("  (%d sec.) new repr. %d\n", (int)elapsed, representant); fflush(stdout);
  }
  // print the medoids
  printf("  after build: medoids are %d", medoids[0]);
  for(h=1;h<clusterK;h++){printf(",%d",medoids[h]);}
  printf("\n"); fflush(stdout);
  // print the cost (objective function)
  time1 = time(0);
  double cost1 = objectivefunc(medoids, fp3);
  double cost2;
  time2 = time(0);
  elapsed = difftime(time2, time1);
  printf("  (%d sec.) objective function: %f\n", (int)elapsed, cost1); fflush(stdout);



  //// SWAP phase ////
  printf("SWAP phase\n"); fflush(stdout);

  float rowmedoid[num_OF_SEQ];
  float rowh[num_OF_SEQ];
  float nearest2[num_OF_SEQ];
  int ii;
  while(1){
   time1 = time(0);
   double contributionMIN = doubleMAX;
   int medoidiMIN = -1;
   int medoidhMIN = -1;
   findnearests(medoids, nearest1, nearest2, fp3);

   // for each medoid
   for(i=0;i<clusterK;i++){
    // i is a medoid case
    int medoid = medoids[i];
    getrowDM(medoid, fp3, rowmedoid);

    // for each no-medoid
    for(h=0;h<num_OF_SEQ;h++){
      // h is no medoids case
      int ismedoid = ismedoidf(medoids,h);
      if(ismedoid==1){ continue; }
      getrowDM(h, fp3, rowh);

      // compute the contribution of swaping between medoid[i] and h
      double contribution = 0.0;
      for(j=0;j<num_OF_SEQ;j++){

        // contribution of j
        float contj = (float)0.0;

        /*
        float dmedoidofj;
        int medoidofj = nearestmedoid(j, medoids, fp3, medoid, &dmedoidofj);
        if(dmedoidofj>=rowmedoid[j]){
          float small = dmedoidofj>rowh[j] ? rowh[j] : dmedoidofj;
          contj = small - rowmedoid[j];
        } else {
          if(rowh[j]<dmedoidofj){
            contj = rowh[j] - dmedoidofj;
          }
        }
        */

        if(rowmedoid[j]==nearest1[j]){
          float small = nearest2[j]>rowh[j] ? rowh[j] : nearest2[j];
          contj = small - nearest1[j];
        } else {
          if(rowh[j] < nearest1[j]) {
            contj = rowh[j] - nearest1[j];
          }
        }

        // total sum of contributions
        contribution = contribution + (double)contj;
      }

      // which j has the minimun contribution value?
      if(contributionMIN>contribution){
        contributionMIN = contribution;
        medoidiMIN = i; // medoid index
        medoidhMIN = h; // no-medoid case
      }

    }
   }
  time2 = time(0);
  elapsed = difftime(time2, time1);

   // end criteria
   if(contributionMIN<0){
     cost2 = cost1 + contributionMIN;
     double cost = objectivefunc(medoids, fp3);
     printf("  (%d sec.) swp new  %d <-> %d old; decreasing diss. %f by %f = %f (%f): ", 
       (int)elapsed, medoidhMIN, medoids[medoidiMIN], cost1, contributionMIN, cost2, cost);
     medoids[medoidiMIN] = medoidhMIN; // SWAP
       printf("%d", medoids[0]);
       for(i=1;i<clusterK;i++){printf(",%d",medoids[i]);}
       printf("\n");
     cost1 = cost2;
   } else {
     double cost = objectivefunc(medoids, fp3);
     printf("  (%d sec.) no swp. (%f)\n", (int)elapsed, cost);
     break;
   }
   fflush(stdout);
  }



  // WRITE the results
  float dist;
  int medoid = nearestmedoid(0, medoids, fp3, -2, &dist);
  printf("clustering: %d", medoid);
  for(i=1;i<num_OF_SEQ;i++){
    medoid = nearestmedoid(i, medoids, fp3, -2, &dist);
    printf(",%d", medoid);
  }
  printf("\n"); fflush(stdout);



  // close the file descriptor
  fclose(fp3);
}


float getDM(int i, int j, FILE *fp){
  int ind = (i*num_OF_SEQ+j)*sizeof(float);
  fseek(fp, ind, SEEK_SET);
  float d;
  fread(&d, sizeof(float), 1, fp);
  return d;
}

void getrowDM(int i, FILE *fp, float row[num_OF_SEQ]){
  int ind = (i*num_OF_SEQ)*sizeof(float);
  fseek(fp, ind, SEEK_SET);
  fread(row, sizeof(float), num_OF_SEQ, fp);
}

int ismedoidf(int *medoids, int caseindex){
  int k;
  int ismedoid = 0;
  for(k=0;k<clusterK;k++){
    if(medoids[k]==-1){ break; }
    if(medoids[k]==caseindex){
      ismedoid = 1;
      break;
    }
  }
  return(ismedoid);
}

int nearestmedoid(int caseindex, int *medoids, FILE *fp, int bannedcase, float *pdist){
  int k;
  float mindisth = floatMAX;
  int medoid = -1;
  for(k=0;k<clusterK;k++){
    if(medoids[k]==-1){ break; }
    if(medoids[k]==bannedcase){ continue; }
    if(mindisth > getDM(medoids[k],caseindex,fp)){
      medoid = medoids[k];
      mindisth = getDM(medoid,caseindex,fp);
    }
  }
  *pdist = mindisth;
  return(medoid);
}

void findnearests(int *medoids, float *nearest1, float *nearest2, FILE *fp){
  float mindisth1, mindisth2;
  int i,j;
  for(i=0;i<num_OF_SEQ;i++){
    mindisth1 = floatMAX;
    mindisth2 = floatMAX;
    for(j=0;j<clusterK;j++){
      float dij = getDM(i,medoids[j],fp);
      if(mindisth1 > dij){
        mindisth2 = mindisth1;
        mindisth1 = dij;
      } else {
        if(mindisth2 > dij){
          mindisth2 = dij;
        }
      }
    }
    nearest1[i] = mindisth1;
    nearest2[i] = mindisth2;
  }
}

void findnearests1(int *medoids, float *nearest1, FILE *fp){
  float mindisth1, mindisth2;
  int i,j;
  for(i=0;i<num_OF_SEQ;i++){
    mindisth1 = floatMAX;
    for(j=0;j<clusterK;j++){
      if(medoids[j]==-1){ break; }
      float dij = getDM(i,medoids[j],fp);
      if(mindisth1 > dij){
        mindisth1 = dij;
      }
    }
    nearest1[i] = mindisth1;
  }
}

double objectivefunc(int *medoids, FILE *fp){
  int i,j;
  double cost = 0.0;
  for(i=0;i<num_OF_SEQ;i++){
    float mindist = floatMAX;
    for(j=0;j<clusterK;j++){
      int medoid = medoids[j];
      float dist = getDM(medoid,i,fp);
      if(mindist>dist){
        mindist = dist;
      }
    }
    cost = cost + (float)mindist;
  }
  return(cost);
}


//// DISTANCES ////

double euclidean(int x[seq_LEN], int lx, int y[seq_LEN], int ly){
  int i;
  double sum = 0;
  int len = lx<ly ? lx : ly;
  for(i=0;i<len;i++){
    sum = sum + pow((double)(x[i]-y[i]), (double)2);
  }
  return(sqrt(sum));
}

double levenshtein(int *word1, int len1, int *word2, int len2){
    int matrix[len1 + 1][len2 + 1];
    int i, j;
    for (i = 0; i <= len1; i++) { matrix[i][0] = i; }
    for (i = 0; i <= len2; i++) { matrix[0][i] = i; }
    for (i = 1; i <= len1; i++) {
        int c1 = word1[i-1];
        for (j = 1; j <= len2; j++) {
            int c2 = word2[j-1];
            if (c1 == c2) {
                matrix[i][j] = matrix[i-1][j-1];
            }
            else {
                int delete;
                int insert;
                int substitute;
                int minimum;

                delete = matrix[i-1][j] + 1;
                insert = matrix[i][j-1] + 1;
                substitute = matrix[i-1][j-1] + 1;

                // min3
                minimum = delete;
                if (insert < minimum) {
                    minimum = insert;
                }
                if (substitute < minimum) {
                    minimum = substitute;
                }
                matrix[i][j] = minimum;
            }
        }
    }

    // distance
    int editDist;
    editDist = matrix[len1][len2];
          //printf("%d",word1[0]); for(i=1;i<len1;i++){ printf(",%d", word1[i]); } printf(" eta ");
          //printf("%d",word2[0]); for(i=1;i<len2;i++){ printf(",%d", word2[i]); } printf(" = ");
          //printf("%d\n", editDist);

    // normalization
    int lenmax = (len1>=len2) ? len1 : len2;
    return((double)editDist/(double)lenmax);
}

