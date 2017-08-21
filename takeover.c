#include<stdio.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<math.h>
#include<python2.7/Python.h>
#include<stdbool.h>
// #include "functions65.h"

struct competition{
  float lambda;
  float sigma;
  float epsilon;
  float coex;
  float dom;
  float tkv;
  float sur;
  float rho1_0;
  float rho2_0;
  float Time;
  float gap;
};

struct competition ReadStabilityFile(char *filename){

  char str[100];
  char *line, *token;
  float rho1, rho2, rho1_0, rho2_0, time_ss, SumTime;
  int domCases = 0, coexCases = 0, tkvCases = 0, surCases = 0;
  int count;
  float deltaRho, SumDeltaRho;
  struct competition Prob;

  Prob.gap = SumDeltaRho = 0;

  filename=strcat(filename,".txt");
  count = 0;
  SumTime = 0.0;
  FILE *fp;
  fp = fopen(filename, "r");
  // printf("%s\n", filename );
  if (fp){
    line=fgets(str,100,fp);
    while (line != NULL){

      token = strtok(line," "); // Reading the first token

      // BEGIN READ HEADER
      while(strcmp(token,"#") == 0){ // If first token is # skip

        token = strtok(NULL," "); //reading following tokens

        if(strcmp(token,"lambda") == 0){
          token = strtok(NULL," "); Prob.lambda = atof(token);
        }else if(strcmp(token,"epsilon0") == 0){
          token = strtok(NULL," "); Prob.epsilon = atof(token);
        }else if(strcmp(token,"sigma") == 0){
          token = strtok(NULL," "); Prob.sigma = atof(token);
        }else if(strcmp(token,"rho1_0") == 0){
          token = strtok(NULL," ");
          rho1_0 = atof(token); Prob.rho1_0 = rho1_0;
        }else if(strcmp(token,"rho2_0") == 0){
          token = strtok(NULL," ");
          rho2_0 = atof(token); Prob.rho2_0 = rho2_0;
        }

        line = fgets(str,100,fp); //Read next line
        token = strtok(line," ");
      }
      // END READ HEADER
      rho1 = atof(token); // Parsing the char pointer to int
      token = strtok(NULL," "); //reading following tokens
      rho2 = atof(token);  // Parsing the char pointer to int

      // In a file, all repetitins will have the same rho1_0, rho2_0

      if(rho1_0 < rho2_0){

        // Network 1 is the fittest, starting with a close to 0 density
        //  A network is considered extinct when it ends with a density lower than
        //  the fittest initial conditions. Then the surviving network dominates
        //
        //  If the fittest ends up with a density higher than its initial, but lower
        //  than the early network it has survived, despite not taken over

        if(rho1 > rho1_0){
          surCases += 1;
          if(rho1 > rho2){ // Takeover
            tkvCases += 1;
            if(rho2 < rho1_0){ // Domination
              domCases += 1;
            }else{
              coexCases += 1; // Coexistence
              SumDeltaRho += rho1 - rho2;
            }
          }
        }

      }else{
        // Network 2 is the fittest, starting with a close to 0 density
        if(rho2 > rho2_0){
          surCases += 1;
          if(rho2 > rho1){ // Takeover
            tkvCases += 1;
            if(rho1 < rho2_0){ // Domination
              domCases += 1;
            }else{
              coexCases += 1; // Coexistence
            }
          }
        }
      }

      token = strtok(NULL," "); //reading following tokens
      time_ss = atof(token);  // Parsing the char pointer to int
      SumTime += time_ss;
      count += 1;
      line = fgets(str,100,fp); //Read next line
    }
  }
  else{printf( "Error!!! file not found !!\n");}


  Prob.dom = domCases/(float)count;
  Prob.coex = coexCases/(float)count;
  Prob.tkv = tkvCases/(float)count;
  Prob.sur = surCases/(float)count;
  Prob.Time = SumTime/(float)count;

  deltaRho = SumDeltaRho/(float)coexCases;

  Prob.gap = deltaRho;

  return Prob;

}

int main(){


  // 1) Open file of evolution
  // 2) store initial densities
  // 3) compare if take over
  char* fileTest = malloc(100*sizeof(char));
  char* folder = malloc(1000);

  struct competition Prob;
  int file, NumFiles;
  FILE *TkoverFile;

  printf("\n");
  printf("      - Enter the network folder name: ");
  scanf("%s", folder);
  printf("\n");
  TkoverFile = fopen("./TkOProbability.txt","w");
  NumFiles = 26; //number of files inside output folder

  for(file=0;file<NumFiles;file++){

    sprintf(fileTest,"./%s/output/final_output_%d",folder,file);
    Prob = ReadStabilityFile(fileTest);
    fprintf(TkoverFile,"%f %f %f %f %f %f %f\n", Prob.epsilon, Prob.sur, Prob.tkv, Prob.coex, Prob.dom, Prob.Time, Prob.gap);

  // printf("  --------------------------------------\n");
  // printf("  - Probability of domination: %f\n", Prob.dom);
  // printf("  - Probability of Coexistence: %f\n", Prob.coex);
  // printf("  - Probability of Take Over: %f\n", Prob.tkv);
  // printf("  --------------------------------------\n");
  // printf("  - Lambda: %f\n", Prob.lambda);
  // printf("  - Epsilon: %f\n", Prob.epsilon);
  // printf("  - Sigma: %f\n", Prob.sigma);
  // printf("  --------------------------------------\n");
  // printf("  - Rho1_0: %f\n", Prob.rho1_0);
  // printf("  - Rho2_0: %f\n", Prob.rho2_0);
  // printf("  --------------------------------------\n");
  }

  fclose(TkoverFile);
}
