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
#include "../../functions/functions.h"
//#include "functions.h"


int main(){

  char *file, *file1;
  char pathFile[4000];
  char *directory = malloc(5000);
  char *folder;
  int* DEGREE_0;
  struct stat dir = {0};

  struct network net1;

  /********************************************************************
  *********************************************************************

  PHASE I: Network file characerization.

    This part of the code is intended to read the source file
    containing the network links in a .txt file with 2 columns,
    each row representing the link and both columns the involved
    nodes. This part returns several vectors containing this
    information.

  *********************************************************************
  *********************************************************************/

  file = malloc(100*sizeof(char));
  file1 = malloc(100*sizeof(char));
  folder = malloc(100*sizeof(char));

  printf("\n");
  printf("\n");
  printf("    ******************************************************\n");
  printf("                N E T W O R K S  TOPOLOGY.\n");
  printf("\n");
  printf("                                                 June 2017\n");
  printf("                                      Author: Jaume Palmer\n");
  printf("    ******************************************************\n");
  printf("\n");
  printf("\n");
  printf("      1) Enter the network filename :");
  scanf("%s", file1);
  printf("\n");

  // printf("      2) Enter the folder path (./foldername/):");
  // scanf("%s", folder);
  strcpy(folder, "../../nets/");

  printf("\n");
  char* file2 = concat(file1,"_topology(rand)");

  strcpy(directory,file2);
  if(stat(directory, &dir) == -1)
  {
      mkdir(directory, 0755);
      printf("created directory testdir successfully! \n");
  }

  file = strcat(folder, file1);

  net1.EN = get_E_N(file);

  net1.NODES = malloc(sizeof(int)*net1.EN.N);
  net1.POINT.INI = malloc(net1.EN.N*sizeof(int));
  net1.POINT.FIN = malloc(net1.EN.N*sizeof(int));
  net1.DEGREE = malloc(net1.EN.N*sizeof(int));
  DEGREE_0 = malloc(net1.EN.maxNODE*sizeof(int));

  printf("...getting pointers\n");
  net1.loops = 0;
  get_pointers(&net1, file, DEGREE_0);

  net1.lenLINKS = 2*(net1.EN.E-net1.loops);
  net1.LINKS = (int*)malloc(net1.lenLINKS*sizeof(int));
  net1.pLINKS = (int*)malloc(net1.lenLINKS*sizeof(int));

  InitializeNet(&net1,file, directory);


  /********************************************************************
  *********************************************************************

  PHASE II: SIS model implementation (06/03/17)

  PHASE III: 2 layer multiplex model (24/03/17)

  PHASE IV: Competition (19/04/2017)

  *********************************************************************
  *********************************************************************/

  int j,r,rep, filecount, netstep, netsize;
  int RUNS = 1; int REPS = 1;
  float t, tau, rand1;
  float step, delta_step;
  bool keep_running;

  /**********************************************************************
    Critical lambdas:

      - Jazz: 0.025
      - Mail: 0.045
      - Brightkite: 0.01
  ***********************************************************************/

  float epsilon, epsilon0 = 0.5;
  epsilon = epsilon0;
  float coupling, coupling0 = 0.75, MAXcoupling = 0.75, delta_coupling = 0.2; // float MAXcoupling;
  coupling = coupling0;
  float lambda,lambda0 = 0.3, MAXlambda = 0.3, delta_lambda = (MAXlambda - lambda0)/(float)REPS; // rate of activation
  lambda = lambda0;
  float lambdaRel, lambdaC;
  // float sigma = 1; // rate of recovery
  float mu = 1; // rate of spontaneous activation
  // float SumTotalRate;
  float OutputAVG[RUNS];
  float rho_avg;
  struct network* LAYERS;
  struct network* LAYERS_0;
  float* G_COEFFS;

  int cl = 1; //Competing networks

  float FRAC[cl];

  FILE* Evolution_File;
  FILE* EpidemicOutput;

  char *SigmaDir = malloc(1000);
  char *EvolutionDir = malloc(5000);
  sprintf(EvolutionDir, "./%s/Evolution",directory); if(stat(EvolutionDir, &dir) == -1){mkdir(EvolutionDir, 0755);}

  LAYERS_0 = malloc(cl*sizeof(struct network));
  LAYERS = malloc(cl*sizeof(struct network));

  printf("...starting epidemic spreading\n");

  /***********************************************************************

    Nsus = nodes still only offline (Susceptible)
    net1.SIS.Eact = online nodes susceptible to be activated (Passive)
    net1.SIS.Nact = Active nodes that can become passive (Active)

    Es2a = offline nodes that have neighbours online and can be activated

  ***********************************************************************/

  /***********************************************************************

  In NS:

    '1' = node already in the online (Active or Passive)
    '0' = node susceptible

  ***********************************************************************/

  sprintf(pathFile, "./%s/epidemic_output_avg.txt", directory);
  EpidemicOutput = fopen(pathFile,"w");

  printf("...running simulation\n");


  LAYERS_0[0] = net1;
  LAYERS_0[0].SIS.Nsus = net1.EN.N;

  LAYERS_0[0].SIS.NI = malloc(net1.EN.N*sizeof(int)); // 1D array of infected nodes
  LAYERS_0[0].SIS.pointNS = malloc(sizeof(int)*net1.EN.N);
  LAYERS_0[0].SIS.NS = malloc(sizeof(int)*net1.EN.N);
  LAYERS_0[0].SIS.pointEA_S = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS_0[0].SIS.pointEA_P = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS_0[0].SIS.pointEA_SP = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS_0[0].SIS.ppS = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS_0[0].SIS.ppP = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS_0[0].SIS.EA2pnt = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS_0[0].SIS.pointerEA = malloc(net1.lenLINKS*sizeof(int)); //same size as LINKS
  LAYERS_0[0].SIS.EA = malloc(2*sizeof(int*)); // 2D array of active links
  for(j=0;j<2;j++){
    LAYERS_0[0].SIS.EA[j] = malloc(net1.EN.E*sizeof(int));
  }


  //------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------

  LAYERS[0].SIS.Nsus = net1.EN.N;

  LAYERS[0].NODES = malloc(net1.EN.N*sizeof(int));
  LAYERS[0].LINKS = malloc(net1.lenLINKS*sizeof(int));
  LAYERS[0].DEGREE = malloc(net1.EN.N*sizeof(int));
  LAYERS[0].POINT.INI = malloc(net1.EN.N*sizeof(int));
  LAYERS[0].POINT.FIN = malloc(net1.EN.N*sizeof(int));

  LAYERS[0].SIS.NI = malloc(net1.EN.N*sizeof(int)); // 1D array of infected nodes
  LAYERS[0].SIS.pointNS = malloc(sizeof(int)*net1.EN.N);
  LAYERS[0].SIS.NS = malloc(sizeof(int)*net1.EN.N);
  LAYERS[0].SIS.pointEA_S = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS[0].SIS.pointEA_P = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS[0].SIS.pointEA_SP = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS[0].SIS.ppS = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS[0].SIS.ppP = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS[0].SIS.EA2pnt = malloc(sizeof(int)*net1.lenLINKS);
  LAYERS[0].SIS.pointerEA = malloc(net1.lenLINKS*sizeof(int)); //same size as LINKS
  LAYERS[0].SIS.EA = malloc(2*sizeof(int*)); // 2D array of active links
  for(j=0;j<2;j++){
    LAYERS[0].SIS.EA[j] = malloc(net1.EN.E*sizeof(int));
  }


  //------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------

  G_COEFFS = malloc((4*cl)*sizeof(float));

  // At the beginning there is only one OSN, so initial population is built
  //  with only this

  printf("\n");
  printf("      %d) Enter the initial fraction of active users in net %d: ",2,1);
  scanf("%f", &FRAC[0]);
  printf("\n");

  initial_population(&LAYERS_0[0], FRAC[0]);


  if(strcmp(file1,"mail_gcc") == 0){
    lambdaC = 0.045; // Known value for previous study
  }else if(strcmp(file1,"jazz") == 0){
    lambdaC = 0.025;
  }else if(strcmp(file1,"bkite_gcc") == 0){
    lambdaC = 0.01;
  }else{
    lambdaC = lambda0;
  }

  lambdaRel = lambda/lambdaC;

  printf("      -------------------------------\n");
  printf("      Parameters:\n");
  printf("\n");
  printf("        - lambdaRel = %f\n", lambdaRel);
  printf("        - v = %f\n", lambda/mu);
  printf("      -------------------------------\n");
  printf("\n");
  printf("...>> spreading infection                             \n");


  clock_t begin = clock(); srand(time(NULL));

  epsilon = epsilon0; coupling = coupling0;

  for(rep=0;rep<REPS;rep++){ ///////////////////////////////////////////////////

    sprintf(SigmaDir, "./%s/sigma_%f",EvolutionDir,coupling); if(stat(SigmaDir, &dir) == -1){mkdir(SigmaDir, 0755);}

    for(r=1;r<=RUNS;r++){ //////////////////////////////////////////////////////

      epsilon = epsilon0;
      sprintf(pathFile, "./%s/%d.txt",SigmaDir,r);
      Evolution_File = fopen(pathFile,"w");
      WriteHeader(Evolution_File,lambdaRel,epsilon,coupling,FRAC[0],FRAC[1]);

      ResetVectors(LAYERS,LAYERS_0,cl);
      Compute_Density(LAYERS, cl);
      Update_Weight(LAYERS, cl, epsilon0, coupling);
      Update_ActRate(LAYERS, cl, lambda, mu);

      filecount = 0;
      step = 0;
      delta_step = 2;

      t = 0;

      netsize = LAYERS[0].EN.N - LAYERS[0].SIS.Nsus;
      netstep = 200;
      // netstep = LAYERS[0].EN.N/300;
      WriteNet_Act(LAYERS[0],EvolutionDir,filecount);
      filecount += 1;

      keep_running = true;

      // GILLESPIE METHOD //////////////////////////////////////////////////////
      //    "keep_running" condition is added to be able to
      //    start Gillespie without any active node
      //    and give enough time for the infection to grow

      while(keep_running){

        rand1 = ((double)rand()+1)/ (double)((unsigned)RAND_MAX);
        // rand2 = ((double)rand()+1)/ (double)((unsigned)RAND_MAX);

        //a_react  | a_deact  | a_mmact  | a_act

        tau = (1/(float)(mu))*log(1/rand1);

        RandomActivation(LAYERS);
        Compute_Density(LAYERS, cl);
        Update_Weight(LAYERS, cl, epsilon, coupling);
        Update_ActRate(LAYERS, cl, lambda, mu);

        if(LAYERS[0].EN.N - LAYERS[0].SIS.Nsus > netsize + netstep){
          WriteNet_Act(LAYERS[0],EvolutionDir,filecount);
          filecount += 1;
          netsize = LAYERS[0].EN.N - LAYERS[0].SIS.Nsus;
        }

        t += tau;

        if(t>=step){
          step += delta_step;

          fprintf(Evolution_File,"%f %f %f\n",
                  t,
                  LAYERS[0].rho.act,
                  epsilon);
        }

        keep_running = (LAYERS[0].SIS.Nsus/(float)LAYERS[0].EN.N > 0.05);

      }// END OF GILLESPIE /////////////////////////////////////////////////////

    fclose(Evolution_File);

    OutputAVG[r-1] = LAYERS[0].rho.act;

    } // END OF RUNS ///////////////////////////////////////////////////////////

    rho_avg = ComputeMean(OutputAVG,RUNS);
    fprintf(EpidemicOutput, "%f %f\n",lambda,rho_avg);

    if(coupling + delta_coupling <= MAXcoupling){
      coupling += delta_coupling;
    }
    if(lambda + delta_lambda <= MAXlambda){
      lambda += delta_lambda;
    }

  } // END OF REPS /////////////////////////////////////////////////////////////

  fclose(EpidemicOutput);
  free(G_COEFFS);

  if(EpidemicDied(LAYERS,cl)){ printf("Epidemic died.                                           \n");}

  printf("End of SIS run.                                                   \n");
  printf("\n");

  printf("...freeing\n");
  free(DEGREE_0);

  printf("...freeing 1 ");free(LAYERS_0[0].SIS.pointerEA);
  printf(" 2");free(LAYERS_0[0].SIS.pointEA_P);
  printf(" 3");free(LAYERS_0[0].SIS.pointEA_S);
  printf(" 4");free(LAYERS_0[0].SIS.NI);
  printf(" 5");free(LAYERS_0[0].SIS.pointNS);
  printf(" 6");free(LAYERS_0[0].SIS.NS);free(LAYERS_0[0].SIS.EA2pnt);
  printf(" 7");

  for(j=0;j<2;j++){
    free(LAYERS_0[0].SIS.EA[j]);
  }free(LAYERS_0[0].SIS.EA);

  printf(" 8");free(LAYERS[0].SIS.pointerEA);
  printf(" 9");free(LAYERS[0].SIS.pointEA_P);
  printf(" 10");//free(LAYERS[0].SIS.pointEA_S);
  printf(" 11");free(LAYERS[0].SIS.pointEA_SP);
  printf(" 12");free(LAYERS[0].POINT.INI);
  printf(" 13");free(LAYERS[0].POINT.FIN);
  printf(" 14");free(LAYERS[0].LINKS);
  printf(" 15");free(LAYERS[0].DEGREE);
  printf(" 16");free(LAYERS[0].NODES);
  printf(" 17");free(LAYERS[0].SIS.NI);
  printf(" 18");//free(LAYERS[0].SIS.pointNS);
  printf(" 19\n");//free(LAYERS[0].SIS.NS);
  printf("...freeing layer %d finished!\n", 1);
  for(j=0;j<2;j++){
    free(LAYERS[0].SIS.EA[j]);
  }free(LAYERS[0].SIS.EA);


  // Memory leaks !!
  //free(folder);
  //free(file);
  //free(file1);
  //free(directory1);
  //free(directory2);
  //free(directory3);
  //free(directory4);


  clock_t end = clock();
  double time_spent = (double)(end-begin)/CLOCKS_PER_SEC;

  printf("(time spent: %f s)\n", time_spent);
  printf("\n");

  return 0;
}
