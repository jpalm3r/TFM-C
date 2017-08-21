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
//#include "functions.h"
#include "../../functions/functions.h"


int main(){

  char *file, *file1;
  int i;
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
  printf("                N E T W O R K S  COEXISTENCE.\n");
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
  char* file2 = concat(file1,"_coex");

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

  int j,r,rep;
  int RUNS = 5; int REPS = 5; int T = 500;
  float t, tau, rand1, rand2;// t_final;
  float step, delta_step;
  bool keep_running, s_state;

  float epsilon,epsilon0,MAXeps; float delta_eps = 0.02;
  epsilon = epsilon0 = 0.5; MAXeps = 1.0;
  float coupling, coupling0 = 0.5, MAXcoupling = 1.0, delta_coupling = 0.1; // float MAXcoupling;
  float SumTotalRate;
  float lambda; // MAXlambda; // rate of activation
  float sigma = 1; // rate of recovery
  float mu = 0.08; // rate of spontaneous activation
  struct network* LAYERS;
  struct network* LAYERS_0;
  //struct competition Prob;
  float* G_COEFFS;
  int lenss; int bkt = 20; int tstep;
  //char* fileTest;

  float MaxEpsAvg;
  float MaxEpsVector[RUNS];

  int cl = 2; //Competing networks

  float FRAC[cl];

  FILE* Evolution_File;
  FILE* SigmaFile;

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


  printf("...running simulation\n");

  for(i=0;i<cl;i++){

    LAYERS_0[i] = net1;
    LAYERS_0[i].SIS.Nsus = net1.EN.N;

    LAYERS_0[i].SIS.NI = malloc(net1.EN.N*sizeof(int)); // 1D array of infected nodes
    LAYERS_0[i].SIS.pointNS = malloc(sizeof(int)*net1.EN.N);
    LAYERS_0[i].SIS.NS = malloc(sizeof(int)*net1.EN.N);
    LAYERS_0[i].SIS.pointEA_S = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS_0[i].SIS.pointEA_P = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS_0[i].SIS.pointEA_SP = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS_0[i].SIS.ppS = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS_0[i].SIS.ppP = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS_0[i].SIS.EA2pnt = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS_0[i].SIS.pointerEA = malloc(net1.lenLINKS*sizeof(int)); //same size as LINKS
    LAYERS_0[i].SIS.EA = malloc(2*sizeof(int*)); // 2D array of active links
    for(j=0;j<2;j++){
      LAYERS_0[i].SIS.EA[j] = malloc(net1.EN.E*sizeof(int));
    }

  }

  //------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------

  for(i=0;i<cl;i++){

    LAYERS[i].SIS.Nsus = net1.EN.N;

    LAYERS[i].NODES = malloc(net1.EN.N*sizeof(int));
    LAYERS[i].LINKS = malloc(net1.lenLINKS*sizeof(int));
    LAYERS[i].pLINKS = malloc(net1.lenLINKS*sizeof(int));
    LAYERS[i].DEGREE = malloc(net1.EN.N*sizeof(int));
    LAYERS[i].POINT.INI = malloc(net1.EN.N*sizeof(int));
    LAYERS[i].POINT.FIN = malloc(net1.EN.N*sizeof(int));

    LAYERS[i].SIS.NI = malloc(net1.EN.N*sizeof(int)); // 1D array of infected nodes
    LAYERS[i].SIS.pointNS = malloc(sizeof(int)*net1.EN.N);
    LAYERS[i].SIS.NS = malloc(sizeof(int)*net1.EN.N);
    LAYERS[i].SIS.pointEA_S = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS[i].SIS.pointEA_P = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS[i].SIS.pointEA_SP = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS[i].SIS.ppS = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS[i].SIS.ppP = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS[i].SIS.EA2pnt = malloc(sizeof(int)*net1.lenLINKS);
    LAYERS[i].SIS.pointerEA = malloc(net1.lenLINKS*sizeof(int)); //same size as LINKS
    LAYERS[i].SIS.EA = malloc(2*sizeof(int*)); // 2D array of active links
    for(j=0;j<2;j++){
      LAYERS[i].SIS.EA[j] = malloc(net1.EN.E*sizeof(int));
    }

    LAYERS[i].SIS.Tss = malloc(sizeof(int)*bkt);
    LAYERS[i].EVOL.real_time = malloc(T*sizeof(float));
    LAYERS[i].EVOL.act_avg = malloc(T*sizeof(float));
    LAYERS[i].EVOL.act = malloc(T*sizeof(float*)); // 2D array of active links
    for(j=0;j<T;j++){
      LAYERS[i].EVOL.act[j] = malloc(RUNS*sizeof(float));
    }
  }

  //------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------

  G_COEFFS = malloc((4*cl)*sizeof(float));


  //MAXlambda = lambda;
  // mu = 0.1;
  // sigma = 1;

  // At the beginning there is only one OSN, so initial population is built
  //  with only this


  if (strcmp(file1,"fbook_gcc") == 0){
    lambda = 0.1;
    FRAC[0] = FRAC[1] = 0.27;
    printf("\n");
    printf("      -------------------------------\n");
    printf("      Network: Facebook\n");
    printf("\n");
    printf("        - rho0_1 = %f\n", FRAC[0]);
    printf("        - rho0_2 = %f\n", FRAC[1]);
    printf("      -------------------------------\n");
    printf("\n");
  }else if (strcmp(file1,"bkite_gcc") == 0){
    lambda = 0.3;
    FRAC[0] = FRAC[1] = 0.24;
    printf("\n");
    printf("      -------------------------------\n");
    printf("      Network: Brightkite\n");
    printf("\n");
    printf("        - rho0_1 = %f\n", FRAC[0]);
    printf("        - rho0_2 = %f\n", FRAC[1]);
    printf("      -------------------------------\n");
    printf("\n");
  }else{
    lambda = 0.3;
    for(i=0;i<cl;i++){
      printf("\n");
      printf("      %d) Enter the initial fraction of active users in net %d: ",2+i,i+1);
      scanf("%f", &FRAC[i]);
      printf("\n");
    }
  }

  for(i=0;i<cl;i++){initial_population(&LAYERS_0[i], FRAC[i]);}


  // printf("\r");
  // printf("\n");
  // printf("\n");
  // printf("      %d) Enter the time you want to compute: ",5+(cl-1));
  // scanf("%f", &t_final);
  // printf("\n");


  printf("      -------------------------------\n");
  printf("      Parameters:\n");
  printf("\n");
  printf("        - lambda = %f\n", lambda);
  printf("        - v = %f\n", lambda/mu);
  printf("        - coupling_0  = %f\n", coupling0);
  printf("\n");
  printf("        - epsilon_0  = %f\n", epsilon);
  printf("      -------------------------------\n");
  printf("\n");
  printf("...>> spreading infection                             \n");


  clock_t begin = clock(); srand(time(NULL));

  /***********************************************************************
    BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP
    BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP
    BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP
    BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP  BIG LOOP
  ***********************************************************************/

  sprintf(pathFile, "%s/MAXepsCoex.txt", directory);
  SigmaFile = fopen(pathFile,"w");
  epsilon = epsilon0; coupling = coupling0;

  for(rep=0;rep<REPS;rep++){ ///////////////////////////////////////////////////

    sprintf(SigmaDir, "./%s/sigma_%f",EvolutionDir,coupling); if(stat(SigmaDir, &dir) == -1){mkdir(SigmaDir, 0755);}

    for(r=1;r<=RUNS;r++){ //////////////////////////////////////////////////////

      epsilon = epsilon0;
      sprintf(pathFile, "./%s/%d.txt",SigmaDir,r);
      Evolution_File = fopen(pathFile,"w");
      WriteHeader(Evolution_File,lambda,epsilon,coupling,FRAC[0],FRAC[1]);

      ResetVectors(LAYERS,LAYERS_0,cl);
      Compute_Density(LAYERS, cl);
      Update_Weight(LAYERS, cl, epsilon0, coupling);
      Update_ActRate(LAYERS, cl, lambda, mu);

      step = 0;
      delta_step = 2;

      t = 0; tstep = 0; lenss = 0;

      keep_running = true;
      s_state = false;

      // GILLESPIE METHOD //////////////////////////////////////////////////////
      //    "keep_running" condition is added to be able to
      //    start Gillespie without any active node
      //    and give enough time for the infection to grow

      while((keep_running)){

        rand1 = ((double)rand()+1)/ (double)((unsigned)RAND_MAX);
        rand2 = ((double)rand()+1)/ (double)((unsigned)RAND_MAX);

        //a_react  | a_deact  | a_mmact  | a_act
        Compute_GllspCoeffs(LAYERS, cl, sigma, G_COEFFS);

        SumTotalRate = 0.0;
        for(i=0;i<cl;i++){
          SumTotalRate += LAYERS[i].SIS.a_i.TOTAL;
        }
        tau = (1/(float)(SumTotalRate))*log(1/rand1);

        // for(i=0;i<cl;i++){
        //   printf("[%d->]  rho_A: %f  ",i+1,LAYERS[i].rho.act);
        // } printf(" - Time: %f -. rep %d/%d. run %d/%d \r",t,rep+1, REPS, r, RUNS);

        Event_Happens(LAYERS, rand2, G_COEFFS);
        Compute_Density(LAYERS, cl);
        Update_Weight(LAYERS, cl, epsilon, coupling);
        Update_ActRate(LAYERS, cl, lambda, mu);

        t += tau;

        if(t>=step){

          //only writing the case for one net each run to save time
          for(i=0;i<cl;i++){
            // i) Printing evolution of one net
            // Computations to check Stationary State
            LAYERS[i].SIS.Tss[lenss] = LAYERS[i].SIS.Nact;
          }

          if(tstep < T){
            for(i=0;i<cl;i++){
              LAYERS[i].EVOL.act[tstep][r-1] = LAYERS[i].rho.act;
              LAYERS[i].EVOL.real_time[tstep] = t;
            }
            tstep += 1;
          }

          step += delta_step; lenss += 1;
          // After saving lenss values, check if they doesn't change --> Stationary State
          if(lenss == bkt){
            s_state = StationaryState(LAYERS,cl,bkt);
            lenss = 0;
          }

          if((s_state /*&& Coex(LAYERS)*/) && ((epsilon+delta_eps)<MAXeps)){
            epsilon +=delta_eps;
            s_state = false;
          }

          fprintf(Evolution_File,"%f %f %f %f\n",
                  t,
                  LAYERS[0].rho.act,
                  LAYERS[1].rho.act,
                  epsilon);
        }

        if(t>1){
          // keep_running = (t<t_final) && !(EpidemicDied(LAYERS,cl));
          keep_running = !(CoexDied(LAYERS,cl)) && (t < 3000);
         }

      }// END OF GILLESPIE /////////////////////////////////////////////////////

    MaxEpsVector[r-1] = epsilon;

    fclose(Evolution_File);

    } // END OF RUNS ///////////////////////////////////////////////////////////

    MaxEpsAvg = ComputeMean(MaxEpsVector,RUNS);

    fprintf(SigmaFile,"%f %f\n",coupling,MaxEpsAvg);

    if(coupling + delta_coupling <= MAXcoupling){
      coupling += delta_coupling;
    }

  } // END OF REPS /////////////////////////////////////////////////////////////

  fclose(SigmaFile);
  free(G_COEFFS);

  if(EpidemicDied(LAYERS,cl)){ printf("Epidemic died.                                           \n");}

  printf("End of SIS run.                                                   \n");
  printf("\n");

  printf("...freeing\n");
  free(DEGREE_0);

  for(i=0;i<cl;i++){

    printf("...freeing 1 ");free(LAYERS_0[i].SIS.pointerEA);
    printf(" 2");free(LAYERS_0[i].SIS.pointEA_P);
    printf(" 3");free(LAYERS_0[i].SIS.pointEA_S);
    printf(" 4");free(LAYERS_0[i].SIS.NI);
    printf(" 5");free(LAYERS_0[i].SIS.pointNS);
    printf(" 6");free(LAYERS_0[i].SIS.NS);free(LAYERS_0[i].SIS.EA2pnt);
    printf(" 7");

    for(j=0;j<2;j++){
      free(LAYERS_0[i].SIS.EA[j]);
    }free(LAYERS_0[i].SIS.EA);

    printf(" 8");free(LAYERS[i].SIS.pointerEA);
    printf(" 9");free(LAYERS[i].SIS.pointEA_P);
    printf(" 10");//free(LAYERS[i].SIS.pointEA_S);
    printf(" 11");free(LAYERS[i].SIS.pointEA_SP);
    printf(" 12");free(LAYERS[i].POINT.INI);
    printf(" 13");free(LAYERS[i].POINT.FIN);
    printf(" 14");free(LAYERS[i].LINKS);
    printf(" 15");free(LAYERS[i].DEGREE);
    printf(" 16");free(LAYERS[i].NODES);
    printf(" 17");free(LAYERS[i].SIS.NI);
    printf(" 18");//free(LAYERS[i].SIS.pointNS);
    printf(" 19\n");//free(LAYERS[i].SIS.NS);
    printf("...freeing layer %d finished!\n", i+1);
    for(j=0;j<2;j++){
      free(LAYERS[i].SIS.EA[j]);
    }free(LAYERS[i].SIS.EA);
  }

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
