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
  int i,net, Nets = 6;
  char pathFile[4000];
  char *directory = malloc(5000);
  char* netfile = malloc(100);
  char *folder;
  int* DEGREE_0;
  struct stat dir = {0};
  int j,r,rep;
  int RUNS = 1; int REPS = 1; int T = 40;
  float t, tau, rand1, rand2, t_final;
  float step, delta_step;
  bool keep_running, s_state;
  float epsilon, epsilon0, MAXeps, delta_eps = 0.02;
  epsilon = epsilon0 = 0.5; MAXeps = 1.0;
  float coupling = 0.75;
  float SumTotalRate;
  float lambda, MAXlambda, sigma = 1, mu = 0.08; // rate of spontaneous activation
  struct network* LAYERS;
  struct network* LAYERS_0;
  float* G_COEFFS;
  float** PREV;
  int lenss, bkt = 10, tstep;

  int cl = 2; //Competing networks

  double time_spent; FILE *TimeFile;

  float FRAC[cl];

  FILE *Evolution_File, *StabilityFile;
  char *EvolutionDir, *StabilityDir;

  EvolutionDir = malloc(5000);
  StabilityDir = malloc(5000);

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
  printf("                O(t)   for   N E T W O R K S\n");
  printf("\n");
  printf("                                                 June 2017\n");
  printf("                                      Author: Jaume Palmer\n");
  printf("    ******************************************************\n");
  printf("\n");
  printf("\n");

  srand(time(NULL));

  TimeFile = fopen("./Ctime.txt","w");

for(net=1;net<=Nets;net++){

    clock_t begin = clock();

    sprintf(netfile, "%d", net);
    strcpy(file1, netfile); strcpy(folder, "../../nets/"); strcpy(directory,file1);
    if(stat(directory, &dir) == -1){mkdir(directory, 0755);}

    file = strcat(folder, file1);

    net1.EN = get_E_N(file);

    net1.NODES = malloc(sizeof(int)*net1.EN.N);
    net1.POINT.INI = malloc(net1.EN.N*sizeof(int));
    net1.POINT.FIN = malloc(net1.EN.N*sizeof(int));
    net1.DEGREE = malloc(net1.EN.N*sizeof(int));
    DEGREE_0 = malloc(net1.EN.maxNODE*sizeof(int));

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


    sprintf(EvolutionDir, "./%s/Evolution",directory); if(stat(EvolutionDir, &dir) == -1){mkdir(EvolutionDir, 0755);}
    sprintf(StabilityDir, "./%s/Stability",directory); if(stat(StabilityDir, &dir) == -1){mkdir(StabilityDir, 0755);}

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
      LAYERS_0[i].SIS.ppS = malloc(sizeof(int)*net1.lenLINKS);
      LAYERS_0[i].SIS.ppP = malloc(sizeof(int)*net1.lenLINKS);
      LAYERS_0[i].SIS.pointEA_SP = malloc(sizeof(int)*net1.lenLINKS);
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
      LAYERS[i].SIS.ppS = malloc(sizeof(int)*net1.lenLINKS);
      LAYERS[i].SIS.ppP = malloc(sizeof(int)*net1.lenLINKS);
      LAYERS[i].SIS.pointEA_SP = malloc(sizeof(int)*net1.lenLINKS);
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
    PREV = malloc(cl*sizeof(float*));
    for(i=0;i<cl;i++){
      PREV[i] = malloc(RUNS*sizeof(float));
    }

    lambda = 1; MAXlambda = lambda;
    // mu = 0.1;
    // sigma = 1;

    // At the beginning there is only one OSN, so initial population is built
    //  with only this
    for(i=0;i<cl;i++){
      FRAC[i] = 0.3; //Initial fraction of active users
    }

    for(i=0;i<cl;i++){initial_population(&LAYERS_0[i], FRAC[i]);}

    t_final = 30;

    printf("...>> spreading infection                             \n");

    sprintf(pathFile, "%s/Epsilons.txt", StabilityDir);
    StabilityFile = fopen(pathFile,"w");
    epsilon = epsilon0;
    WriteHeader(StabilityFile,lambda,epsilon,coupling,FRAC[0],FRAC[1]);

    for(rep=0;rep<REPS;rep++){ ///////////////////////////////////////////////////

      sprintf(pathFile, "./%s/%d.txt",EvolutionDir,rep);
      Evolution_File = fopen(pathFile,"w");
      WriteHeader(Evolution_File,lambda,epsilon,coupling,FRAC[0],FRAC[1]);

      for(r=1;r<=RUNS;r++){ //////////////////////////////////////////////////////

        ResetVectors(LAYERS,LAYERS_0,cl);
        Compute_Density(LAYERS, cl);
        Update_Weight(LAYERS, cl, epsilon0, coupling);
        Update_ActRate(LAYERS, cl, lambda, mu);

        step = 0; delta_step = t_final/(float)50;
        t = 0; tstep = 0; lenss = 0;
        keep_running = true; s_state = false;

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

          Event_Happens(LAYERS, rand2, G_COEFFS);
          Compute_Density(LAYERS, cl);
          Update_Weight(LAYERS, cl, epsilon, coupling);
          Update_ActRate(LAYERS, cl, lambda,mu);

          t += tau;

          if(t>=step){
            //only writing the case for one net each run to save time
            for(i=0;i<cl;i++){
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

            if((s_state && Coex(LAYERS)) && ((epsilon+delta_eps)<MAXeps)){
              epsilon +=delta_eps;
            }
          }

          if(t>1){
            keep_running = (t<t_final) && !(EpidemicDied(LAYERS,cl));
           }
        }// END OF GILLESPIE /////////////////////////////////////////////////////


      } // END OF RUNS ///////////////////////////////////////////////////////////

      for(i=0;i<cl;i++){
        for(j=0;j<tstep-1;j++){
          LAYERS[i].EVOL.act_avg[j] = ComputeMean(LAYERS[i].EVOL.act[j],RUNS);
        }
      }

      if((lambda+MAXlambda/(float)REPS)<MAXlambda){
        lambda += MAXlambda/(float)REPS;
      }
      // Printing Evolution_File: it contains the averaged evolution of both nets
      if(cl == 2){
        for(j=0;j<tstep-1;j++){
            fprintf(Evolution_File,"%f %f %f\n",
                    LAYERS[0].EVOL.real_time[j],
                    LAYERS[0].EVOL.act_avg[j],
                    LAYERS[1].EVOL.act_avg[j]);
        }

        fprintf(StabilityFile,"%f %f %f\n",
                LAYERS[0].rho.act,
                LAYERS[1].rho.act,
                epsilon);
      }


    } // END OF REPS /////////////////////////////////////////////////////////////

    fclose(Evolution_File);
    fclose(StabilityFile);
    free(G_COEFFS);

    if(EpidemicDied(LAYERS,cl)){ printf("Epidemic died.                                           \n");}

    free(DEGREE_0);

    for(i=0;i<cl;i++){

      free(LAYERS_0[i].SIS.pointerEA); free(LAYERS_0[i].SIS.pointEA_P);
      free(LAYERS_0[i].SIS.pointEA_S); free(LAYERS_0[i].SIS.NI); free(LAYERS_0[i].SIS.pointEA_SP);
      free(LAYERS_0[i].SIS.pointNS); free(LAYERS_0[i].SIS.NS); free(LAYERS_0[i].SIS.EA2pnt);

      for(j=0;j<2;j++){free(LAYERS_0[i].SIS.EA[j]);}free(LAYERS_0[i].SIS.EA);

      free(LAYERS[i].SIS.pointerEA);free(LAYERS[i].SIS.pointEA_P);
      free(LAYERS[i].SIS.pointEA_S);
      free(LAYERS[i].SIS.ppS);free(LAYERS[i].SIS.ppP);
      free(LAYERS[i].SIS.pointEA_SP);free(LAYERS[i].POINT.INI);
      free(LAYERS[i].POINT.FIN);free(LAYERS[i].LINKS);
      // free(LAYERS[i].pLINKS);
      free(LAYERS[i].DEGREE);free(LAYERS[i].NODES);
      free(LAYERS[i].SIS.NI);//free(LAYERS[i].SIS.pointNS);
      free(LAYERS[i].SIS.NS);
      for(j=0;j<2;j++){free(LAYERS[i].SIS.EA[j]);}free(LAYERS[i].SIS.EA);
    }

    free(net1.NODES);free(net1.POINT.INI);
    free(net1.POINT.FIN);free(net1.DEGREE);
    free(net1.LINKS);//free(net1.pLINKS);

    free(LAYERS);free(LAYERS_0);
    // Memory leaks !!
    //free(folder);
    //free(file);
    //free(file1);
    //free(directory1);
    //free(directory2);
    //free(directory3);
    //free(directory4);


    clock_t end = clock();
    time_spent = (double)(end-begin)/CLOCKS_PER_SEC;

    fprintf(TimeFile, "%d %f\n",net1.EN.N,time_spent);

  }
  fclose(TimeFile);
  printf("DONE !! \n");

  return 0;
}
