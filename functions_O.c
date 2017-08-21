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
#include "functions.h"

/*
struct result {
  int E;
  int N;
  int* NODES;
  int maxNODE;
};

struct state {
  int s2a;
  int p2a;
};

struct network {
  struct result EN;
  struct activation SIS;
  struct pointers POINT;
  int* NODES;
  int* LINKS;
  int* DEGREE;
  int* WEIGHT;
};

struct pointers {
  int* INI;
  int* FIN;
};

struct activation {
  int Nact;
  int Eact;
  int* NI;
  int** EA;
  int* pointerEA;
};

*/


bool inVECTOR(int k, int size, int* VECTOR){

  int j;
  bool found;

  j=0;
  found = false;
  while((j< size)&&(!found)){
    found = (k == VECTOR[j]);
    j += 1;
  }

  return found;

}

struct result get_E_N(char *filename){

  char str[100]; // 100 is an arbitrary number representing maximum line length to read
  //IMPORTANT: to be large enough to fit two large numbers (~10⁶) in one line
  char *line, *token; // line and token are pointers to strings of characters
  int i, node_i,node_j, count;
  struct result EN; // result = {E,N} with E = # of links, N = # of nodes
  bool found_i, found_j;
  float progress;

  FILE* fp;

  filename=strcat(filename,".txt");
  fp = fopen(filename,"r+");
  if (fp == NULL){
    printf("ERROR!!! File path not found!!!\n");
    printf("\n");
    exit(0);
  }

  EN.E = 0;

  line=fgets(str,100,fp);

  printf("...>> getting E\n");
  while (line){ // Canviar line != NULL per line m'ha solucionat un error
    EN.E += 1;
    line = fgets(str,100,fp);
  }

  fclose(fp);

  fp = fopen(filename,"r");

  line=fgets(str,100,fp);
  EN.NODES = malloc((EN.E+1)*sizeof(int)); //max possible number of nodes is E+1 if there is E links
  for (i=0;i<(EN.E+1);i++){
    EN.NODES [i]= 0;
  }
  i = 0;
  EN.N = 0; EN.maxNODE = 0;
  printf("...>> getting N\n");
  count = 0;
  while (line != NULL){
    token = strtok(line," "); // Reading the first token
    node_i = atoi(token); // Parsing the char pointer to int
    found_i = inVECTOR(node_i,(EN.E+1),EN.NODES);
    if (!found_i){
      EN.NODES[i] = node_i;
      if (node_i >= EN.maxNODE){
        EN.maxNODE = node_i;
      }
      EN.N += 1;
      i += 1;
    }

    token = strtok(NULL," "); //reading following tokens
    node_j = atoi(token);  // Parsing the char pointer to int
    found_j = inVECTOR(node_j,(EN.E+1),EN.NODES);
    if (!found_j){
      EN.NODES[i] = node_j;
      if (node_j >= EN.maxNODE){
        EN.maxNODE = node_j;
      }
      EN.N += 1;
      i += 1;
    }
    count += 1;
    progress = count*100/(float)EN.E;
    printf("Progress %f %c\r", progress, '%');
    line = fgets(str,100,fp); //Read next line
  }


  rewind(fp); //free(EN.NODES);
  fclose(fp);
  printf("\n");

  return EN;
}

void get_pointers(struct network *net, char* filename, int* DEGREE_0){

  // IMPORTANT: this function don't need the network txt to be ordered

  char str[100];
  char* line;
  char* token;
  int node_i, node_j, i, index;

  FILE* fp;
  fp = fopen(filename,"r");

  // Preparing degree for later filling
  for (i=0;i<net->EN.maxNODE;i++){
    DEGREE_0[i] = 0;
  }

  line=fgets(str,100,fp); // Setting first
  net->loops = 0;
  while (line){
    token = strtok(line," "); // Reading the first token
    node_i = atoi(token); // Parsing the char pointer to int
    token = strtok(NULL," "); //reading following tokens
    node_j = atoi(token);  // Parsing the char pointer to int

    if (node_i != node_j){ // This if works to avoid counting loops two times
      DEGREE_0[node_i-1] +=1;
      DEGREE_0[node_j-1] +=1;
    }
    else if(node_i == node_j){
      net->loops += 1;
    }
    line = fgets(str,100,fp); // Read next line
  }
  fclose(fp);


  net->POINT.INI[0] = 0;
  index=1;
  for (i=1;i<net->EN.maxNODE;i++){
    if (DEGREE_0[i-1] != 0){
      net->POINT.INI[index] = net->POINT.INI[index-1] + DEGREE_0[i-1]; //If there are missing nodes, we don't want to sum them
      index +=1;
    }
  }

  for (i=0;i<net->EN.N;i++){
    net->POINT.FIN[i] = net->POINT.INI[i]-1;
  }

}

int cmpfunc_c(const void * a, const void * b){
   return ( *(int*)a - *(int*)b );
 } // Function to be used for qsort()

int get_index(int ind, int* VECTOR, int len){

  // NEED VECTOR SORTED!!
  int k, middle, first, last;
  // int i;

  first = 0;
  last = len - 1;
  middle = (first+last)/2;

  while (first <= last) {
    if (VECTOR[middle] < ind){
      first = middle + 1;
    }else if (VECTOR[middle] == ind){
      k = middle;
      break;
    }else{
      last = middle - 1;
    }
    middle = (first + last)/2;
  }

  if (first > last){
    printf("ERROR!!!! Not found! %d is not present in the list.\n", ind);
    k = INT_MAX;

    // for(i=0;i<len;i++){
    //   printf("%d ",VECTOR[i]);
    // }printf("\n");
    exit(1);
  }

  return k;
}

void last_link(struct network *net, char* filename){

  char str[100];
  char* line;
  char* token;
  int i, j, node_j, node_i, index_i, index_j;
  int count;
  float progress;
  bool found_i, found_j;

  FILE* fp;
  fp = fopen(filename,"r");

  // Counters to keep track of the index for filling LINKS

  for (i=0;i<net->EN.N;i++){
    net->DEGREE[i] = 0;
  }

  net->repetitions = 0;
  count = 0;
  line=fgets(str,100,fp); // Setting first
  while (line != NULL){
    token = strtok(line," "); // Reading the first token
    node_i = atoi(token); // Parsing the char pointer to int
    token = strtok(NULL," "); //reading following tokens
    node_j = atoi(token);  // Parsing the char pointer to int
    if(node_i == 0){printf("here: node_j = %d and line is around %d \n", node_j, count);exit(1);}
    if(node_j == 0){printf("no, here\n");exit(1);}
    index_i = get_index(node_i, net->NODES, net->EN.N);
    if (node_i != node_j){
      index_j = get_index(node_j, net->NODES, net->EN.N);
      // INDEX i
      if (net->POINT.FIN[index_i]>=net->POINT.INI[index_i]){ // there is already something written
        i = net->POINT.INI[index_i];
        found_i = false;
        while((!found_i) && (i<=net->POINT.FIN[index_i])){
          found_i = (node_j == net->LINKS[i]); //look if the link is already in place
          i += 1;
        }
        if(!found_i){
          net->POINT.FIN[index_i] += 1;
          net->LINKS[net->POINT.INI[index_i]+net->DEGREE[index_i]] = node_j;
          net->DEGREE[index_i] +=1;
        }
        else{
          net->repetitions +=1;
        }
      }
      else{
        net->POINT.FIN[index_i] += 1;
        net->LINKS[net->POINT.INI[index_i]+net->DEGREE[index_i]] = node_j;
        net->DEGREE[index_i] +=1;
      }
      // INDEX j
      if (net->POINT.FIN[index_j]>=net->POINT.INI[index_j]){
        j = net->POINT.INI[index_j];
        found_j = false;
        while((!found_j) && (j<=net->POINT.FIN[index_j])){
          found_j = (node_i == net->LINKS[j]); //look if the link is already in place
          j += 1;
        }
        if(!found_j){
          net->POINT.FIN[index_j] += 1;
          net->LINKS[net->POINT.INI[index_j]+net->DEGREE[index_j]] = node_i;
          net->DEGREE[index_j] +=1;
        }
      }
      else{
        net->POINT.FIN[index_j] += 1;
        net->LINKS[net->POINT.INI[index_j]+net->DEGREE[index_j]] = node_i;
        net->DEGREE[index_j] +=1;
      }
    }
    else{
      //If node_i = node_j there is a loop, we just skip it
      //printf("ERROR!!!! Node i is the same as node j!!!!!!!\n");
    }

    count += 1;
    progress = count*100/(float)net->EN.E;
    printf("Progress %f %c\r", progress, '%');
    line = fgets(str,100,fp);
  }
  fclose(fp);
}

bool link_is_once(int* link, int** MATRIX, int len){

  int i, count;
  bool is_once;

  count = 0;
  for(i=0;i<len;i++){
    if(((MATRIX[0][i] == link[0]) && (MATRIX[1][i] == link[1]))
      || ((MATRIX[1][i] == link[0]) && (MATRIX[0][i] == link[1]))){
        count +=1;
      }
  }

  is_once = (count == 1);

  return is_once;
}

bool check_file(char* filename, struct network *net){

  int i, j, node_i, node_j, count;
  char str[100];
  char *line, *token;
  int** MATRIXfile;
  int* link;
  bool sf_sg;
  float progress;

  FILE* fr = fopen(filename, "r");

  line=fgets(str,100,fr); // Setting first

  MATRIXfile = malloc(  net->EN.E*sizeof(int*));
  for(i=0;i<net->EN.E;i++){
    MATRIXfile[i] = malloc(2*sizeof(int));
  }

  // 1) Writing a matrix with the links in the file
  printf("...storing the file in array\n");
  i=0;
  while (line != NULL){
    token = strtok(line," "); // Reading the first token
    node_i = atoi(token); // Parsing the char pointer to int
    token = strtok(NULL," "); //reading following tokens
    node_j = atoi(token);  // Parsing the char pointer to int
    MATRIXfile[0][i] = node_i;
    MATRIXfile[1][i] = node_j;
    line = fgets(str,100,fr);
    i+=1;
    progress = i*100/(float)net->EN.E;
    printf("Progress %f %c\r", progress, '%');
  }

  fclose(fr);

  // 2) Checking each link in the generated vectors are once and only once in the matrix
  printf("...checking links in stored array\n");
  count = 0;
  sf_sg = true;
  for(i=0;i<net->EN.N;i++){
    for(j=net->POINT.INI[i];j<net->POINT.FIN[i]+1;j++){
      link = malloc(2*sizeof(int));
      link[0] = net->NODES[i];
      link[1] = net->LINKS[j];
      sf_sg = link_is_once(link, MATRIXfile, net->EN.E);
      free(link);
      count += 1;
      progress = count*100/(float)(2*net->EN.E); // IMPORTANT: if loops, result won't add to 100
      printf(".\r");
      printf("..\r");
      printf("...\r");
    }
  }
  printf("\n");
  for(i=0;i<net->EN.E;i++){
    free(MATRIXfile[i]);
  }
  free(MATRIXfile);

  return sf_sg;
}

void OverwritePointer(struct network *net, int af,int lenLINKS,int tmp,int neighbor){

  int i,l, index_act, index_act1, index_act2;

  index_act = get_index(net->LINKS[neighbor],net->NODES,net->EN.N);
  for(l=net->POINT.INI[index_act];l<=net->POINT.FIN[index_act];l++){
    if(net->LINKS[l] == af){
      net->SIS.pointerEA[l] = 0;
    }
  }
  index_act1 = get_index(net->SIS.EA[0][tmp-1],net->NODES,net->EN.N);
  for(l=net->POINT.INI[index_act1];l<=net->POINT.FIN[index_act1];l++){
    if(net->LINKS[l] == net->SIS.EA[1][tmp-1]){
      net->SIS.pointerEA[l] = tmp;
    }
  }
  index_act2 = get_index(net->SIS.EA[1][tmp-1],net->NODES,net->EN.N);
  for(l=net->POINT.INI[index_act2];l<=net->POINT.FIN[index_act2];l++){
    if(net->LINKS[l] == net->SIS.EA[0][tmp-1]){
      net->SIS.pointerEA[l] = tmp;
    }
  }
  for(i=0;i<lenLINKS;i++){
    if(net->SIS.pointerEA[i] == net->SIS.Eact){
      net->SIS.pointerEA[i] = 0;
    }
  }
}

void Event_Happens(struct network *LAYERS, float inf, float* COEFFS){ //ONLY WORKING FOR 2 NETS

  /********************************************************************************
  function intended to replace modify NI and EA depending on the event:

    - inf = 1 if activation
    - inf = 0 if recovery

  *********************************************************************************/

  int i,j;

  for(i = 0; i<LAYERS[0].SIS.Eact; i++){
    if(LAYERS[0].SIS.EA[0][i] == 0){
      printf("WARNING!!!!!!\n");
      printf("%f\n", inf);
      for(j=0;j<8;j++){
        printf("%f ", COEFFS[j]);
      }
      printf("\n");
      for(j=0;j<8;j++){
        printf("%f ", COEFFS[j]);
      }
      printf("\n");printf("\n");
      break;
    }
  }

    if((inf <= COEFFS[0]) && ((LAYERS[0].SIS.count.p2a > 0) && (LAYERS[0].SIS.Nact < LAYERS[0].EN.N))){ // viral reactivation

      Reactivation(&LAYERS[0]);

    }
    else if (((COEFFS[0] < inf) && (inf <= COEFFS[1])) && (LAYERS[0].SIS.Nact > 0)){ // deactivation

      Deactivation(&LAYERS[0]);

    }
    else if (((COEFFS[1] < inf)&&(inf <= COEFFS[2])) && ((LAYERS[0].SIS.Nsus > 0)&&(LAYERS[0].SIS.Nact<LAYERS[0].EN.N))){ // m.m. activation

      MMactivation(&LAYERS[0]);

    }
    else if (((COEFFS[2]<inf)&&(inf <= COEFFS[3])) && ((LAYERS[0].SIS.count.s2a > 0) && (LAYERS[0].SIS.Nact < LAYERS[0].EN.N))){ // viral activation

      Activation(&LAYERS[0]);

    }
    else if(((COEFFS[3]<inf) && (inf <= COEFFS[4])) &&((LAYERS[1].SIS.count.p2a > 0) &&(LAYERS[1].SIS.Nact < LAYERS[1].EN.N))){ // viral reactivation

      Reactivation(&LAYERS[1]);

    }
    else if (((COEFFS[4] < inf) && (inf <= COEFFS[5])) && (LAYERS[1].SIS.Nact > 0)){ // deactivation

      Deactivation(&LAYERS[1]);

    }
    else if (((COEFFS[5] < inf)&&(inf <= COEFFS[6])) && ((LAYERS[1].SIS.Nsus > 0)&&(LAYERS[1].SIS.Nact < LAYERS[1].EN.N))){ // m.m. activation

      MMactivation(&LAYERS[1]);

    }
    else if (((COEFFS[6]<inf)&&(inf <= COEFFS[7])) && ((LAYERS[1].SIS.count.s2a > 0) && (LAYERS[1].SIS.Nact < LAYERS[1].EN.N))){ // viral activation

      Activation(&LAYERS[1]);

    }


}

void initial_population(struct network *net, float fraction){

    ///////////////////////////////////////////////////////////////
    // Initial infected population

  int i, j, k, index_inf, index_sus, tmp, count;
  int infected, act;
  bool found2;

  net->SIS.Nact = 0;
  net->SIS.Nsus = net->EN.N;
  net->SIS.Eact = 0; // Active Edges
  net->SIS.count.s2a = 0;
  net->SIS.count.p2a = 0;

  for(i=0;i<net->EN.N;i++){
    net->SIS.NI[i] = 0; // I do this bc I don't know what is there in the allocated memory
  }

  for(i=0;i<net->lenLINKS;i++){
    net->SIS.pointerEA[i] = net->SIS.pointEA_SP[i] = 0; // I do this bc I don't know what is there in the allocated memory
    net->SIS.ppS[i] = net->SIS.ppP[i] = 0;
  }
  for(i=0;i<net->EN.E;i++){
    net->SIS.EA[0][i] = 0;
    net->SIS.EA[1][i] = 0;
  }

  for(i=0;i<net->EN.N;i++){
    net->SIS.NS[i] = 0;
    net->SIS.pointNS[i] = i+1;
  }


  for(k=0;k<fraction*net->EN.N;k++){

    printf("Progress: %f %c  \r",k*100/(float)(fraction*net->EN.N),'%');
  //printf("Progress: %f  %c           \r",k*100/(float)net->SIS.Nact,'%');
    index_sus = rand()%net->SIS.Nsus;
    index_inf = net->SIS.pointNS[index_sus]-1;
    infected = net->NODES[index_inf];

    // printf("Infected: %d \n",infected);
    net->SIS.Nact += 1;
    net->SIS.Nsus -= 1;

    net->SIS.pointNS[index_sus] = net->SIS.pointNS[net->SIS.Nsus];

    net->SIS.NI[k] = infected;

    for(j=net->POINT.INI[index_inf];j<=net->POINT.FIN[index_inf];j++){

      speedUpdateActivation(net,j,infected);

    }

    net->SIS.NS[index_inf] = 1; // by definition of NS only nodes still non active can be chosen

  }

}

void WriteNet_Act(struct network net, char* directory, int filecount){

  int i,j,node, index;
  char pathFile[1000];


  snprintf(pathFile,1000,"%s/net%d.txt",directory, filecount);
  FILE* f = fopen(pathFile,"w");

  int* WRITTEN = malloc(net.SIS.Nact*sizeof(int));
  for(i=0;i<net.SIS.Nact;i++){
    WRITTEN[i] = 0;
  }

  for(i=0;i<net.SIS.Nact;i++){
    node = net.SIS.NI[i];
    index = get_index(node,net.NODES,net.EN.N);
    for(j=net.POINT.INI[index];j<=net.POINT.FIN[index];j++){
      if(((!inVECTOR(net.LINKS[j],net.SIS.Nact,WRITTEN)) && (inVECTOR(net.LINKS[j],net.SIS.Nact,net.SIS.NI)))
      && ((!inVECTOR(node,net.SIS.Nact,WRITTEN)) && (inVECTOR(node,net.SIS.Nact,net.SIS.NI)))){
        //if LINKS[j] is in vector it means the link has been written previously
        // We only print the net of active nodes connected
        fprintf(f,"%d %d\n", node, net.LINKS[j]);
      }
    }
    WRITTEN[i] = node;
  }

  fclose(f);
  free(WRITTEN);
}

void WriteNet_Total(struct network net, char* directory, int filecount){

  int i,j,node, index, index_ngbr, total_nodes;
  char pathFile[1000];


  snprintf(pathFile,1000,"%s/net%d.txt",directory, filecount);
  FILE* f = fopen(pathFile,"w");

  total_nodes = net.EN.N-net.SIS.Nsus;

  int* WRITTEN = malloc(total_nodes*sizeof(int));
  for(i=0;i<total_nodes;i++){
    WRITTEN[i] = 0;
  }

  for(i=0;i<net.SIS.Nact;i++){
    node = net.SIS.NI[i];
    index = get_index(node,net.NODES,net.EN.N);
    for(j=net.POINT.INI[index];j<=net.POINT.FIN[index];j++){
      index_ngbr = get_index(net.LINKS[j],net.NODES,net.EN.N);
      if(((!inVECTOR(net.NODES[index_ngbr],total_nodes,WRITTEN)) && (net.SIS.pointNS[index_ngbr] != 0))
        && (((!inVECTOR(node,total_nodes,WRITTEN)) && (net.SIS.pointNS[index] != 0)))){
        //if LINKS[j] is in vector it means the link has been written previously
        // We only print the net of active nodes connected
        fprintf(f,"%d %d\n", node, net.NODES[index_ngbr]);
      }
    }
    WRITTEN[i] = node;
  }

  fclose(f);
  free(WRITTEN);
}

struct state ComputeEs2A(struct network *net){

  int i;
  int index, neighbor;
  int Es2a, Ep2a;
  struct state online;

  Es2a = 0;
  Ep2a = 0;
  for (i=0;i<net->SIS.Eact;i++){
    neighbor = net->SIS.EA[1][i];
    index = get_index(neighbor,net->NODES,net->EN.N);
    if(net->SIS.NS[index]==0){
      Es2a += 1;
    }
    else{
      Ep2a += 1;
    }
  }

  //using result structure only for convenience
  online.s2a = Es2a;
  online.p2a = Ep2a;

  return online;

}

void Update_pointEA(struct network *LAYERS, int cl){

  int i, index, j,ks,kp;

   for(i=0;i<cl;i++){

    ks = 0;
    kp = 0;

    for(j=0;j<LAYERS[i].SIS.Eact;j++){

      index = LAYERS[i].SIS.EA2pnt[j];
      if(LAYERS[i].SIS.pointEA_SP[index] == 1){

        LAYERS[i].SIS.pointEA_S[ks] = index;
        ks += 1;


      } else if(LAYERS[i].SIS.pointEA_SP[index] == 2){

        LAYERS[i].SIS.pointEA_P[kp] = index;
        kp += 1;
      }
    }

    LAYERS[i].SIS.count.s2a = ks;
    LAYERS[i].SIS.count.p2a = kp;

    if(LAYERS[i].SIS.count.p2a + LAYERS[i].SIS.count.s2a != LAYERS[i].SIS.Eact){

      printf("\n");
      printf("------------------------------------------------\n");
      printf("WARNING!!!! Not counting active links properly\n");
      printf("\n");
      printf("NET %d \n",i+1);
      printf("\n");
      printf("  Es2a: %d \n", LAYERS[i].SIS.count.s2a);
      printf("  Ep2a: %d \n", LAYERS[i].SIS.count.p2a);
      printf(" _____________\n");
      printf("  Eact: %d \n", LAYERS[i].SIS.Eact);
      printf("------------------------------------------------\n");
      printf("\n");

      exit(1);
    }

  }

}

float ComputeMean(float* VECTOR,int len){

  int i;
  float mean;

  mean = 0;
  for(i=0;i<len;i++){
    mean += VECTOR[i];
  }

  return mean/(float)len;
}

void ResetVectors(struct network *LAYERS, struct network *LAYERS_0, int cl){

  int i,k;

  for (i=0;i<cl;i++){

    LAYERS[i].EN.N = LAYERS_0[i].EN.N;
    LAYERS[i].EN.E = LAYERS_0[i].EN.E;
    LAYERS[i].lenLINKS = LAYERS_0[i].lenLINKS;
    LAYERS[i].SIS.Nact = LAYERS_0[i].SIS.Nact;
    LAYERS[i].SIS.Eact = LAYERS_0[i].SIS.Eact;
    LAYERS[i].SIS.Nsus = LAYERS_0[i].SIS.Nsus;
    LAYERS[i].SIS.Nact = LAYERS_0[i].SIS.Nact;
    LAYERS[i].loops = LAYERS_0[i].loops;
    LAYERS[i].zeros = LAYERS_0[i].zeros;
    LAYERS[i].repetitions = LAYERS_0[i].repetitions;

    LAYERS[i].SIS.count.s2a = LAYERS_0[i].SIS.count.s2a;
    LAYERS[i].SIS.count.p2a = LAYERS_0[i].SIS.count.p2a;


    for(k=0;k<LAYERS_0[i].EN.E;k++){
      LAYERS[i].SIS.EA[1][k] = LAYERS_0[i].SIS.EA[1][k];
      LAYERS[i].SIS.EA[0][k] = LAYERS_0[i].SIS.EA[0][k];
    }
    for(k=0;k<LAYERS_0[i].EN.N;k++){
      LAYERS[i].NODES[k] = LAYERS_0[i].NODES[k];
      LAYERS[i].DEGREE[k] = LAYERS_0[i].DEGREE[k];
      LAYERS[i].POINT.INI[k] = LAYERS_0[i].POINT.INI[k];
      LAYERS[i].POINT.FIN[k] = LAYERS_0[i].POINT.FIN[k];
      LAYERS[i].SIS.NS[k] = LAYERS_0[i].SIS.NS[k];
      LAYERS[i].SIS.NI[k] = LAYERS_0[i].SIS.NI[k];
      LAYERS[i].SIS.pointNS[k] = LAYERS_0[i].SIS.pointNS[k];
    }


    for(k=0;k<LAYERS_0[i].lenLINKS;k++){
      LAYERS[i].LINKS[k] = LAYERS_0[i].LINKS[k];
      LAYERS[i].pLINKS = LAYERS_0[i].pLINKS;
      LAYERS[i].SIS.pointerEA[k] = LAYERS_0[i].SIS.pointerEA[k];
      LAYERS[i].SIS.pointEA_S[k] = LAYERS_0[i].SIS.pointEA_S[k];
      LAYERS[i].SIS.pointEA_P[k] = LAYERS_0[i].SIS.pointEA_P[k];
      LAYERS[i].SIS.pointEA_SP[k] = LAYERS_0[i].SIS.pointEA_SP[k];
      LAYERS[i].SIS.ppS[k] = LAYERS_0[i].SIS.ppS[k];
      LAYERS[i].SIS.ppP[k] = LAYERS_0[i].SIS.ppP[k];
      LAYERS[i].SIS.EA2pnt[k] = LAYERS_0[i].SIS.EA2pnt[k];
    }

  }
}

void InitializeNet(struct network *net, char* file, char* directory){

  int i, j, maxDEGREE, sumDEGREE, Kcount;
  float cummulateP;
  bool same_file;
  char pathFile[4000];
  char *answer;

  answer = malloc(10);

  FILE *fw;

  net->zeros = net->EN.maxNODE - net->EN.N; // Count the number of nodes with degree 0

  printf("\n");
  printf("      ---------------------------\n");
  printf("      Number of nodes (N) is %d\n", net->EN.N);
  printf("      Number of links (E) is %d\n", net->EN.E);
  printf("\n");
  printf("      (%d missing nodes)\n", net->zeros); //A missing node is a node with degree 0
  printf("      ---------------------------\n");
  printf("\n");


  /*
    For the last_link I need a vector, sorted from min to max index with the
    different nodes in the network, appearing only once each of them. These nodes
    are already in EN.NODES however EN.NODES is a larger vector, as it is defined
    before knowing the exact number of N. For this reason I define another vector
    of length N for accuracy.
  */

  for (i=0;i<net->EN.N;i++){
    net->NODES[i] = net->EN.NODES[i];
  }

  free(net->EN.NODES);
  /*
  newDEGREE contains the same information that DEGREE but without the zeros,
  also it has the repetitions corrected
  */

  qsort(net->NODES, net->EN.N, sizeof(int),cmpfunc_c);

  net->repetitions = 0;
  printf("...filling LINKS\n");
  last_link(net,file);
  printf("\n");
  if (net->repetitions != 0){
    for(i=0;i<net->EN.N;i++){
      net->DEGREE[i] = net->POINT.FIN[i] - net->POINT.INI[i]+1;
    }
  }

  if (net->repetitions != 0){
    printf("\n");
    printf("      -------------------------------------------\n");
    printf("      %d repetitions have been found and cleaned\n", net->repetitions);
    printf("      -------------------------------------------\n");
    printf("\n");
  }

  maxDEGREE = net->DEGREE[0]; // Looking for the max degree
  for(i=1;i<net->EN.N;i++){
    if(net->DEGREE[i] > maxDEGREE){
      maxDEGREE = net->DEGREE[i];
    }
  }

  net->DEG.P = malloc((maxDEGREE+1)*sizeof(float));
  net->DEG.Pc = malloc((maxDEGREE+1)*sizeof(float));
  net->DEG.P[0] = 0;

  for(i=1;i<=maxDEGREE;i++){
    Kcount = 0;
    for(j=0;j<net->EN.N;j++){
      if(net->DEGREE[j] == i){
        Kcount += 1; // Counting all nodes, that have a degree = i
      }
    }
    net->DEG.P[i] = Kcount/(float)net->EN.N;
  }

  sprintf(pathFile,"%s/degree.txt",directory);
  fw = fopen(pathFile, "w");
  if(fw == NULL){
    printf("Error opening file!\n");
    exit(1);
  }
  for(i=0;i<=maxDEGREE;i++){ // i is the degree
    cummulateP = 0.0;
    for(j=i;j<=maxDEGREE;j++){
      cummulateP += net->DEG.P[j];
    }
    if (net->DEG.P[i] != 0){ //Excluding the degrees with probability 0
      net->DEG.Pc[i] = cummulateP;
      fprintf(fw,"%d %f %f\n", i, net->DEG.P[i], net->DEG.Pc[i]);
    }
  }

  fclose(fw);


  sumDEGREE = 0;
  for (i=0;i<net->EN.N;i++){
    sumDEGREE += net->DEGREE[i];
  }

  if (sumDEGREE != (net->lenLINKS-2*net->repetitions)){
    printf("\n");
    printf("ERROR!!! DEGREE does not match size of LINKS !!\n");
    printf("sumDEGREE = %d || lenLINKS - 2*repetitions = %d \n", sumDEGREE, (net->lenLINKS-2*net->repetitions));

  }

  strcpy(answer, "N");

  if(strcmp(answer, "Y") == 0){
    printf("\n");
    printf("sit back and relax, it can take a while...\n");
    printf(" \n");
    printf("...checking vector output matches initial file\n");
    same_file = check_file(file,net);
    if (same_file){
      printf("\n");
      printf("      CONGRATULATIONS!!! File read and created are the same!!\n");
    }
    else{
      printf(" \n");
      printf("      Too bad :( ... Better luck next time\n");
    }
  }
  else{
    printf(" \n");
  }


  free(answer);

}

void Compute_GllspCoeffs(struct network *LAYERS,int cl, float lambda, float sigma, float mu, float* COEFFS){

  int i;
  int A_0;

  for(i=0;i<4*cl;i++){
    COEFFS[i] = 0;
  }

  A_0 = 0;
  for(i=0;i<cl;i++){
    LAYERS[i].SIS.a_i.TOTAL = LAYERS[i].lambda*LAYERS[i].SIS.count.s2a
                              + LAYERS[i].lambda*LAYERS[i].SIS.count.p2a
                              + sigma*LAYERS[i].SIS.Nact
                              + mu*LAYERS[i].SIS.Nsus;

    A_0 += LAYERS[i].SIS.a_i.TOTAL;

  }


  for (i=0;i<cl;i++){

    LAYERS[i].SIS.a_i.react = LAYERS[i].lambda*LAYERS[i].SIS.count.p2a/(float)A_0;
    LAYERS[i].SIS.a_i.deact = sigma*LAYERS[i].SIS.Nact/(float)A_0;
    LAYERS[i].SIS.a_i.mmact = mu*LAYERS[i].SIS.Nsus/(float)A_0;
    LAYERS[i].SIS.a_i.act = LAYERS[i].lambda*LAYERS[i].SIS.count.s2a/(float)A_0;

    COEFFS[4*(cl-1)*i] = LAYERS[i].SIS.a_i.react;
    COEFFS[4*(cl-1)*i+1] = LAYERS[i].SIS.a_i.deact;
    COEFFS[4*(cl-1)*i+2] = LAYERS[i].SIS.a_i.mmact;
    COEFFS[4*(cl-1)*i+3] = LAYERS[i].SIS.a_i.act;

  }

  for(i=1;i<4*cl;i++){
    COEFFS[i] = COEFFS[i] + COEFFS[i-1];
  }


}

void Compute_Density(struct network *LAYERS, int cl){

  int i;

  for(i=0;i<cl;i++){
    LAYERS[i].rho.act = LAYERS[i].SIS.Nact/(float)LAYERS[i].EN.N;
    LAYERS[i].rho.susc = LAYERS[i].SIS.Nsus/(float)LAYERS[i].EN.N;
    LAYERS[i].rho.pas = 1 - LAYERS[i].rho.act - LAYERS[i].rho.susc;



    if((LAYERS[i].rho.act >1) || (LAYERS[i].rho.act < -1/(float)LAYERS[i].EN.N)){
      printf("----------------------------------\n");
      printf("  ERROR in layer %d !!!!!\n", i+1);
      printf("    rho_act =  %f\n", LAYERS[i].rho.act);
      printf("    rho_pas =  %f\n", LAYERS[i].rho.pas);
      printf("    rho_susc =  %f\n", LAYERS[i].rho.susc);
      printf("----------------------------------\n");

      exit(1);
    }
    if((LAYERS[i].rho.pas > 1) || (LAYERS[i].rho.pas < -1/(float)LAYERS[i].EN.N)){
      printf("----------------------------------\n");
      printf("  ERROR in layer %d !!!!!\n", i+1);
      printf("    rho_act =  %f\n", LAYERS[i].rho.act);
      printf("    rho_pas =  %f\n", LAYERS[i].rho.pas);
      printf("    rho_susc =  %f\n", LAYERS[i].rho.susc);
      printf("----------------------------------\n");

      exit(1);
    }
    if((LAYERS[i].rho.susc > 1) || (LAYERS[i].rho.susc < -1/(float)LAYERS[i].EN.N)){
      printf("----------------------------------\n");
      printf("  ERROR in layer %d !!!!!\n", i+1);
      printf("    rho_act =  %f\n", LAYERS[i].rho.act);
      printf("    rho_pas =  %f\n", LAYERS[i].rho.pas);
      printf("    rho_susc =  %f\n", LAYERS[i].rho.susc);
      printf("----------------------------------\n");

      exit(1);
    }

  }

}

void Update_ActRate(struct network *LAYERS, int cl, float lambda){

  int i;

  for(i=0;i<cl;i++){
    LAYERS[i].lambda = lambda*LAYERS[i].weight;

  }
}

void Update_Weight(struct network *LAYERS, int cl, float epsilon, float coupling){

  int i;
  float total_weight;

  if (cl == 2){

    LAYERS[0].fitness = epsilon;
    LAYERS[1].fitness = 1-epsilon;

    total_weight = FI_function(LAYERS[0].rho.act, coupling)*LAYERS[0].fitness
                  + FI_function(LAYERS[1].rho.act, coupling)*LAYERS[1].fitness;

    if (total_weight <= 0){total_weight = 1;} // This parameter is needed for the start if rho_act_i = 0

    LAYERS[0].weight = FI_function(LAYERS[0].rho.act, coupling)*LAYERS[0].fitness/(float)total_weight;
    LAYERS[1].weight = FI_function(LAYERS[1].rho.act, coupling)*LAYERS[1].fitness/(float)total_weight;

  }else {

    // printf("\r");
    // printf("   >>> FITNESS MODEL NOT APPLIED <<<\r");
    // printf("\r");

    total_weight = 0.0;

    for(i=0;i<cl;i++){
      total_weight += FI_function(LAYERS[i].rho.act, coupling);
    }

    if (total_weight <= 0){total_weight = 1.0;}

    for(i=0;i<cl;i++){
      LAYERS[i].weight = (FI_function(LAYERS[i].rho.act, coupling))/(float)total_weight;
    }
  }
}

void Reactivation(struct network *net){

  int j, index, index_inf;
  int infected, rand_ind;

  // which susceptible link becomes infected [I,S] -> [I,I]?
  rand_ind = rand()%(net->SIS.count.p2a);
  index_inf = net->SIS.pointEA_P[rand_ind] - 1;
  infected = net->SIS.EA[1][index_inf];
  net->SIS.NI[net->SIS.Nact] = infected;
  // printf("Reactivation of %d   !!!     \n",infected);
  // Which new neighbors become susceptible?
  index = get_index(infected,net->NODES,net->EN.N);
  net->SIS.Nact += 1;

  for(j=net->POINT.INI[index];j<=net->POINT.FIN[index];j++){

      speedUpdateActivation(net,j,infected);

  }

  net->SIS.NS[index] = 1;

}

void Deactivation(struct network *net){

  int j, index, index_rec, recovered;

  index_rec = rand()%(net->SIS.Nact);
  recovered = net->SIS.NI[index_rec];

  net->SIS.Nact -= 1;
  net->SIS.NI[index_rec] = net->SIS.NI[net->SIS.Nact];
  net->SIS.NI[net->SIS.Nact] = 0;

  index = get_index(recovered,net->NODES,net->EN.N);

  for(j=net->POINT.INI[index];j<=net->POINT.FIN[index];j++){

    speedUpdateDeactivation(net,j,recovered);

  }

}

void MMactivation(struct network *net){

  int index_sus, index_inf, infected;
  int j;

  index_sus = rand()%(net->SIS.Nsus);
  index_inf = net->SIS.pointNS[index_sus]-1;
  infected = net->NODES[index_inf];

  net->SIS.Nsus -= 1;
  net->SIS.pointNS[index_sus] = net->SIS.pointNS[net->SIS.Nsus];
  net->SIS.NI[net->SIS.Nact] = infected;
  net->SIS.Nact += 1;

  // printf("MMactivation of %d   !!!     \n",infected);

  for(j=net->POINT.INI[index_inf];j<=net->POINT.FIN[index_inf];j++){

      speedUpdateActivation(net,j,infected);

  }

  net->SIS.NS[index_inf] = 1; // by definition of NS only nodes still non active can be chosen
}

void Activation(struct network *net){

  int index_sus, index_pNS, index, infected;
  int rand_ind, j;

  rand_ind = rand()%(net->SIS.count.s2a);
  index_sus = net->SIS.pointEA_S[rand_ind] - 1;
  infected = net->SIS.EA[1][index_sus];
  // printf("\n");
  // printf("Activation  of %d   !!!     \n",infected);
  // printf("\n");
  net->SIS.NI[net->SIS.Nact] = infected;
  index = get_index(infected,net->NODES,net->EN.N);

  qsort(net->SIS.pointNS, net->SIS.Nsus, sizeof(int),cmpfunc_c);

  index_pNS = get_index(index+1,net->SIS.pointNS,net->SIS.Nsus);

  net->SIS.Nsus -= 1;
  net->SIS.pointNS[index_pNS] = net->SIS.pointNS[net->SIS.Nsus];
  net->SIS.Nact += 1;  // Which new neighbors become susceptible?

  for(j=net->POINT.INI[index];j<=net->POINT.FIN[index];j++){

        speedUpdateActivation(net,j,infected);

  }
  net->SIS.NS[index] = 1; // by definition of NS only nodes still non active can be chosen
}

float FI_function(float rho, float coupling){

  float fi;

  fi = pow(rho,coupling);

  return fi;
}

bool EpidemicDied(struct network *LAYERS, int cl){

  int i;
  bool dead;

  dead = true;
  i=0;

  while((i<cl)&& dead){
    dead = (LAYERS[i].SIS.Nact == 0);
    i += 1;
  }

  return dead;
}

bool Coex(struct network *LAYERS){

  return (LAYERS[0].rho.act > 0.05) && (LAYERS[1].rho.act > 0.05);

}

bool StationaryState(struct network *LAYERS, int cl, int bkt){

  int i,j;
  int SS[cl];
  int sstate;

  sstate = 0;
  for(j=0;j<cl;j++){

    SS[j] = true;
    i = 1;

    while ((i<bkt) && (SS[j])){
      SS[j] = (abs(LAYERS[j].SIS.Tss[i] - LAYERS[j].SIS.Tss[i-1]) <= 0.01*(float)LAYERS[j].EN.N); // 3% tolerance
      i += 1;
    }
    sstate += SS[j];
  }

  return (sstate == cl);

}

char* concat(const char *s1, const char *s2){

    char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

void WriteHeader(FILE* stability, float lambda, float epsilon, float sigma, float rho1_0, float rho2_0){

  fprintf(stability, "# lambda %f\n", lambda);
  fprintf(stability, "# epsilon0 %f\n", epsilon);
  fprintf(stability, "# sigma %f\n", sigma);
  fprintf(stability, "# rho1_0 %f\n", rho1_0);
  fprintf(stability, "# rho2_0 %f\n", rho2_0);
  fprintf(stability, "# \n");

}

struct competition ReadStabilityFile(char *filename){

  char str[100];
  char *line, *token;
  float rho1, rho2, rho1_0, rho2_0;
  float tol;
  bool extinction;
  int domCases = 0, coexCases = 0, tkvCases = 0;
  int count;
  struct competition Prob;

  filename=strcat(filename,".txt");
  tol = 0.01; count = 0;

  FILE *fp;
  fp = fopen(filename, "r");
  if (fp){
    line=fgets(str,100,fp);
    while (line != NULL){

      token = strtok(line," "); // Reading the first token

      // BEGINN READ HEADER
      while(strcmp(token,"#") == 0){ // If first token is # skip

        token = strtok(NULL," "); //reading following tokens
        if(strcmp(token,"lambda") == 0){
          token = strtok(NULL," ");
          Prob.lambda = atof(token);
        }else if(strcmp(token,"epsilon") == 0){
          token = strtok(NULL," ");
          Prob.epsilon = atof(token);
        }else if(strcmp(token,"sigma") == 0){
          token = strtok(NULL," ");
          Prob.sigma = atof(token);
        }else if(strcmp(token,"rho1_0") == 0){
          token = strtok(NULL," ");
          rho1_0 = atof(token);
          Prob.rho1_0 = rho1_0;
        }else if(strcmp(token,"rho2_0") == 0){
          token = strtok(NULL," ");
          rho2_0 = atof(token);
          Prob.rho2_0 = rho2_0;
        }

        line = fgets(str,100,fp); //Read next line
        token = strtok(line," ");
      }
      // END READ HEADER

      rho1 = atof(token); // Parsing the char pointer to int
      token = strtok(NULL," "); //reading following tokens
      rho2 = atof(token);  // Parsing the char pointer to int
      extinction = (abs(rho1-0) < tol) || (abs(rho2-0) < tol);  //extinction when one or the other net end with less than 1% of active nodes

      if(extinction){ // Dominance
        domCases += 1;
      }else{ // Coexistence
        coexCases += 1;
      }

      if(((rho1_0 > rho2_0) && (rho1 < rho2)) || ((rho1_0 < rho2_0) && (rho1 > rho2))){ //take over
        tkvCases += 1;
      }
      count += 1;
      line = fgets(str,100,fp); //Read next line
    }
  }
  else{printf( "Error!!! file not found !!\n");}


  Prob.dom = domCases/(float)count;
  Prob.coex = coexCases/(float)count;
  Prob.tkv = tkvCases/(float)count;

  return Prob;

}

void UpdatePointerValue(struct network *net,int node_i, int node_j, int newValue, int act){

  /*
    0) there has been an infection/recovery of, j iterates over its links (neighbors)
    1) Writes in position j the newValue, the position+1 in the matrix of edges EA (= 0 if link desactivated)
    2) Look for the infected/recovered in the neighbors of j (l) and writes new value
  */


    int index0_i, index0_j, index_i, index_j, n;
    bool susceptible, found;

    index0_i = get_index(node_i,net->NODES,net->EN.N); // = 1

    n = net->POINT.INI[index0_i]; found = false;
    while((!found)&&(n<=net->POINT.FIN[index0_i])){
      found = (net->LINKS[n] == node_j);
      if(found){index_j = n;}
      n += 1;
    }
    if((net->SIS.Nsus == 0)&&(net->SIS.count.s2a > 0)){
      printf(">>> WARNING!! Nsus = 0 but edges S2A = %d  \n",net->SIS.count.s2a);
    }
    ///////////////////////////////////////////////////////////////////////////////////7
    index0_j = get_index(node_j,net->NODES,net->EN.N); // = 1

    susceptible = (net->SIS.NS[index0_j] == 0)&&(net->SIS.Nsus > 0);
    n = net->POINT.INI[index0_j]; found = false;
    while((!found)&&(n<=net->POINT.FIN[index0_j])){
      found = (net->LINKS[n] == node_i);
      if(found){index_i = n;}
      n += 1;
    }

    net->SIS.pointerEA[index_i] = net->SIS.pointerEA[index_j] = newValue;

    if(act == 0){ //addition
      net->SIS.EA2pnt[net->SIS.Eact-1] = index_i;
      if(susceptible){ //active edge with susceptible node
        net->SIS.pointEA_SP[index_i] = net->SIS.pointEA_SP[index_j] = 1;
      }else{ //active edge with passive node
        net->SIS.pointEA_SP[index_i] = net->SIS.pointEA_SP[index_j] = 2;
      }
    }else if(act == 1){ //substraction
      net->SIS.pointEA_SP[index_i] = net->SIS.pointEA_SP[index_j] = 0;
    }else if(act == 2){ //relocation
      net->SIS.EA2pnt[newValue-1] = index_i;
    }
}
    //  net->SIS.count.s2a -= 1;

void FillpLINKS(struct network *net){

  int i,j,k;
  int actualnode, index;
  bool found;

  k=0;

  for(i=0;i<net->EN.N;i++){
    actualnode = net->NODES[i];
    for(j=net->POINT.INI[i];j<=net->POINT.FIN[i];j++){
      index = get_index(net->LINKS[j], net->NODES, net->EN.N);
      found = false;
      k=net->POINT.INI[index];
      while(k<=net->POINT.FIN[index] && !found){
        found = (net->LINKS[k] == actualnode);
        k += 1;
      }
      net->pLINKS[j] = k-1; //because the last += 1 for k
    }
  }

}

void speedUpdateActivation(struct network *net, int j, int infected){

  int i, indexs2a, indexp2a, indexAddEA, neighbor, indexNeighborNodes;
  int index2rmv, InfWasSus, indexInfectedNodes;
  int i_reloc, j_reloc;
  int lastInpEA_S, lastInpEA_P;
  int indexInpEA_S, indexInpEA_P;
  int i_LastInpEA, j_LastInpEA;
  int ind2rlc_InpEA_S, ind2rmv_InpEA_S, indLast_InpEA_S;
  int ind2rlc_InpEA_P, ind2rmv_InpEA_P, indLast_InpEA_P;
  bool lastisSusceptible;
  // int k;

  neighbor = net->LINKS[j];
  i = net->pLINKS[j];

  indexNeighborNodes = get_index(neighbor, net->NODES, net->EN.N);
  indexInfectedNodes = get_index(infected, net->NODES, net->EN.N);

  // index_inf depends on the type of activation:
  //   - For Activation is the index in pointEA_S
  //   - For Reactivation is the index in pointEA_P
  //   - For MM Activation is the index in pointNS


  if(net->SIS.pointerEA[j] == 0){

    // printf("    ...addition of an edge\n");

      // Activation of an edge:
      // - If there is an activation and pEA[j] was 0 means that neighbor is inactive.
      //   In such case a new active edge will be added.
      // - The addition is at the end of EA matrix, index = Eact
      // - Then the pointerEA is updated with the new position in EA (+1)

    indexAddEA = net->SIS.Eact;
    net->SIS.EA[0][indexAddEA] = infected;
    net->SIS.EA[1][indexAddEA] = neighbor;
    net->SIS.pointerEA[i] = net->SIS.pointerEA[j] = indexAddEA + 1;

    if(net->SIS.NS[indexNeighborNodes] == 0){ // neighbor is susceptible

      net->SIS.pointEA_SP[i] = net->SIS.pointEA_SP[j] = 1;
      indexs2a = net->SIS.count.s2a;
      net->SIS.pointEA_S[indexs2a] = indexAddEA + 1;

      net->SIS.ppS[i] = net->SIS.ppS[j] = indexs2a + 1;

      net->SIS.count.s2a += 1;

    }else if(net->SIS.NS[indexNeighborNodes] == 1){ // neighbor is passive

      net->SIS.pointEA_SP[i] = net->SIS.pointEA_SP[j] = 2;
      indexp2a = net->SIS.count.p2a;
      net->SIS.pointEA_P[indexp2a] = indexAddEA + 1;

      net->SIS.ppP[i] = net->SIS.ppP[j] = indexp2a + 1;

      net->SIS.count.p2a += 1;

    }

    net->SIS.EA2pnt[indexAddEA] = i;
    net->SIS.Eact += 1;

  }else if(net->SIS.pointerEA[j] != 0){

    // printf("    ...removal of an edge\n");

      // Deactivation of an edge:
      //  - If the activated node was part of an active edge it become inactive
      //  - This means remove it from EA. It is substituted by the last edge
      //    of the matrix
      //  - Then the last row of the matrix is deleted
      //  - The relocation also needs to be updated in the pointer Vectors
      //  - Now it is necessary to check whether the activated node was passive or susceptible
      //  - Last link in EA can be susceptible or passive

    index2rmv = net->SIS.pointerEA[j] - 1;
    net->SIS.EA[0][index2rmv] = net->SIS.EA[0][net->SIS.Eact-1];
    net->SIS.EA[1][index2rmv] = net->SIS.EA[1][net->SIS.Eact-1];
    net->SIS.EA[0][net->SIS.Eact - 1] = 0;
    net->SIS.EA[1][net->SIS.Eact - 1] = 0;
    net->SIS.pointerEA[i] = net->SIS.pointerEA[j] = 0;

    i_reloc = net->SIS.EA2pnt[net->SIS.Eact - 1];
    j_reloc = net->pLINKS[i_reloc];

    if(index2rmv != net->SIS.Eact - 1){
      net->SIS.pointerEA[i_reloc] = net->SIS.pointerEA[j_reloc] = index2rmv + 1;
      net->SIS.EA2pnt[index2rmv] = i_reloc;
    }

    lastisSusceptible = (net->SIS.pointEA_SP[i_reloc] == 1);
    net->SIS.EA2pnt[net->SIS.Eact - 1] = 0;

    InfWasSus = (net->SIS.NS[indexInfectedNodes] == 0); // infected was susceptible

    if(InfWasSus && lastisSusceptible){

      // Remove the max value from pointEA_S
      // qsort(net->SIS.pointEA_S, net->SIS.count.s2a, sizeof(int), cmpfunc_c);
      // net->SIS.pointEA_S[net->SIS.count.s2a - 1] = 0;

      ind2rmv_InpEA_S = net->SIS.ppS[i] - 1;
      ind2rlc_InpEA_S = net->SIS.ppS[i_reloc] - 1;
      indLast_InpEA_S = net->SIS.count.s2a - 1;

      lastInpEA_S = net->SIS.pointEA_S[indLast_InpEA_S] - 1;

      i_LastInpEA = net->SIS.EA2pnt[lastInpEA_S];
      j_LastInpEA = net->pLINKS[i_LastInpEA];

      net->SIS.pointEA_S[ind2rlc_InpEA_S] = net->SIS.pointEA_S[ind2rmv_InpEA_S];
      net->SIS.pointEA_S[ind2rmv_InpEA_S] = net->SIS.pointEA_S[indLast_InpEA_S];
      net->SIS.pointEA_S[indLast_InpEA_S] = 0;

      if(indLast_InpEA_S != ind2rlc_InpEA_S){
        net->SIS.ppS[i_LastInpEA] = net->SIS.ppS[j_LastInpEA] = ind2rmv_InpEA_S + 1;
        net->SIS.ppS[i_reloc] = net->SIS.ppS[j_reloc] = ind2rlc_InpEA_S + 1;
      }else{
        net->SIS.ppS[i_reloc] = net->SIS.ppS[j_reloc] = ind2rmv_InpEA_S + 1;
      }

      net->SIS.ppS[i] = net->SIS.ppS[j] = 0;

      net->SIS.count.s2a -= 1;

    }else if(InfWasSus && !lastisSusceptible){ // infected was passive

      // 1) Include in pointEA_P the new index and remove the max value
      // 2) Remover from pointEA_S the new index

      // qsort(net->SIS.pointEA_P, net->SIS.count.p2a, sizeof(int), cmpfunc_c);
      // net->SIS.pointEA_P[net->SIS.count.p2a - 1] = index2rmv + 1;

      indexInpEA_P = net->SIS.ppP[i_reloc] - 1;
      net->SIS.pointEA_P[indexInpEA_P] = index2rmv + 1;
      net->SIS.ppP[i_reloc] = net->SIS.ppP[j_reloc] = indexInpEA_P + 1;
      net->SIS.ppP[i] = net->SIS.ppP[j] = 0;

      // qsort(net->SIS.pointEA_S, net->SIS.count.s2a, sizeof(int), cmpfunc_c);
      // indd = get_index(index2rmv+1, net->SIS.pointEA_S, net->SIS.count.s2a);
      // net->SIS.pointEA_S[indd] = net->SIS.pointEA_S[net->SIS.count.s2a - 1];
      // net->SIS.pointEA_S[net->SIS.count.s2a - 1] = 0;

      indexInpEA_S = net->SIS.ppS[i] - 1;
      lastInpEA_S = net->SIS.pointEA_S[net->SIS.count.s2a - 1];

      i_LastInpEA = net->SIS.EA2pnt[lastInpEA_S-1];
      j_LastInpEA = net->pLINKS[i_LastInpEA];
      net->SIS.ppS[i_LastInpEA] = net->SIS.ppS[j_LastInpEA] = net->SIS.ppS[i];

      net->SIS.pointEA_S[indexInpEA_S] = lastInpEA_S;
      net->SIS.pointEA_S[net->SIS.count.s2a - 1] = 0;
      net->SIS.ppS[i] = net->SIS.ppS[j] = 0;
      net->SIS.ppS[i_reloc] = net->SIS.ppS[j_reloc] = 0;

      net->SIS.count.s2a -= 1;


    }else if(!InfWasSus && lastisSusceptible){

      // 1) Include in pointEA_S the new index and remove the max value
      // 2) Remover from pointEA_P the new index

      // qsort(net->SIS.pointEA_S, net->SIS.count.s2a, sizeof(int), cmpfunc_c);
      // net->SIS.pointEA_S[net->SIS.count.s2a - 1] = index2rmv + 1;

      indexInpEA_S = net->SIS.ppS[i_reloc] - 1;
      net->SIS.pointEA_S[indexInpEA_S] = index2rmv + 1;
      net->SIS.ppS[i_reloc] = net->SIS.ppS[j_reloc] = indexInpEA_S + 1;
      net->SIS.ppS[i] = net->SIS.ppS[j] = 0;

      // qsort(net->SIS.pointEA_P, net->SIS.count.p2a, sizeof(int), cmpfunc_c);
      // indd = get_index(index2rmv+1, net->SIS.pointEA_P, net->SIS.count.p2a);
      // net->SIS.pointEA_P[indd] = net->SIS.pointEA_P[net->SIS.count.p2a - 1];
      // net->SIS.pointEA_P[net->SIS.count.p2a - 1] = 0;

      indexInpEA_P = net->SIS.ppP[i] - 1;
      lastInpEA_P = net->SIS.pointEA_P[net->SIS.count.p2a - 1];

      i_LastInpEA = net->SIS.EA2pnt[lastInpEA_P-1];
      j_LastInpEA = net->pLINKS[i_LastInpEA];
      net->SIS.ppP[i_LastInpEA] = net->SIS.ppP[j_LastInpEA] = net->SIS.ppP[i];

      net->SIS.pointEA_P[indexInpEA_P] = lastInpEA_P;
      net->SIS.pointEA_P[net->SIS.count.p2a - 1] = 0;
      net->SIS.ppP[i] = net->SIS.ppP[j] = 0;
      net->SIS.ppP[i_reloc] = net->SIS.ppP[j_reloc] = 0;

      net->SIS.count.p2a -= 1;


    }else if(!InfWasSus && !lastisSusceptible){

      // Remove the max value from pointEA_P
      // qsort(net->SIS.pointEA_P, net->SIS.count.p2a, sizeof(int), cmpfunc_c);
      // net->SIS.pointEA_P[net->SIS.count.p2a - 1] = 0;

      ind2rmv_InpEA_P = net->SIS.ppP[i] - 1;
      ind2rlc_InpEA_P = net->SIS.ppP[i_reloc] - 1;
      indLast_InpEA_P = net->SIS.count.p2a - 1;

      lastInpEA_P = net->SIS.pointEA_P[indLast_InpEA_P] - 1;

      i_LastInpEA = net->SIS.EA2pnt[lastInpEA_P];
      j_LastInpEA = net->pLINKS[i_LastInpEA];

      net->SIS.pointEA_P[ind2rlc_InpEA_P] = net->SIS.pointEA_P[ind2rmv_InpEA_P];
      net->SIS.pointEA_P[ind2rmv_InpEA_P] = net->SIS.pointEA_P[indLast_InpEA_P];
      net->SIS.pointEA_P[indLast_InpEA_P] = 0;

      if(indLast_InpEA_P != ind2rlc_InpEA_P){
        net->SIS.ppP[i_LastInpEA] = net->SIS.ppP[j_LastInpEA] = ind2rmv_InpEA_P + 1;
        net->SIS.ppP[i_reloc] = net->SIS.ppP[j_reloc] = ind2rlc_InpEA_P + 1;
      }else{
        net->SIS.ppP[i_reloc] = net->SIS.ppP[j_reloc] = ind2rmv_InpEA_P + 1;
      }

      net->SIS.ppP[i] = net->SIS.ppP[j] = 0;

      net->SIS.count.p2a -= 1;

    }

    net->SIS.pointEA_SP[i] = net->SIS.pointEA_SP[j] = 0;
    net->SIS.Eact -= 1;

    // DebugPP(net);

  }


}

void speedUpdateDeactivation(struct network *net, int j, int recovered){

  int i, indexp2a, indexAddEA, neighbor;
  int index2rmv, indexNeighborNodes, NeighIsSus;
  int lastInpEA_S, lastInpEA_P;
  int i_reloc, j_reloc, indexInpEA_P, indexInpEA_S;
  int i_LastInpEA, j_LastInpEA;
  int ind2rlc_InpEA_S, ind2rmv_InpEA_S, indLast_InpEA_S;
  int ind2rlc_InpEA_P, ind2rmv_InpEA_P, indLast_InpEA_P;
  bool lastisSusceptible;
  // int k;

  neighbor = net->LINKS[j];
  i = net->pLINKS[j];

  indexNeighborNodes = get_index(neighbor, net->NODES, net->EN.N);

  if(net->SIS.pointerEA[j] == 0){

    // printf("    ...addition of an edge\n");

    // Activation of an edge:
    // - If there is an activation and pEA[j] was 0 means that neighbor is active.
    //   In such case a new active edge will be added.
    // - The recovered will also affect P

    indexAddEA = net->SIS.Eact;
    net->SIS.EA[0][indexAddEA] = neighbor;
    net->SIS.EA[1][indexAddEA] = recovered;
    net->SIS.pointerEA[i] = net->SIS.pointerEA[j] = indexAddEA + 1;

    net->SIS.pointEA_SP[i] = net->SIS.pointEA_SP[j] = 2;
    indexp2a = net->SIS.count.p2a;
    net->SIS.pointEA_P[indexp2a] = indexAddEA + 1;

    net->SIS.ppP[i] = net->SIS.ppP[j] = indexp2a + 1;

    net->SIS.count.p2a += 1;

    net->SIS.EA2pnt[indexAddEA] = i;
    net->SIS.Eact += 1;

  }else if(net->SIS.pointerEA[j] != 0){

    // printf("    ...removal of an edge\n");
    // Deactivation of an edge:
    //  - If the deactivated node was part of an active edge it become inactive

    index2rmv = net->SIS.pointerEA[j] - 1;
    net->SIS.EA[0][index2rmv] = net->SIS.EA[0][net->SIS.Eact-1];
    net->SIS.EA[1][index2rmv] = net->SIS.EA[1][net->SIS.Eact-1];
    net->SIS.EA[0][net->SIS.Eact - 1] = 0;
    net->SIS.EA[1][net->SIS.Eact - 1] = 0;
    net->SIS.pointerEA[i] = net->SIS.pointerEA[j] = 0;

    i_reloc = net->SIS.EA2pnt[net->SIS.Eact - 1];
    j_reloc = net->pLINKS[i_reloc];

    if(index2rmv != net->SIS.Eact - 1){
      net->SIS.pointerEA[i_reloc] = net->SIS.pointerEA[j_reloc] = index2rmv + 1;
      net->SIS.EA2pnt[index2rmv] = i_reloc;
    }

    lastisSusceptible = (net->SIS.pointEA_SP[i_reloc] == 1);
    net->SIS.EA2pnt[net->SIS.Eact - 1] = 0;

    NeighIsSus = (net->SIS.NS[indexNeighborNodes] == 0); // infected was susceptible

    if(NeighIsSus && lastisSusceptible){

      // Remove the max value from pointEA_S
      // qsort(net->SIS.pointEA_S, net->SIS.count.s2a, sizeof(int), cmpfunc_c);
      // net->SIS.pointEA_S[net->SIS.count.s2a - 1] = 0;

      ind2rmv_InpEA_S = net->SIS.ppS[i] - 1;
      ind2rlc_InpEA_S = net->SIS.ppS[i_reloc] - 1;
      indLast_InpEA_S = net->SIS.count.s2a - 1;

      lastInpEA_S = net->SIS.pointEA_S[indLast_InpEA_S] - 1;

      i_LastInpEA = net->SIS.EA2pnt[lastInpEA_S];
      j_LastInpEA = net->pLINKS[i_LastInpEA];

      net->SIS.pointEA_S[ind2rlc_InpEA_S] = net->SIS.pointEA_S[ind2rmv_InpEA_S];
      net->SIS.pointEA_S[ind2rmv_InpEA_S] = net->SIS.pointEA_S[indLast_InpEA_S];
      net->SIS.pointEA_S[indLast_InpEA_S] = 0;

      if(indLast_InpEA_S != ind2rlc_InpEA_S){
        net->SIS.ppS[i_LastInpEA] = net->SIS.ppS[j_LastInpEA] = ind2rmv_InpEA_S + 1;
        net->SIS.ppS[i_reloc] = net->SIS.ppS[j_reloc] = ind2rlc_InpEA_S + 1;
      }else{
        net->SIS.ppS[i_reloc] = net->SIS.ppS[j_reloc] = ind2rmv_InpEA_S + 1;
      }

      net->SIS.ppS[i] = net->SIS.ppS[j] = 0;

      net->SIS.count.s2a -= 1;

    }else if(NeighIsSus && !lastisSusceptible){ // infected was passive

      // 1) Include in pointEA_P the new index and remove the max value
      // 2) Remover from pointEA_S the new index

      // qsort(net->SIS.pointEA_P, net->SIS.count.p2a, sizeof(int), cmpfunc_c);
      // net->SIS.pointEA_P[net->SIS.count.p2a - 1] = index2rmv + 1;

      indexInpEA_P = net->SIS.ppP[i_reloc] - 1; //Index for the last item in EA, which was passive
      net->SIS.pointEA_P[indexInpEA_P] = index2rmv + 1; // I substitute this last item for the new index in EA
      net->SIS.ppP[i_reloc] = net->SIS.ppP[j_reloc] = indexInpEA_P + 1;
      net->SIS.ppP[i] = net->SIS.ppP[j] = 0;

      // qsort(net->SIS.pointEA_S, net->SIS.count.s2a, sizeof(int), cmpfunc_c);
      // indd = get_index(index2rmv+1, net->SIS.pointEA_S, net->SIS.count.s2a);
      // net->SIS.pointEA_S[indd] = net->SIS.pointEA_S[net->SIS.count.s2a - 1];
      // net->SIS.pointEA_S[net->SIS.count.s2a - 1] = 0;

      indexInpEA_S = net->SIS.ppS[i] - 1;
      lastInpEA_S = net->SIS.pointEA_S[net->SIS.count.s2a - 1];

      i_LastInpEA = net->SIS.EA2pnt[lastInpEA_S-1];
      j_LastInpEA = net->pLINKS[i_LastInpEA];
      net->SIS.ppS[i_LastInpEA] = net->SIS.ppS[j_LastInpEA] = net->SIS.ppS[i];

      net->SIS.pointEA_S[indexInpEA_S] = lastInpEA_S;
      net->SIS.pointEA_S[net->SIS.count.s2a - 1] = 0;
      net->SIS.ppS[i] = net->SIS.ppS[j] = 0;
      net->SIS.ppS[i_reloc] = net->SIS.ppS[j_reloc] = 0;

      net->SIS.count.s2a -= 1;

    }else if(!NeighIsSus && lastisSusceptible){

      // 1) Include in pointEA_S the new index and remove the max value
      // 2) Remove from pointEA_P the new index

      // qsort(net->SIS.pointEA_S, net->SIS.count.s2a, sizeof(int), cmpfunc_c);
      // net->SIS.pointEA_S[net->SIS.count.s2a - 1] = index2rmv + 1;

      indexInpEA_S = net->SIS.ppS[i_reloc] - 1;
      net->SIS.pointEA_S[indexInpEA_S] = index2rmv + 1;
      net->SIS.ppS[i_reloc] = net->SIS.ppS[j_reloc] = indexInpEA_S + 1;
      net->SIS.ppS[i] = net->SIS.ppS[j] = 0;

      // qsort(net->SIS.pointEA_P, net->SIS.count.p2a, sizeof(int), cmpfunc_c);
      // indd = get_index(index2rmv+1, net->SIS.pointEA_P, net->SIS.count.p2a);
      // net->SIS.pointEA_P[indd] = net->SIS.pointEA_P[net->SIS.count.p2a - 1];
      // net->SIS.pointEA_P[net->SIS.count.p2a - 1] = 0;

      indexInpEA_P = net->SIS.ppP[i] - 1;
      lastInpEA_P = net->SIS.pointEA_P[net->SIS.count.p2a - 1];

      i_LastInpEA = net->SIS.EA2pnt[lastInpEA_P-1];
      j_LastInpEA = net->pLINKS[i_LastInpEA];
      net->SIS.ppP[i_LastInpEA] = net->SIS.ppP[j_LastInpEA] = net->SIS.ppP[i];

      net->SIS.pointEA_P[indexInpEA_P] = lastInpEA_P;
      net->SIS.pointEA_P[net->SIS.count.p2a - 1] = 0;
      net->SIS.ppP[i] = net->SIS.ppP[j] = 0;
      net->SIS.ppP[i_reloc] = net->SIS.ppP[j_reloc] = 0;

      net->SIS.count.p2a -= 1;

    }else if(!NeighIsSus && !lastisSusceptible){

      // Remove the max value from pointEA_P
      // qsort(net->SIS.pointEA_P, net->SIS.count.p2a, sizeof(int), cmpfunc_c);
      // net->SIS.pointEA_P[net->SIS.count.p2a - 1] = 0;

      ind2rmv_InpEA_P = net->SIS.ppP[i] - 1;
      ind2rlc_InpEA_P = net->SIS.ppP[i_reloc] - 1;
      indLast_InpEA_P = net->SIS.count.p2a - 1;

      lastInpEA_P = net->SIS.pointEA_P[indLast_InpEA_P] - 1;

      i_LastInpEA = net->SIS.EA2pnt[lastInpEA_P];
      j_LastInpEA = net->pLINKS[i_LastInpEA];

      net->SIS.pointEA_P[ind2rlc_InpEA_P] = net->SIS.pointEA_P[ind2rmv_InpEA_P];
      net->SIS.pointEA_P[ind2rmv_InpEA_P] = net->SIS.pointEA_P[indLast_InpEA_P];
      net->SIS.pointEA_P[indLast_InpEA_P] = 0;

      if(indLast_InpEA_P != ind2rlc_InpEA_P){
        net->SIS.ppP[i_LastInpEA] = net->SIS.ppP[j_LastInpEA] = ind2rmv_InpEA_P + 1;
        net->SIS.ppP[i_reloc] = net->SIS.ppP[j_reloc] = ind2rlc_InpEA_P + 1;
      }else{
        net->SIS.ppP[i_reloc] = net->SIS.ppP[j_reloc] = ind2rmv_InpEA_P + 1;
      }

      net->SIS.ppP[i] = net->SIS.ppP[j] = 0;

      net->SIS.count.p2a -= 1;

    }

    net->SIS.pointEA_SP[i] = net->SIS.pointEA_SP[j] = 0;
    net->SIS.Eact -= 1;

  }

  // DebugPP(net);

}

void DebugPointers(struct network *net){

  int i,j;
  int count, max, countSus, index_sus, totalEactSum;
  int *positions;

  count = 0;
  countSus = 0;
  max = 0;
  for(i=0;i<net->EN.N;i++){
    for(j=net->POINT.INI[i];j<=net->POINT.FIN[i];j++){
      count += net->SIS.pointerEA[j];
      if(max < net->SIS.pointerEA[j]){
        max = net->SIS.pointerEA[j];
      }
    }
  }

  totalEactSum = 0;
  for(i=1;i<=net->SIS.Eact;i++){
    totalEactSum += i;
  }
  if(max != net->SIS.Eact){
    printf("\n");
    printf("------------------------------------\n");
    printf(" ERROR #1: \n");
    printf("  Max value of pEA doesn't match Eact\n");
    printf("  pEA (Eact = %d) : ", net->SIS.Eact);
    for(i=0;i<net->lenLINKS;i++){
      if( net->SIS.pointerEA[i] != 0){printf("%d ", net->SIS.pointerEA[i]);}
    }printf("\n");
    printf("------------------------------------\n");
    printf("\n");
    exit(1);
  }

  if(count != 2*totalEactSum){
    printf("\n");
    printf("------------------------------------\n");
    printf(" ERROR #2: \n");
    printf("  pEA values not twice each\n");
    printf("  pEA (Eact = %d) : ", net->SIS.Eact);
    for(i=0;i<net->lenLINKS;i++){
      if( net->SIS.pointerEA[i] != 0){printf("%d ", net->SIS.pointerEA[i]);}
    }printf("\n");
    printf("------------------------------------\n");
    printf("\n");
    exit(1);
  }


  positions = malloc(sizeof(int)*net->SIS.Eact);
  for(i=0;i<net->SIS.count.s2a;i++){
    positions[i] = net->SIS.pointEA_S[i];
  }
  for(i=0;i<net->SIS.count.p2a;i++){
    positions[i+net->SIS.count.s2a] = net->SIS.pointEA_P[i];
  }
  qsort(positions, net->SIS.Eact, sizeof(int), cmpfunc_c);
  for(i=0;i<net->SIS.Eact;i++){
    if(positions[i] != i+1){
      printf("\n");
      printf("------------------------------------\n");
      printf(" ERROR #3: \n");
      printf("  Unbalanced pEA_S, pEA_P\n");
      printf("  Eact = %d \n", net->SIS.Eact);
      printf("  pEA_S : ");
      for(i=0;i<net->SIS.count.s2a;i++){
        printf("%d ", net->SIS.pointEA_S[i]);
      }printf("\n");
      printf("  pEA_P : ");
      for(i=0;i<net->SIS.count.p2a;i++){
        printf("%d ", net->SIS.pointEA_P[i]);
      }printf("\n");
      printf("------------------------------------\n");
      printf("\n");
      exit(1);
    }
  }
  free(positions);

  for(i=0;i<net->EN.N;i++){
    countSus += net->SIS.NS[i];
    if(net->SIS.NS[i] == 0){
      qsort(net->SIS.pointNS, net->SIS.Nsus, sizeof(int), cmpfunc_c);
      index_sus = get_index(i+1,net->SIS.pointNS,net->SIS.Nsus);
      if(index_sus > net->SIS.Nsus){
        printf("\n");
        printf("------------------------------------\n");
        printf(" ERROR #4: \n");
        printf("  pNS not working\n");
        printf("------------------------------------\n");
        printf("\n");
        exit(1);
      }
    }
  }

  if(countSus != net->EN.N - net->SIS.Nsus){
    printf("\n");
    printf("------------------------------------\n");
    printf(" ERROR #5: \n");
    printf("  NS not working\n");
    printf("  NS (Nsus = %d) : ", net->SIS.Nsus);
    printf("  countSus: %d \n", countSus);
    for(i=0;i<net->EN.N;i++){
      printf("%d ", net->SIS.NS[i]);
    }printf("\n");
    printf("  NODES: ");
    for(i=0;i<net->EN.N;i++){
      printf("%d ", net->NODES[i]);
    }printf("\n");
    printf("------------------------------------\n");
    printf("\n");
    exit(1);
  }
}

void DebugPP(struct network *net){

  int i,j;
  int maxS,maxP;
  int count;

  maxS = maxP = 0;

  for(i=1;i<=net->SIS.count.s2a;i++){
    count = 0;
    for(j=0;j<net->lenLINKS;j++){
      if(net->SIS.ppS[j] == i){ count += 1; }
      if(maxS < net->SIS.ppS[j]){
        maxS = net->SIS.ppS[j];
      }
    }
  }
  if((net->SIS.count.s2a != 0) && ((count != 2) || (maxS != net->SIS.count.s2a))){
    printf("Error in ppS; \n");
    printf("  pEA_S : ");
    for(i=0;i<net->SIS.count.s2a;i++){
      printf("%d ", net->SIS.pointEA_S[i]);
    }printf("\n");
    printf("  ppS : ");
    for(i=0;i<net->lenLINKS;i++){
      if(net->SIS.ppS[i] != 0){printf("%d ", net->SIS.ppS[i]);}
    }printf("\n");

    exit(1);
  }

  for(i=1;i<=net->SIS.count.p2a;i++){
    count = 0;
    for(j=0;j<net->lenLINKS;j++){
      if(net->SIS.ppP[j] == i){ count += 1; }
      if(maxP < net->SIS.ppP[j]){
        maxP = net->SIS.ppP[j];
      }
    }

  }

  if((net->SIS.count.p2a != 0) && ((count != 2) || (maxP != net->SIS.count.p2a))){
    printf("Error in ppP: \n");
    printf("  pEA_P : ");
    for(i=0;i<net->SIS.count.p2a;i++){
      printf("%d ", net->SIS.pointEA_P[i]);
    }printf("\n");
    printf("  ppP : ");
    for(i=0;i<net->lenLINKS;i++){
      if(net->SIS.ppP[i] != 0){printf("%d ", net->SIS.ppP[i]);}
    }printf("\n");

    exit(1);

  }


}
