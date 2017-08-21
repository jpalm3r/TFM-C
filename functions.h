

struct result {
  int E;
  int N;
  int* NODES;
  int maxNODE;
};

struct pointers {
  int* INI;
  int* FIN;
};

struct state {
  int s2a;
  int p2a;
};

struct competition{
  float lambda;
  float sigma;
  float epsilon;
  float coex;
  float dom;
  float tkv;
  float rho1_0;
  float rho2_0;

};

struct gllsp_coef {
  float react;
  float mmact;
  float act;
  float deact;
  float TOTAL;
};

struct activation {
  int Nact;
  int Nsus;
  int Eact;
  int* NS;
  int* pointNS;
  int* NI;
  int** EA;
  int* pointerEA;
  int* pointEA_S;
  int* pointEA_P;
  int* pointEA_SP;
  int* ppS;
  int* ppP;
  int* EA2pnt;
  int* Tss;
  struct gllsp_coef a_i;
  struct state count;
};

struct degree{
  float* P;
  float* Pc;
};

struct density {
  float act;
  float pas;
  float susc;
};

struct evolution {
  float** act;
  float* act_avg;
  float** pas;
  float* act_pas;
  float* real_time;
};

struct network {
  struct result EN;
  struct degree DEG;
  struct activation SIS;
  struct pointers POINT;
  struct density rho;
  struct evolution EVOL;
  int* NODES;
  int* LINKS;
  int* pLINKS;
  int* DEGREE;
  float weight;
  float lambda;
  float mu;
  int loops;
  int zeros;
  int repetitions;
  int lenLINKS;
  float fitness;
};


/*****************************************************************************/

bool inVECTOR(int k, int size, int* VECTOR);
struct result get_E_N(char *filename);
void get_pointers(struct network *net, char* filename, int* DEGREE_0);
int cmpfunc_c(const void * a, const void * b);
int get_index(int ind, int* VECTOR, int len);
void last_link(struct network *net, char* filename);
bool link_is_once(int* link, int** MATRIX, int len);
bool check_file(char* filename,struct network *net);
void OverwritePointer(struct network *net, int af,int lenLINKS,int tmp,int neighbor);
void Event_Happens(struct network *LAYERS, float inf, float* COEFFS);
void initial_population(struct network *net, float fraction);
void WriteNet_Act(struct network net, char* directory, int filecount);
void WriteNet_Total(struct network net, char* directory, int filecount);
struct state ComputeEs2A(struct network *net);
//void Update_pointNS(struct network *net, int cl);
void Update_pointEA(struct network *LAYERS, int cl);
float ComputeMean(float* VECTOR,int len);
void ResetVectors(struct network *LAYERS, struct network *LAYERS_0, int cl);
void InitializeNet(struct network *net, char* file, char* directory);
void Compute_GllspCoeffs(struct network *LAYERS, int cl, float sigma, float* COEFFS);
void RandomActivation(struct network *LAYERS);
void Compute_Density(struct network *LAYERS, int cl);
void Update_ActRate(struct network *LAYERS, int cl, float lambda, float mu);
void Update_Weight(struct network *LAYERS, int cl, float epsilon, float coupling);
void Reactivation(struct network *net);
void Deactivation(struct network *net);
void MMactivation(struct network *net);
void Activation(struct network *net);
float FI_function(float rho, float coupling);
bool EpidemicDied(struct network *LAYERS, int cl);
bool StationaryState(struct network *LAYERS, int cl, int bkt);
char* concat(const char *s1, const char *s2);
void WriteHeader(FILE* stabiltiy, float lambda, float epsilon, float sigma, float rho1_0, float rho2_0);
struct competition ReadStabilityFile(char *filename);
bool Coex(struct network *LAYERS);
bool CoexDied(struct network *LAYERS, int cl);
void speedUpdateActivation(struct network *net, int j, int infected);
void speedUpdateDeactivation(struct network *net, int j, int recovered);
void FillpLINKS(struct network *net);
void DebugPointers(struct network *net);
void DebugPP(struct network *net);
