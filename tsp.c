/*****************************************************************************
  A template program for 2-dimensional euclidean symmetric TSP solver. 
  Subroutines to read instance data and compute the objective value of a given 
  tour (solution) are included. The one to output the computed tour in the
  TSPLIB format is also included. 

  The URL of TSPLIB is:
       http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/

  NOTE: Indices of nodes range from 0 to n-1 in this program,
        while it does from 1 to n in "README.md" and the data files of
        instances and of tours. 

  If you would like to use various parameters, it might be useful to modify
  the definition of struct "Param" and mimic the way the default value of
  "timelim" is given and how its value is input from the command line.
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "cpu_time.c"


/***** constants *************************************************************/
#define MAX_STR    1024

/***** macros ****************************************************************/
#define dist(k,l) ( (int)( sqrt( (tspdata->x[k]-tspdata->x[l])*(tspdata->x[k]-tspdata->x[l]) + (tspdata->y[k]-tspdata->y[l])*(tspdata->y[k]-tspdata->y[l]) ) + 0.5 ) )

/***** default values of parameters ******************************************/
#define TIMELIM    300 /* the time limit for the algorithm in seconds */
#define GIVESOL    0   /* 1: input a solution; 0: do not give a solution */
#define OUTFORMAT  2   /* 0: do not output the tour;
			  1: output the computed tour in TSPLIB format;
			  2: output the computed tour in TSP_VIEW format */
#define TOURFILE   "result.tour"   /* the output file of computed tour */

#define POPULATION 20  /* population of gene */


typedef struct {
  int    timelim;              /* the time limit for the algorithm in secs. */
  int    givesol;              /* give a solution (1) or not (0) */
  int    outformat;            /* 1: output the computed tour in TSPLIB format;
				  0: do not output it */
  char   tourfile[MAX_STR];    /* the output file of computed tour */
  /* NEVER MODIFY THE ABOVE VARIABLES.  */
  /* You can add more components below. */

} Param;                /* parameters */


typedef struct{
  char     name[MAX_STR];         /* name of the instance */
  int      n;                     /* number of nodes */
  double   *x;                    /* x-coordinates of nodes */
  double   *y;                    /* y-coordinates of nodes */
  int      min_node_num;          /* minimum number of nodes the solution contains */
} TSPdata;              /* data of TSP instance */

typedef struct {
  double        timebrid;       /* the time before reading the instance data */
  double        starttime;      /* the time the search started */
  double        endtime;        /* the time the search ended */
  int           *bestsol;       /* the best solution found so far */
  /* NEVER MODIFY THE ABOVE FOUR VARIABLES. */
  /* You can add more components below. */

} Vdata;                /* various data often necessary during the search */

/************************ declaration of functions ***************************/
FILE *open_file( char *fname, char *mode );
void *malloc_e( size_t size );

void copy_parameters( int argc, char *argv[], Param *param );
void prepare_memory( TSPdata *tspdata, Vdata *vdata );
void read_header( FILE *in, TSPdata *tspdata );
void read_tspfile( FILE *in, TSPdata *tspdata, Vdata *vdata );
void read_tourfile( FILE *in, TSPdata *tspdata, int *tour );
void output_tour( FILE *out, TSPdata *tspdata, int *tour );
void output_tour_for_tsp_view( FILE *out, TSPdata *tspdata, int *tour );
void recompute_obj( Param *param, TSPdata *tspdata, Vdata *vdata );

int compute_distance( double x1, double y1, double x2, double y2 );
int compute_cost( TSPdata *tspdata, int *tour );
int is_feasible( TSPdata *tspdata, int *tour );

/***** open the file with given mode *****************************************/
FILE *open_file( char *fname, char *mode ){
  FILE *fp;
  fp=fopen(fname,mode);
  if(fp==NULL){
    fprintf(stderr,"file not found: %s\n",fname);
    exit(EXIT_FAILURE);
  }
  return fp;
}

/***** malloc with error check ***********************************************/
void *malloc_e( size_t size ){
  void *s;
  if ( (s=malloc(size)) == NULL ) {
    fprintf( stderr, "malloc : not enough memory.\n" );
    exit(EXIT_FAILURE);
  }
  return s;
}


/***** copy and read the parameters ******************************************/
/***** Feel free to modify this subroutine. **********************************/
void copy_parameters( int argc, char *argv[], Param *param ){
    
  /**** copy the default parameters ****/
  param->timelim    = TIMELIM;
  param->givesol    = GIVESOL;
  param->outformat  = OUTFORMAT;
  strcpy(param->tourfile,TOURFILE);
  
  /**** read the parameters ****/
  if(argc>0 && (argc % 2)==0){
    printf("USAGE: ./%s [param_name, param_value] [name, value]...\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  else{
    int i;
    for(i=1; i<argc; i+=2){
      if(strcmp(argv[i],"timelim")==0)    param->timelim    = atoi(argv[i+1]);
      if(strcmp(argv[i],"givesol")==0)    param->givesol    = atoi(argv[i+1]);
      if(strcmp(argv[i],"outformat")==0)  param->outformat  = atoi(argv[i+1]);
      if(strcmp(argv[i],"tourfile")==0)   strcpy(param->tourfile,argv[i+1]);
    }
  }
}


/***** prepare memory space **************************************************/
/***** Feel free to modify this subroutine. **********************************/
void prepare_memory( TSPdata *tspdata, Vdata *vdata ){
  int k,n;
  n=tspdata->n;
  tspdata->x       = (double*)malloc_e(n*sizeof(double));
  tspdata->y       = (double*)malloc_e(n*sizeof(double));
  vdata->bestsol   = (int*)malloc_e(n*sizeof(int));
  /* the next line is just to give an initial solution */
  for(k=0;k<n;k++)
    vdata->bestsol[k]=k;
}

/***** reading the header of a file in TSPLIB format *************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_header( FILE *in, TSPdata *tspdata ){
  char str[MAX_STR],name[MAX_STR],dim[MAX_STR],type[MAX_STR],edge[MAX_STR],min[MAX_STR];
  int flag=0;

  for(;;){
    char *w,*u;
    /* error */
    if(fgets(str,MAX_STR,in)==NULL){
      fprintf(stderr,"error: invalid data input.\n");
      exit(EXIT_FAILURE);
    }
    /* halt condition */
    if(strcmp(str,"NODE_COORD_SECTION\n")==0){ break; }
    if(strcmp(str,"TOUR_SECTION\n")==0){ flag=1; break; }
    /* data input */
    w = strtok(str," :\n");
    u = strtok(NULL," :\n");
    if(w==NULL || u==NULL) continue;
    if(strcmp("NAME",w)==0)                  strcpy(name,u);
    if(strcmp("DIMENSION",w)==0)             strcpy(dim,u);
    if(strcmp("TYPE",w)==0)                  strcpy(type,u);
    if(strcmp("EDGE_WEIGHT_TYPE",w)==0)      strcpy(edge,u);
    if(strcmp("MIN_NODE_NUM",w)==0)          strcpy(min,u);
  }

  /* read a TSP instance */
  if(flag==0){
    strcpy(tspdata->name,name);
    tspdata->min_node_num=atoi(min);
    tspdata->n=atoi(dim);
    if(strcmp("TSP",type)!=0 || strcmp("EUC_2D",edge)!=0 ){
      fprintf(stderr,"error: invalid instance.\n");
      exit(EXIT_FAILURE);
    }
  }
  /* read a tour */
  else{
    if(strcmp("TOUR",type)!=0){
      fprintf(stderr,"error: invalid tour.\n");
      exit(EXIT_FAILURE);
    }
  }
}

/***** reading the file of TSP instance **************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_tspfile( FILE *in, TSPdata *tspdata, Vdata *vdata ){
  char str[MAX_STR];
  int k;

  /* reading the instance */
  read_header(in,tspdata);
  prepare_memory(tspdata,vdata);
  for(k=0;k<tspdata->n;k++){
    int dummy;
    if(fgets(str,MAX_STR,in)==NULL) break;
    if(strcmp(str,"EOF\n")==0) break;
    sscanf(str,"%d%lf%lf",&dummy, &(tspdata->x[k]),&(tspdata->y[k]));
  }
  if(k!=tspdata->n){
    fprintf(stderr,"error: invalid instance.\n");
    exit(EXIT_FAILURE);
  }
}

/***** read the tour in the TSPLIB format with feasibility check *************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_tourfile( FILE *in, TSPdata *tspdata, int *tour ){
  int k;

  read_header(in,tspdata);
  for(k=0;k<tspdata->n;k++){
    int val;
    if(fscanf(in,"%d",&val)==EOF) break;
    if(val==-1) break;
    tour[k]=val-1;
  }
  if(k!=tspdata->n){
    fprintf(stderr,"error: invalid tour.\n");
    exit(EXIT_FAILURE);
  }
}


/***** output the tour in the TSPLIB format **********************************/
/***** note: the output tour starts from the node "1" ************************/
void output_tour( FILE *out, TSPdata *tspdata, int *tour ){
  int k,idx=0;

  fprintf(out,"NAME : %s\n",tspdata->name);
  fprintf(out,"COMMENT : tour_length=%d\n",compute_cost(tspdata,tour));
  fprintf(out,"TYPE : TOUR\n");
  fprintf(out,"DIMENSION : %d\n",tspdata->n);
  fprintf(out,"TOUR_SECTION\n");
  for(k=0;k<tspdata->n;k++)
    if(tour[k]==0)
      { idx=k; break; }
  for(k=idx;k<tspdata->n;k++){
    if(tour[k]<0) break;
    fprintf(out,"%d\n",tour[k]+1);
  }
  for(k=0;k<idx;k++)
    fprintf(out,"%d\n",tour[k]+1);
  fprintf(out,"-1\n");
  fprintf(out,"EOF\n");
}

/***** output the tour in the TSP_VIEW format **********************************/
/***** note: the indices of the tour starts from "0" ***************************/
void output_tour_for_tsp_view( FILE *out, TSPdata *tspdata, int *tour ){
  int k;

  fprintf(out, "%d\n", tspdata->n);
  for(k=0; k<tspdata->n; k++){
    fprintf(out, "%g %g\n", tspdata->x[k], tspdata->y[k]);
  }
  for(k=0; k<tspdata->n; k++){
    if(tour[k]<0) break;
    fprintf(out,"%d\n", tour[k]);
  }
}

/***** check the feasibility and recompute the cost **************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void recompute_obj( Param *param, TSPdata *tspdata, Vdata *vdata ){
  if(!is_feasible(tspdata,vdata->bestsol)){
    fprintf(stderr,"error: the computed tour is not feasible.\n");
    exit(EXIT_FAILURE);
  }
  printf("recomputed tour length = %d\n",
	 compute_cost(tspdata,vdata->bestsol));
  printf("time for the search:   %7.2f seconds\n",
         vdata->endtime - vdata->starttime);
  printf("time to read the instance: %7.2f seconds\n",
         vdata->starttime - vdata->timebrid);
}

/***** cost of the tour ******************************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
int compute_cost( TSPdata *tspdata, int *tour ){
  int k,cost=0,n;
  n=tspdata->n;

  for(k=0;k<n-1;k++){
    if(tour[k+1]<0) break;
    cost += dist(tour[k],tour[k+1]);
  }
  cost += dist(tour[k],tour[0]);
  return cost;
}

/***** check the feasibility of the tour *************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
int is_feasible( TSPdata *tspdata, int *tour ){
  int k,n,*visited,flag=1,num_visited=0;
  n=tspdata->n;
  visited=(int*)malloc_e(n*sizeof(int));

  for(k=0;k<n;k++)
    visited[k]=0;
  for(k=0;k<n;k++){
    if(tour[k]<0)
      { break; }
    if(tour[k]>=n)
      { flag=0; break; }
    if(visited[tour[k]])
      { flag=0; break; }
    else{
      visited[tour[k]]=1;
      num_visited++;
    }
  }
  if(num_visited < tspdata->min_node_num) flag=0;

  free(visited);
  /* if tour is feasible (resp., not feasible), then flag=1 (resp., 0) */
  return flag;
}

/* my function and algorithm ********************************************************/
void print_array(int *a, int n){
    int i, s[n];

    
    for (i = 0; i < n; i++)
    {
        s[i] = a[i];
        printf("%d ", s[i]);
    }
    printf("\n");
}

void print_matrix(int n, int a[][n]){
    int i, j;
    
    for (i = 0; i < POPULATION; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%d ", a[i][j]);
        }
        printf("\n");
    }
}

void create_matrix(int n, int a[][n])
{
    int i, j;

    for (i = 0; i < POPULATION; i++)
    {
        for (j = 0; j < n; j++)
        {
            a[i][j] = rand() % (n - j) + 1;
        }
    }
}

void copy_array(int *a, int n, int j)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (i == j - 1)
        {
            a[i] = a[i + 1];
            for (i = j; i < n - 1; i++)
            {
                a[i] = a[i + 1];
            }
            a[n - 1] = 0;
        }
    }  
}

int min(int *a)
{
  int i, min;
  min = a[POPULATION-1];

  for (i = 0; i < POPULATION; i++)
  {
    if (min > a[i])
    {
    min = a[i];
    }
  }
  
  return min;
}

void order_representation(int n, int a[][n], int route[][n])
{
    int i, j, k;
    int order_list[n];

    for (i = 0; i < POPULATION; i++)
    {
        for (j = 0; j < n; j++)
        {
            order_list[j] = j + 1;
        }

        for (j = 0; j < n; j++)
        {
            k = a[i][j];
            route[i][j] = order_list[k - 1];
            copy_array(order_list, n, k);
        }
        
    }
    
}


void evaluate_route(int n, int route[][n], int *a, TSPdata *tspdata)
{
  int i, j, b[n];

  for (i = 0; i < POPULATION; i++)
  {
      a[i] = 0;
  }
   
  for (i = 0; i < POPULATION; i++)
  {
    for (j = 0; j < n; j++)
    {
      b[j] = route[i][j] - 1;
    }
    a[i] = compute_cost(tspdata, b);
  }

}

void ranking_selection(int *a, int n, int b[][n], int c[][n]) 
{
    int i, j, m, idx, d[POPULATION];

    for (i = 0; i < POPULATION; i++)
    {
        d[i] = a[i];
    }
    

    for (i = 0; i < POPULATION; i++)
    {
        m = min(d);
        for (j = 0; j < POPULATION; j++)
        {
            if (d[j] == m)
            {
                idx = j;
            }
        }

        for (j = 0; j < n; j++)
        {
            c[i][j] = b[idx][j];
        }
        d[idx] = 1000000;
    }

}


void tow_point_crossover(int n, int a[][n], int b[][n])   
{
  int i, j, point;
    
  if ((n % 2) == 0)
  {
    point = (n / 2) - 1;
  }else
  {
    point = ((n + 1)/2) - 1;
  }  
    
  for (i = 0; i < POPULATION; i=i+2)
  {
    for (j = 0; j < (point/2) + 1; j++)
    {
      a[i][j] = b[i][j];
      a[i + 1][j] = b[i + 1][j];
    }
    
    for (j = (point/2) + 1; j < (point + (point/2)) + 1; j++)
    {
      a[i][j] = b[i + 1][j];
      a[i + 1][j] = b[i][j];
    }

    for (j = (point + (point/2)) + 1; j < n; j++)
    {
      a[i][j] = b[i][j];
      a[i + 1][j] = b[i + 1][j];
    }
  }
    
}


void mutation(int n, int a[][n])
{
  int i, point, r = rand() % (POPULATION - 1) + 1;

  if ((n % 2) == 0)
  {
    point = (n / 2);
  }else
  {
    point = ((n + 1)/2);
  }  
  
  if (r > (POPULATION/2))
  {
    for (i = 0; i < point; i++)
  {
    a[r][i] = 1;
  }
  a[r][point] = 2;

  }else
  {
    for (i = point; i < n-10; i++)
    {
      a[r][i] = 1;
    }
    a[r][point+10] = 2;
  }
  
}


void rand_crossover(int n, int a[][n])
{
  int i, p;
  p = 3;

  for (i = 0; i < n; i++)
  {
    a[POPULATION-p][i] = a[p][i];
  }
  
}


void genetic_algorithm( Param *param, TSPdata *tspdata, Vdata *vdata )
{
  srand((unsigned int)time(NULL));

  int i, len, r1, r2;

  len = tspdata->n;

  int gene[POPULATION][len], route[POPULATION][len], fitness[POPULATION], gene_tmp[POPULATION][len];

  create_matrix(len, gene);

  for (i = 0; i < len; i++)
  {
    gene[0][i] = i;
  }

  for(i=0; i<tspdata->n; i++){
    vdata->bestsol[i] = i;
    }
  while(cpu_time() - vdata->starttime < param->timelim){
  order_representation(len, gene, route);
  evaluate_route(len, route, fitness, tspdata);
  ranking_selection(fitness, len, gene, gene_tmp);

  r1 = rand() % (20);
  r2 = rand() % (20);

  if (r1 == 7)
  {
    rand_crossover(len, gene_tmp);
  }


  tow_point_crossover(len, gene, gene_tmp);

  if (r1 ==  r2)
  {
    mutation(len, gene);
  }
  
  }

  order_representation(len, gene_tmp, route);
  for (i = 0; i < len; i++)
  {
    vdata->bestsol[i] = route[0][i] - 1;
  }

}



/***** main ******************************************************************/
int main(int argc, char *argv[]){

  Param     param;     /* parameters */
  TSPdata   tspdata;   /* data of TSP instance */
  Vdata     vdata;     /* various data often needed during search */

  vdata.timebrid = cpu_time();
  copy_parameters(argc, argv, &param);
  read_tspfile(stdin,&tspdata,&vdata);
  if(param.givesol==1) read_tourfile(stdin,&tspdata,vdata.bestsol);
  vdata.starttime = cpu_time();

  /*****

    Write your algorithm here. 
    Of course you can add your subroutines outside main(). 
    At this point, the instance data is stored in the structure "tspdata". 

       tspdata.name :     the name of the instance
       tspdata.n :        the number of nodes n
       tspdata.x[k] :     x-coordinate of node k (k = 0,1,...,n-1)
       tspdata.y[k] :     y-coordinate of node k (k = 0,1,...,n-1)
       tspdata.min_node_num :
                          the number of nodes that the solution must contain

   You can use the following macro to get the distance between node k and l.
       dist(k,l)    :  the distance between node k and l
                       (k,l = 0,1,...,n-1), where
                       dist(k,k) = 0 and
                       dist(k,l) = dist(l,k)
                       for any k and l

    Store your best tour (solution) in vdata.bestsol. If number of nodes in the tour
    is less than tspdata.n, you have to fill from bestsol[m] to bestsol[n-1] with -1.
    The "compute_cost()" will compute its objective value. The "is_feasible()" will
    return whether the solution is feasible or not. The format of vdata.bestsol is:

       vdata.bestsol[i] = k    if the i-th visited node of the tour is node k,
       vdata.bestsol[i] = -1   if i is greater than the number of visited nodes of the tour.
                               where i,k = 0,1,...,n-1.

  *****/

  genetic_algorithm(&param,&tspdata,&vdata);

  vdata.endtime = cpu_time();
  recompute_obj(&param,&tspdata,&vdata);
  if(param.outformat==1){
    output_tour(open_file(param.tourfile,"w"),&tspdata,vdata.bestsol);}
  else if(param.outformat==2){
    output_tour_for_tsp_view(open_file(param.tourfile,"w"),&tspdata,vdata.bestsol);}

  return EXIT_SUCCESS;
}
