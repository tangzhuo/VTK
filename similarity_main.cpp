#define MAXDocNum 1024
# include <stdio.h>
#include <math.h>
# include "common.h"
using namespace std;

char docfile[MAXDocNum];           /* file with training examples */
char outfile[MAXDocNum];         /* file for resulting kernel compute */
//char modelfile[200];         /* file for resulting classifier */

void   read_input_parameters(int, char **, char *,  char *,long *, long *, 
			      KERNEL_PARM *);
void   wait_any_key();
void   print_help();
void compute_p_s(DOC *docs);

int main (int argc, char* argv[])
{ 
  DOC *docs;  /* training examples */
  long max_docs,max_words_doc;
  long totwords,totdoc,ll,i;
  long kernel_cache_size;
  long runtime_start, runtime_end;
  double *target;
  
  int wordNum;
  int size;
  
  KERNEL_PARM kernel_parm;
//  MODEL model;
 
  FILE *fp;
  read_input_parameters(argc,argv,docfile,outfile,&verbosity,
			&kernel_cache_size,&kernel_parm);
  /*
  fopen("","rb");
  if(f==NULL){
	  printf("Embedding file not found\n");
  }
  else{
	fscanf(f, "%d", &wordNum);
	fscanf(f, "%d", &size);
	for (int b = 0; b < wordNum; b++) {
		char ch;
		fscanf(f, "%s%c", str, &ch);
  }
  */
  /* STANDARD SVM KERNELS */
      
  LAMBDA = kernel_parm.lambda; // to make faster the kernel evaluation

  if(verbosity>=1) {
    printf("Scanning examples..."); 
	fflush(stdout);
  }
  nol_ll(docfile,&max_docs,&max_words_doc,&ll); /* scan size of input file */
 // printf("done\n"); 
  max_words_doc+=10;
  ll+=10;
  max_docs+=2;
  if(verbosity>=1) {
    printf("done\n");
	fflush(stdout);
  }

  docs = (DOC *)my_malloc(sizeof(DOC)*max_docs);         /* feature vectors */
 // target = (double *)my_malloc(sizeof(double)*max_docs); /* target values */

  /* start timer, exclude I/O */
	runtime_start = get_runtime();

 //printf("\nMax docs: %ld, approximated number of feature occurences %ld, maximal length of a line %ld\n\n",max_docs,max_words_doc,ll);
  read_documents(docfile,docs,max_words_doc,ll,&totwords,&totdoc,&kernel_parm);
  
  printf("\nNumber of examples: %ld, linear space size: %ld\n\n",totdoc,totwords);
 
 //if(kernel_parm.kernel_type==5) totwords=totdoc; // The number of features is proportional to the number of parse-trees, i.e. totdoc 
  				                 // or should we still use totwords to approximate svm_maxqpsize for the Tree Kernel (see hideo.c) ???????

  

  /* Warning: The model contains references to the original data 'docs'.
     If you want to free the original data, and only keep the model, you 
     have to make a deep copy of 'model'. */
 // write_model(modelfile,&model);

 // free(model.supvec);
  //free(model.alpha);
  //free(model.index);
  /* end timer */
   runtime_end = get_runtime();
	
   printf("Training [%s] time in cpu seconds (excluding I/O): %.2f\n", docfile,((float) runtime_end - (float) runtime_start)/100.0);  
	
   printf("===========================\n");

   //write the kernel score to file
  // compute_p_s(docs);
   fp = fopen(outfile,"w");
   for(i=0;i<totdoc;i++){
	 //  printf("写入 %s %3f\n",outfile,docs[i].k_score);
	   fprintf(fp,"%.3f\n",docs[i].k_score*((double)5));
	   
	  // printf("\n\k_score = %.4f\n\n",docs[i].k_score);
   }
   fclose(fp);
  
   for(i=0;i<totdoc;i++){
     freeExample(&docs[i]);
   }
  
  free(docs);
 // free(target);

  system("pause");
  return(0);
}


/*---------------------------------------------------------------------------*/

void read_input_parameters(int argc,char *argv[],char *docfile,char *outfile,
			   long *verbosity,long *kernel_cache_size,
			   KERNEL_PARM *kernel_parm)
{
  long i;
  char type[100];
  
  /* set default */
  strcpy (outfile, "k_output");
  (*verbosity)=1;
  (*kernel_cache_size)=40;
  kernel_parm->kernel_type=0;
  kernel_parm->poly_degree=3;
  kernel_parm->rbf_gamma=1.0;
  kernel_parm->coef_lin=1;
  kernel_parm->coef_const=1;
  kernel_parm->lambda=0.2;
  kernel_parm->tree_constant=1;
  kernel_parm->second_kernel=1;
  kernel_parm->first_kernel=1; 
  kernel_parm->normalization=3;
  kernel_parm->combination_type='T'; //no combination
  kernel_parm->vectorial_approach_standard_kernel='S';
  kernel_parm->vectorial_approach_tree_kernel='S';
  kernel_parm->mu=0.1; // Default Duffy and Collins Kernel 
  kernel_parm->tree_kernel_param=0; // Default no params
  strcpy(kernel_parm->custom,"empty");
  strcpy(type,"c");
  
  for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) 
      { 
      case '?': print_help(); exit(0);
      case 'z': i++; strcpy(type,argv[i]); break;
      case 'v': i++; (*verbosity)=atol(argv[i]); break;
      case 'm': i++; (*kernel_cache_size)=atol(argv[i]); break;
      case 't': i++; kernel_parm->kernel_type=atol(argv[i]); break;
      case 'd': i++; kernel_parm->poly_degree=atol(argv[i]); break;
      case 'g': i++; kernel_parm->rbf_gamma=atof(argv[i]); break;
      case 's': i++; kernel_parm->coef_lin=atof(argv[i]); break;
      case 'r': i++; kernel_parm->coef_const=atof(argv[i]); break;
      case 'u': i++; strcpy(kernel_parm->custom,argv[i]); break;
      case 'L': i++; kernel_parm->lambda=atof(argv[i]); break;
      case 'T': i++; kernel_parm->tree_constant=atof(argv[i]); break;
      case 'C': i++; kernel_parm->combination_type=*argv[i]; break;
      case 'F': i++; kernel_parm->first_kernel=atoi(argv[i]); break;
      case 'S': i++; kernel_parm->second_kernel=atoi(argv[i]);  break;
      case 'V': i++; kernel_parm->vectorial_approach_standard_kernel=*argv[i]; break;
      case 'W': i++; kernel_parm->vectorial_approach_tree_kernel=*argv[i]; break;
      case 'M': i++; kernel_parm->mu=atof(argv[i]); break; 
	  case 'N': i++; kernel_parm->normalization=atoi(argv[i]); break; 
	  case 'U': i++; kernel_parm->tree_kernel_param=atoi(argv[i]); break; // user defined parameters 


      default: printf("\nUnrecognized option %s!\n\n",argv[i]);
	       print_help();
	       exit(0);
      }
  }

  LAMBDA = kernel_parm->lambda; // to make faster the kernel evaluation 
  LAMBDA2 = LAMBDA*LAMBDA;
  MU= kernel_parm->mu;
  TKGENERALITY=kernel_parm->first_kernel;
  PARAM_VECT=kernel_parm->tree_kernel_param;
  if(PARAM_VECT == 1) 
	  read_input_tree_kernel_param(); // if there is the file tree_kernel.param load paramters
	

  if(i>=argc) {
    printf("\nNot enough input parameters!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  
  char s1[5];
  strcpy (docfile, argv[i]);
  if((i+1)<argc) {
    strcpy (outfile, argv[i+1]);
  }
//  strcpy(outfile,"TKW2V_sim_300_30_");
/*  
  sprintf(outfile,"%d",r1);
  sprintf(s1,"%d",r2);
  strcat(outfile,"_");
  strcat(outfile,s1);
  strcat(outfile,"_");
  strcat(outfile,"TK_");
  strcat(outfile,docfile);
  */
  printf("进来了！！ MU=%lf ,LAMDA=%lf",MU,LAMBDA);
}

void wait_any_key()
{
  printf("\n(more)\n");
  (void)getc(stdin);
}

void print_help()
{
  printf("\nTree Kernels in SVM-light %s : SVM Learning module %s\n",VERSION,VERSION_DATE);
  printf("by Alessandro Moschitti, moschitti@info.uniroma2.it\n");
  printf("University of Rome \"Tor Vergata\"\n\n");

  //copyright_notice();
  printf("   usage: svm_learn [options] example_file model_file\n\n");
  printf("Arguments:\n");
  printf("         example_file-> file with training data\n");
  printf("         model_file  -> file to store learned decision rule in\n");

  printf("General options:\n");
  printf("         -?          -> this help\n");
  printf("         -v [0..3]   -> verbosity level (default 1)\n");
  printf("Learning options:\n");
  printf("         -z {c,r,p}  -> select between classification (c), regression (r),\n");
  printf("                        and preference ranking (p) (default classification)\n");
  printf("         -c float    -> C: trade-off between training error\n");
  printf("                        and margin (default [avg. x*x]^-1)\n");
  printf("         -w [0..]    -> epsilon width of tube for regression\n");
  printf("                        (default 0.1)\n");
  printf("         -j float    -> Cost: cost-factor, by which training errors on\n");
  printf("                        positive examples outweight errors on negative\n");
  printf("                        examples (default 1) (see [4])\n");
  printf("         -b [0,1]    -> use biased hyperplane (i.e. x*w+b>0) instead\n");
  printf("                        of unbiased hyperplane (i.e. x*w>0) (default 1)\n");
  printf("         -i [0,1]    -> remove inconsistent training examples\n");
  printf("                        and retrain (default 0)\n");
  printf("Performance estimation options:\n");
  printf("         -x [0,1]    -> compute leave-one-out estimates (default 0)\n");
  printf("                        (see [5])\n");
  printf("         -o ]0..2]   -> value of rho for XiAlpha-estimator and for pruning\n");
  printf("                        leave-one-out computation (default 1.0) (see [2])\n");
  printf("         -k [0..100] -> search depth for extended XiAlpha-estimator \n");
  printf("                        (default 0)\n");
  printf("Transduction options (see [3]):\n");
  printf("         -p [0..1]   -> fraction of unlabeled examples to be classified\n");
  printf("                        into the positive class (default is the ratio of\n");
  printf("                        positive and negative examples in the training data)\n");

  printf("Kernel options:\n");
  printf("         -t int      -> type of kernel function:\n");
  printf("                        0: linear (default)\n");
  printf("                        1: polynomial (s a*b+c)^d\n");
  printf("                        2: radial basis function exp(-gamma ||a-b||^2)\n");
  printf("                        3: sigmoid tanh(s a*b + c)\n");
  printf("                        4: user defined kernel from kernel.h\n");

  printf("                        5: combination of forest and vector sets according to W, V, S, C options\n");
  printf("                        11: re-ranking based on trees (each instance must have two trees),\n");
  printf("                        12: re-ranking based on vectors (each instance must have two vectors)\n");
  printf("                        13: re-ranking based on both tree and vectors (each instance must have\n");
  printf("                            two trees and two vectors)  \n");
  printf("         -W [S,A]    -> with an 'S', a tree kernel is applied to the sequence of trees of two input\n");
  printf("                        forests and the results are summed;  \n");
  printf("                     -> with an 'A', a tree kernel is applied to all tree pairs from the two forests\n");
  printf("                        (default 'S')\n");
  printf("         -V [S,A]    -> same as before but regarding sequences of vectors are used (default 'S' and\n");
  printf("                        the type of vector-based kernel is specified by the option -S)\n");
  printf("         -S [0,4]    -> kernel to be used with vectors (default polynomial of degree 3,\n");
  printf("                        i.e. -S = 1 and -d = 3)\n");
  printf("         -C [*,+,T,V]-> combination operator between forests and vectors (default 'T')\n");
  printf("                     -> 'T' only the contribution from trees is used (specified by option -W)\n");
  printf("                     -> 'V' only the contribution from vectors is used (specified by option -V)\n");
  printf("                     -> '+' or '*' sum or multiplication of the contributions from vectors and \n");
  printf("                            trees (default T) \n");
  printf("         -D [0,1]    -> 0, SubTree kernel or 1, SubSet Tree kernels (default 1)\n");
  printf("         -L float    -> decay factor in tree kernel (default 0.4)\n");
  printf("         -S [0,4]    -> kernel to be used with vectors (default polynomial of degree 3, \n");
  printf("                        i.e. -S = 1 and -d = 3)\n");
  printf("         -T float    -> multiplicative constant for the contribution of tree kernels when -C = '+'\n");
  printf("         -N float    -> 0 = no normalization, 1 = tree normalization, 2 = vector normalization and \n");
  printf("                        3 = tree normalization of both trees and vectors. The normalization is applied \n");
  printf("                        to each individual tree or vector (default 3).\n");

  printf("         -u string   -> parameter of user defined kernel\n");
  printf("         -d int      -> parameter d in polynomial kernel\n");
  printf("         -g float    -> parameter gamma in rbf kernel\n");
  printf("         -s float    -> parameter s in sigmoid/poly kernel\n");
  printf("         -r float    -> parameter c in sigmoid/poly kernel\n");
  printf("         -u string   -> parameter of user defined kernel\n");
 
  printf("Optimization options (see [1]):\n");
  printf("         -q [2..]    -> maximum size of QP-subproblems (default 10)\n");
  printf("         -n [2..q]   -> number of new variables entering the working set\n");
  printf("                        in each iteration (default n = q). Set n<q to prevent\n");
  printf("                        zig-zagging.\n");
  printf("         -m [5..]    -> size of cache for kernel evaluations in MB (default 40)\n");
  printf("                        The larger the faster...\n");
  printf("         -e float    -> eps: Allow that error for termination criterion\n");
  printf("                        [y [w*x+b] - 1] >= eps (default 0.001)\n");
  printf("         -h [5..]    -> number of iterations a variable needs to be\n"); 
  printf("                        optimal before considered for shrinking (default 100)\n");
  printf("         -f [0,1]    -> do final optimality check for variables removed\n");
  printf("                        by shrinking. Although this test is usually \n");
  printf("                        positive, there is no guarantee that the optimum\n");
  printf("                        was found if the test is omitted. (default 1)\n");
  printf("Output options:\n");
  printf("         -l string   -> file to write predicted labels of unlabeled\n");
  printf("                        examples into after transductive learning\n");
  printf("         -a string   -> write all alphas to this file after learning\n");
  printf("                        (in the same order as in the training set)\n");
  wait_any_key();
  printf("\nMore details in:\n");
  printf("[1] T. Joachims, Making Large-Scale SVM Learning Practical. Advances in\n");
  printf("    Kernel Methods - Support Vector Learning, B. Sch鰈kopf and C. Burges and\n");
  printf("    A. Smola (ed.), MIT Press, 1999.\n");
  printf("[2] T. Joachims, Estimating the Generalization performance of an SVM\n");
  printf("    Efficiently. International Conference on Machine Learning (ICML), 2000.\n");
  printf("[3] T. Joachims, Transductive Inference for Text Classification using Support\n");
  printf("    Vector Machines. International Conference on Machine Learning (ICML),\n");
  printf("    1999.\n");
  printf("[4] K. Morik, P. Brockhausen, and T. Joachims, Combining statistical learning\n");
  printf("    with a knowledge-based approach - A case study in intensive care  \n");
  printf("    monitoring. International Conference on Machine Learning (ICML), 1999.\n");
  printf("[5] T. Joachims, Learning to Classify Text Using Support Vector\n");
  printf("    Machines: Methods, Theory, and Algorithms. Dissertation, Kluwer,\n");
  printf("    2002.\n\n");
  printf("\nFor Tree-Kernel details:\n");
  printf("[6] A. Moschitti, A study on Convolution Kernels for Shallow Semantic Parsing.\n");
  printf("    In proceedings of the 42-th Conference on Association for Computational\n");
  printf("    Linguistic, (ACL-2004), Barcelona, Spain, 2004.\n\n");
  printf("[7] A. Moschitti, Making tree kernels practical for natural language learning.\n");
  printf("    In Proceedings of the Eleventh International Conference for Computational\n");
  printf("    Linguistics, (EACL-2006), Trento, Italy, 2006.\n\n");
  
}
