#define _CRT_SECURE_NO_WARNINGS
# include "common.h"
#include <string>
#include <stdlib.h>
#include "tchar.h"
#include <math.h>
#include <map>
#include <vector>
#include <ctype.h>
#include <algorithm>
#include <iostream>
using namespace std;

#define DPS(i,j) (*(DPS+(i)*(m+1)+j))
#define DP(i,j) (*(DP+(i)*(m+1)+j))

map<string,double*> dict;
double LAMBDA2;
double LAMBDA;
double SIGMA;
double MU;
double REMOVE_LEAVES = 0; // if equal to MU*LAMBDA2, it removes the leaves contribution;
const double eps = 1e-8;

short  TKGENERALITY;  //store the generality of the kernel PT is 2 whereas SST and ST is 1
                      // used to load the opportune data structures (with SST and ST uses faster approach)
short  PARAM_VECT; // if 1 the vector of parameters is defined

double delta_matrix[MAX_NUMBER_OF_NODES][MAX_NUMBER_OF_NODES];


// local functions
void determine_sub_lists(FOREST *a, FOREST *b, nodePair *intersect, int *n);
double Delta_SST_ST( TreeNode * Nx, TreeNode * Nz); // delta function SST and ST kernels
double Delta_SK(TreeNode **Sx, TreeNode ** Sz, int n, int m); // delta over children
double Delta_PT( TreeNode * Nx, TreeNode * Nz,int flag); // Delta for PT
//double Delta_SK(TreeNode **Sx, TreeNode ** Sz, char * sxpname,char * szpname,int n, int m); // delta over children
//double Delta_PT( TreeNode * Nx, TreeNode * Nz,char * NxPsname,char * NzPsname,int flag); // Delta for PT
double Delta( TreeNode * Nx, TreeNode * Nz);


void evaluateNorma(KERNEL_PARM * kernel_parm, DOC * d); // evaluate norm of trees and vectors

double basic_kernel(KERNEL_PARM *kernel_parm, DOC *a, DOC *b, int i, int j);  // svm-light kernels
//double tree_kernel(KERNEL_PARM *kernel_parm, DOC * a, DOC * b, int i, int j); // ST and SST kernels
double SK(TreeNode **Sx, TreeNode ** Sz, int n, int m); // string kernel
double string_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b, int i, int j);

double choose_second_kernel(KERNEL_PARM *kernel_parm, DOC *a, DOC *b); //choose a standard kernel
double choose_tree_kernel(KERNEL_PARM *kernel_parm, DOC *a, DOC *b); //choose a tree kernel


// kernel combinations
double sequence(KERNEL_PARM * kernel_parm, DOC * a, DOC * b);
double AVA(KERNEL_PARM * kernel_parm, DOC * a, DOC * b);
double AVA_tree_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b);
double sequence_tree_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b);
double evaluateParseTreeKernel(nodePair *pairs, int n);
double advanced_kernels(KERNEL_PARM * kernel_parm, DOC * a, DOC * b);
double sequence_ranking(KERNEL_PARM * kernel_parm, DOC * a, DOC * b, int memberA, int memberB);//all_vs_all vectorial kernel
double SRL2008(KERNEL_PARM *kernel_parm, DOC *a, DOC *b);
double ACL2008_Entailment_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b);
double ACL2007_Entailment_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b);
double ACL2008(KERNEL_PARM *kernel_parm, DOC *a, DOC *b);

// Tree Kernels
double hasdigit(char * Nx_sname,char *Nz_sname);
int isantonym(char *word1, char *word2);
int ReadEmbedding(const char *file_name);
double getW2V_similarity(char * Nx_sname,char *Nz_sname);

double Delta( TreeNode * Nx, TreeNode * Nz){
    int i;
    double prod=1;
   
//printf("Delta Matrix: %1.30lf node1:%s node2:%s, LAMBDA %lf, SIGMA %lf\n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nz->sName,LAMBDA,SIGMA);

	if(delta_matrix[Nx->nodeID][Nz->nodeID]>=0) {
//printf("Delta Matrix: %1.30lf ?diverso da -1 boh \n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nz->sName,LAMBDA,SIGMA);

	      return delta_matrix[Nx->nodeID][Nz->nodeID]; // Case 0 (Duffy and Collins 2002);
       }
	else 
	   if(Nx->pre_terminal || Nz->pre_terminal)
	      return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA); // case 1 
	   else{
		   for(i=0;i<Nx->iNoOfChildren;i++)
	   	       if(strcmp(Nx->pChild[i]->production,Nz->pChild[i]->production)==0)
		          prod*= (SIGMA+Delta( Nx->pChild[i], Nz->pChild[i])); // case 2

       	   return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA*prod);
	   }
}

double evaluateParseTreeKernel(nodePair *pairs, int n){

    int i;
    double sum=0,contr;
	   
	   for(i=0;i<n;i++){
		contr=Delta(pairs[i].Nx,pairs[i].Nz);
//		printf("Score for the pairs (%s , %s): %f\n",pairs[i].Nx->sName,pairs[i].Nz->sName,contr);fflush(stdout);
		sum+=contr;
	   }
// printf("\n\nFORM EVALUATE KERNEL = %f \n",sum); 

	   return sum;
}

//-------------------------------------------------------------------------------------------------------
// SST and ST KERNELS see Moschitti - ECML 2006

/*
double Delta_SST_ST( TreeNode * Nx, TreeNode * Nz){
	int i;
	double prod=1;
	
	//printf("Delta Matrix: %1.30lf node1:%s node2:%s, LAMBDA %lf, SIGMA %lf\n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nz->sName,LAMBDA,SIGMA);
	
	if(delta_matrix[Nx->nodeID][Nz->nodeID]>=0)
		return delta_matrix[Nx->nodeID][Nz->nodeID]; // Case 0 (Duffy and Collins 2002);
	else 
		if(Nx->pre_terminal || Nz->pre_terminal)
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA); // case 1 
		else{
			for(i=0;i<Nx->iNoOfChildren;i++)
				if(Nx->pChild[i]->production!=NULL && Nz->pChild[i]->production!= NULL)
					if(strcmp(Nx->pChild[i]->production, Nz->pChild[i]->production)==0)
						prod*= (SIGMA+Delta_SST_ST( Nx->pChild[i], Nz->pChild[i])); // case 2
					else prod*=SIGMA; 	          
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA*prod);
		}
}

double getW2V_similarity(char *Nx_sname,char *Nz_sname){
	if(strcmp(Nx_sname,Nz_sname)==0)
		return 1.0;
	else 
		return 0.0;
}*/


double Delta_SST_ST( TreeNode * Nx, TreeNode * Nz){
	int i;
	double prod=1;
	
	//printf("Delta Matrix: %1.30lf node1:%s node2:%s, LAMBDA %lf, SIGMA %lf\n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nz->sName,LAMBDA,SIGMA);
	
	if(delta_matrix[Nx->nodeID][Nz->nodeID]>=0)
		return delta_matrix[Nx->nodeID][Nz->nodeID]; // Case 0 (Duffy and Collins 2002);
	else 
		if((Nx->pre_terminal==-1) || (Nz->pre_terminal==-1)){
			printf("%s and %s 有一个为-1！\n",Nx->sName,Nz->sName);
		//	return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA); // case 1 
			return (delta_matrix[Nx->nodeID][Nz->nodeID] = getW2V_similarity(Nx->sName,Nz->sName));
		}
		//else if((Nx->pre_terminal==-1) || (Nz->pre_terminal==-1))
		//	return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA);
		else{
			//printf("%s and %s 进来递归了！\n",Nx->sName,Nz->sName);
			for(i=0;i<Nx->iNoOfChildren;i++)
				if(Nx->pChild[i]->production!=NULL && Nz->pChild[i]->production!= NULL)
					if(strcmp(Nx->pChild[i]->production, Nz->pChild[i]->production)==0)
						prod*= (SIGMA+Delta_SST_ST( Nx->pChild[i], Nz->pChild[i])); // case 2
					else prod*=SIGMA; 	          
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA*prod);
		}
}

double evaluate_SST_ST(nodePair *pairs, int n){
	
	int i;
	double sum=0,contr;
	
	for(i=0;i<n;i++){
		//printf("Score for the pairs (%s , %s)",pairs[i].Nx->sName,pairs[i].Nz->sName);fflush(stdout);
		//printf("\ntree 1: "); writeTreeString(pairs[i].Nx); printf("\ntree 2: "); writeTreeString(pairs[i].Nz); printf("\n");
		//fflush(stdout);
		// 
		contr=Delta_SST_ST(pairs[i].Nx,pairs[i].Nz);
		//printf("%f\n",contr);fflush(stdout);
		sum+=contr;
	}
	//printf("FINAL KERNEL = %f \n\n\n",sum); 
	
	return sum;
}


//-------------------------------------------------------------------------------------------------------
// Kernel with slot used in ACL07 and ECIR07 to simulate a fast PT (within the first tree level)
//-------------------------------------------------------------------------------------------------------

double Delta_ACL07( TreeNode * Nx, TreeNode * Nz){
  int i;
  double prod=1;
	
  //printf("Delta Matrix: %1.30lf node1:%s node2:%s, LAMBDA %lf, SIGMA %lf\n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nz->sName,LAMBDA,SIGMA);
	
  if(delta_matrix[Nx->nodeID][Nz->nodeID]>=0) {
		//printf("Delta Matrix: %1.30lf diverso da -1 boh \n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nz->sName,LAMBDA,SIGMA);
		
		return delta_matrix[Nx->nodeID][Nz->nodeID]; // Case 0 (Duffy and Collins 2002);
  }
  else 
		if(Nx->pre_terminal || Nz->pre_terminal)
			//Alessandro
			if (strcmp(Nx->pChild[0]->sName,"null")==0 ||  strcmp(Nz->pChild[0]->sName,"null")==0) return 0; //don't consider null slots
			else
				return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA); // case 1 
			else{
				for(i=0;i<Nx->iNoOfChildren;i++)
					if(strcmp(Nx->pChild[i]->production,Nz->pChild[i]->production)==0)
						prod*= (SIGMA+Delta_ACL07( Nx->pChild[i], Nz->pChild[i])); // case 2
				
				return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA*prod);
			}
}

double evaluate_SST_ACL07(nodePair *pairs, int n){
	
	int i;
	double sum=0,contr;
	
	for(i=0;i<n;i++){
		//printf("Score for the pairs (%s , %s)",pairs[i].Nx->sName,pairs[i].Nz->sName);fflush(stdout);
		//printf("\ntree 1: "); writeTreeString(pairs[i].Nx); printf("\ntree 2: "); writeTreeString(pairs[i].Nz); printf("\n");
		//fflush(stdout);
		
		// Don't consider BOX label in the Bag of X trees
		if(strcmp(pairs[i].Nx->sName,"BOX")!= 0 && strcmp(pairs[i].Nz->sName,"BOX")!= 0) // to be removed as it is verified that is useless
		{
			contr=Delta_ACL07(pairs[i].Nx,pairs[i].Nz);
			//printf("%f\n",contr);fflush(stdout);
			sum+=contr;
		}
	}
	// printf("FINAL KERNEL = %f \n\n",sum); 
	
	return sum;
}



//-------------------------------------------------------------------------------------------------------
// A more general Collins' based kernel (leaves are considered as features)
//-------------------------------------------------------------------------------------------------------



double Delta_GSST( TreeNode * Nx, TreeNode * Nz){
	int i;
	double prod=1;
	if(delta_matrix[Nx->nodeID][Nz->nodeID]>=0) return delta_matrix[Nx->nodeID][Nz->nodeID]; // cashed
	else{
		for(i=0;i<Nx->iNoOfChildren && i<Nz->iNoOfChildren ;i++)
			if(Nx->pChild[i]->production != NULL && Nz->pChild[i]->production != NULL &&
				 strcmp(Nx->pChild[i]->production,Nz->pChild[i]->production)==0 
				 && Nx->pChild[i]->pre_terminal != -1 && Nz->pChild[i]->pre_terminal != -1)
				prod*= (1+Delta_GSST( Nx->pChild[i], Nz->pChild[i])); // case 2
		return (delta_matrix[Nx->nodeID][Nz->nodeID]=LAMBDA*prod);
	}
}

double evaluate_GSST(nodePair *pairs, int n){
	
	int i;
	double sum=0,contr;
	
	for(i=0;i<n;i++){
		//printf("Score for the pairs (%s , %s)",pairs[i].Nx->sName,pairs[i].Nz->sName);fflush(stdout);
		//printf("\ntree 1: "); writeTreeString(pairs[i].Nx); printf("\ntree 2: "); writeTreeString(pairs[i].Nz); printf("\n");
		//fflush(stdout);
		//        if(pairs[i].Nx->iNoOfChildren && pairs[i].Nz->iNoOfChildren) 
		contr=Delta_GSST(pairs[i].Nx,pairs[i].Nz);
		//        else contr=0;
		//printf("%f\n",contr);fflush(stdout);
		
		sum+=contr;
	}
	//printf("FINAL KERNEL = %f \n\n\n",sum); 
	
	return sum;
}

//-------------------------------------------------------------------------------------------------------
// Partial Tree Kernel - see Moschitti - ECML 2006
//-------------------------------------------------------------------------------------------------------


#ifdef FAST

//double Delta_SK(TreeNode **Sx, TreeNode ** Sz, int n, int m,char * sxpname,char * szpname,int flag){
double Delta_SK(TreeNode **Sx, TreeNode ** Sz, int n, int m,int flag){
	double DPS[MAX_NUMBER_OF_CHILDREN_PT][MAX_NUMBER_OF_CHILDREN_PT];
	double DP[MAX_NUMBER_OF_CHILDREN_PT][MAX_NUMBER_OF_CHILDREN_PT];
	double kernel_mat[MAX_NUMBER_OF_CHILDREN_PT];
  
	int i,j,l,p;
	double K;
  
  
	p = n;
	if (m<n) 
	  p=m;
	if (p>MAX_CHILDREN) 
	  p=MAX_CHILDREN;
  
	//  if(n==0 || m==0 || m!=n) return 0;
  
	for (j=0; j<=m; j++)
		for (i=0; i<=n; i++)
			DPS[i][j]=DP[i][j]=0;
	
	kernel_mat[0]=0;
	for (i=1; i<=n; i++)
		for (j=1; j<=m; j++)
			
			if((strcmp((*(Sx+i-1))->sName,(*(Sz+j-1))->sName)==0) )
			//	|| (((*(Sx+i-1))->iNoOfChildren==0) && (*(Sz+j-1))->iNoOfChildren==0)) 
			{
				//printf("jinlaile. Nx(%s)= %d,Nz(%s)=%d\n",(*(Sx+i-1))->sName,(*(Sx+i-1))->nodeID,(*(Sz+j-1))->sName,(*(Sz+j-1))->nodeID);
				DPS[i][j]=Delta_PT(*(Sx+i-1),*(Sz+j-1),flag);
				//DPS[i][j]=Delta_PT(*(Sx+i-1),*(Sz+j-1),sxpname,szpname,flag);
				kernel_mat[0]+=DPS[i][j];
			}
			
			else
				DPS[i][j]=0;

	for(l=1;l<p;l++){
		kernel_mat[l]=0;
		for (j=0; j<=m; j++) 
			DP[l-1][j]=0;
		for (i=0; i<=n; i++) 
			DP[i][l-1]=0;
		
		for (i=l; i<=n; i++)
			for (j=l; j<=m; j++){
				DP[i][j] = DPS[i][j]+LAMBDA*DP[i-1][j]
						+ LAMBDA*DP[i][j-1]
						- LAMBDA2*DP[i-1][j-1];
				
				if( (strcmp((*(Sx+i-1))->sName,(*(Sz+j-1))->sName)==0))
				//	 || (((*(Sx+i-1))->iNoOfChildren==0) && (*(Sz+j-1))->iNoOfChildren==0)) 
				{
					 // printf("dierci\n");
				    	DPS[i][j] = Delta_PT(*(Sx+i-1),*(Sz+j-1),flag)* DP[i-1][j-1];
				//		DPS[i][j] = Delta_PT(*(Sx+i-1),*(Sz+j-1),sxpname,szpname,flag)* DP[i-1][j-1];
						kernel_mat[l] += DPS[i][j];
				}
				
						// else DPS[i][j] = 0;             
			}
				//      printf("\n----------------------------------\n"); printf("DPS i:%d, j:%d, l:%d\n",n,m,l+1);stampa_math(DPS,n,m);printf("DP\n");stampa_math(DP,n,m); 
	}
				
	//  K=kernel_mat[p-1];
	K=0;
	for(l=0;l<p;l++){
		K+=kernel_mat[l];
		//printf("String kernel of length %d: %1.7f \n\n",l+1,kernel_mat[l]);
    }
    return K;
	
	//   printf("\nDPS\n"); stampa_math(DPS,n,m); printf("DP\n");  stampa_math(DP,n,m); 
	
	//   printf("kernel: n=%d m=%d, %s %s \n\n",n,m,(*(Sx))->sName,(*(Sz))->sName);
	
	
}

#endif

/*
#ifdef FAST

double Delta_SK(TreeNode **Sx, TreeNode ** Sz, int n, int m,int flag){
	
	double DPS[MAX_NUMBER_OF_CHILDREN_PT][MAX_NUMBER_OF_CHILDREN_PT];
	double DP[MAX_NUMBER_OF_CHILDREN_PT][MAX_NUMBER_OF_CHILDREN_PT];
	double kernel_mat[MAX_NUMBER_OF_CHILDREN_PT];
  
	int i,j,l,p;
	double K;
  
  
	p = n;
	if (m<n) 
	  p=m;
	if (p>MAX_CHILDREN) 
	  p=MAX_CHILDREN;
  
	//  if(n==0 || m==0 || m!=n) return 0;
  
	for (j=0; j<=m; j++)
		for (i=0; i<=n; i++)
			DPS[i][j]=DP[i][j]=0;
	
	kernel_mat[0]=0;
	for (i=1; i<=n; i++)
		for (j=1; j<=m; j++)
			
			if((strcmp((*(Sx+i-1))->sName,(*(Sz+j-1))->sName)==0) ) 
			{
				//printf("jinlaile. Nx(%s)= %d,Nz(%s)=%d\n",(*(Sx+i-1))->sName,(*(Sx+i-1))->nodeID,(*(Sz+j-1))->sName,(*(Sz+j-1))->nodeID);
				DPS[i][j]=Delta_PT(*(Sx+i-1),*(Sz+j-1),flag);
				kernel_mat[0]+=DPS[i][j];
			}
			
			else
				DPS[i][j]=0;

	for(l=1;l<p;l++){
		kernel_mat[l]=0;
		for (j=0; j<=m; j++) 
			DP[l-1][j]=0;
		for (i=0; i<=n; i++) 
			DP[i][l-1]=0;
		
		for (i=l; i<=n; i++)
			for (j=l; j<=m; j++){
				DP[i][j] = DPS[i][j]+LAMBDA*DP[i-1][j]
						+ LAMBDA*DP[i][j-1]
						- LAMBDA2*DP[i-1][j-1];
				
				if( strcmp((*(Sx+i-1))->sName,(*(Sz+j-1))->sName)==0  )
				{
					 // printf("dierci\n");
						DPS[i][j] = Delta_PT(*(Sx+i-1),*(Sz+j-1),flag)* DP[i-1][j-1];
						kernel_mat[l] += DPS[i][j];
				}
				
						// else DPS[i][j] = 0;             
			}
				//      printf("\n----------------------------------\n"); printf("DPS i:%d, j:%d, l:%d\n",n,m,l+1);stampa_math(DPS,n,m);printf("DP\n");stampa_math(DP,n,m); 
	}
				
	//  K=kernel_mat[p-1];
	K=0;
	for(l=0;l<p;l++){
		K+=kernel_mat[l];
		//printf("String kernel of length %d: %1.7f \n\n",l+1,kernel_mat[l]);
    }
    return K;
	
	//   printf("\nDPS\n"); stampa_math(DPS,n,m); printf("DP\n");  stampa_math(DP,n,m); 
	
	//   printf("kernel: n=%d m=%d, %s %s \n\n",n,m,(*(Sx))->sName,(*(Sz))->sName);
	
	
}

#endif
*/
#ifndef FAST


void stampa_math(double *DPS,int n,int m){
	int i,j;
	
	printf("\n");  
	for (i=0; i<=n; i++){
		for (j=0; j<=m; j++)
			printf("%1.8f\t",DPS(i,j));
		printf("\n");  
	}
  printf("\n");  
}

// SLOW SOLUTION BUT ABLE TO DEAL WITH MORE DATA 

double Delta_SK(TreeNode **Sx, TreeNode ** Sz, int n, int m){
	
	
  double *DPS =(double*) malloc((m+1)*(n+1)*sizeof(double));
  double *DP = (double*) malloc((m+1)*(n+1)*sizeof(double));
  double *kernel_mat = (double*) malloc((n+1)*sizeof(double));
  
  int i,j,l,p;
  double K;
	
  p = n; if (m<n) p=m;if (p>MAX_CHILDREN) p=MAX_CHILDREN;
	
	//  if(n==0 || m==0 || m!=n) return 0;
  
  for (j=0; j<=m; j++)
		for (i=0; i<=n; i++) DPS(i,j) = DP(i,j) =0;
	
  
	//printf("\nDPS(%d,%d)\n",n,m); fflush(stdout);
	//stampa_math(DPS,n,m); fflush(stdout);
	
  kernel_mat[0]=0;
  for (i=1; i<=n; i++)
		for (j=1; j<=m; j++)
			if(strcmp((*(Sx+i-1))->sName,(*(Sz+j-1))->sName)==0) 
			{
				DPS(i,j)=Delta_PT(*(Sx+i-1),*(Sz+j-1));
				kernel_mat[0]+=DPS(i,j);
			}
			else DPS(i,j)=0;
	
	
	//  printf("\n\nDPS(%d,%d)\n",n,m); fflush(stdout);
	//  stampa_math(DPS,n,m); fflush(stdout);
	//  printf("\n\nDP(%d,%d)\n",n,m);  fflush(stdout);
	//  stampa_math(DPS,n,m); fflush(stdout);
	//  printf("\n\nKernel: n=%d m=%d, %s %s \n\n",n,m,(*(Sx))->sName,(*(Sz))->sName);fflush(stdout);
	
	for(l=1;l<p;l++){
		kernel_mat[l]=0;
		for (j=0; j<=m; j++)DP(l-1,j)=0;
		for (i=0; i<=n; i++)DP(i,l-1)=0;
		
		for (i=l; i<=n; i++)
			for (j=l; j<=m; j++){
				DP(i,j) = DPS(i,j)+LAMBDA*DP(i-1,j)
				+ LAMBDA*DP(i,j-1)
				- LAMBDA2*DP(i-1,j-1);
				
				if(strcmp((*(Sx+i-1))->sName,(*(Sz+j-1))->sName)==0){
					DPS(i,j) = Delta_PT(*(Sx+i-1),*(Sz+j-1))* DP(i-1,j-1);
					kernel_mat[l] += DPS(i,j);
				}
				// else DPS[i][j] = 0;             
			}
		//      printf("\n----------------------------------\n"); printf("DPS i:%d, j:%d, l:%d\n",n,m,l+1);stampa_math(DPS,n,m);printf("DP\n");stampa_math(DP,n,m); 
	}
	//  K=kernel_mat[p-1];
	K=0;
	for(l=0;l<p;l++){K+=kernel_mat[l];
		//printf("String kernel of legnth %d: %1.7f \n\n",l+1,kernel_mat[l]);
	}
	
  
  free(kernel_mat);
  free(DPS);
  free(DP);
  
  return K;
}

#endif

// DELTA FUNCTION Moschitti's Partial Tree

double Delta_PT( TreeNode * Nx, TreeNode * Nz,int flag){
	double sum=0;
	//printf("flag= %d\n",flag);
	//flag=1;
	
	if(delta_matrix[Nx->nodeID][Nz->nodeID]!=-1) {
		//printf("similarity(%s,%s)is %f\n",Nx->sName,Nz->sName,delta_matrix[Nx->nodeID][Nz->nodeID]);
	  return delta_matrix[Nx->nodeID][Nz->nodeID]; // already there
	}
	
	if(strcmp(Nx->sName,Nz->sName)!=0)
		return (delta_matrix[Nx->nodeID][Nz->nodeID]=0);
	
	else if(Nx->iNoOfChildren==0 || Nz->iNoOfChildren==0){
		if(flag) //PT remove leaves
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=0);
		else  //PT
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=MU*LAMBDA2);
	}
	
	else{
		sum = MU*(LAMBDA2+Delta_SK(Nx->pChild, Nz->pChild,Nx->iNoOfChildren, Nz->iNoOfChildren,flag)); 
		//sum=1;
		//printf("\n sum = %f \n",sum);
		
		 //            printf("\n (node1:%s node2:%s) -----------------------> %1.30f \n\n\n",Nx->sName,Nz->sName, sum);
		return (delta_matrix[Nx->nodeID][Nz->nodeID]=sum);
	}
	
	return 0;
}

/*
double Delta_PT( TreeNode * Nx, TreeNode * Nz,int flag){
	double sum=0;
	double similarity = 1.0;
	flag = 1;

	//printf("\nwbcs:delta_matrix = %f ,Nx(%s)->ID = %d ,Nz(%s)->ID = %d\n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nx->nodeID,Nz->sName,Nz->nodeID);
	if(delta_matrix[Nx->nodeID][Nz->nodeID]!=-1) {
		if( (delta_matrix[Nx->nodeID][Nz->nodeID]==0) &&(Nx->iNoOfChildren==0 && Nz->iNoOfChildren==0)){
			//leaves but not same,then use the word2vec similarity
				//use the word2vec similarity
			//printf("\ndelta_matrix = %f ,Nx(%s)->ID = %d ,Nz(%s)->ID = %d\n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nx->nodeID,Nz->sName,Nz->nodeID);
			char buf1[50] ="null";
			char buf2[50] ="null";
			sscanf(strlwr(Nx->sName),"%[^a-z]",buf1);
			sscanf(strlwr(Nz->sName),"%[^a-z]",buf2);
			//printf("\nbuf1=%s,buf2=%s",buf1,buf2);
			//printf("\njinlaile %s %s %s %s\n",Nx->sName,Nz->sName,buf1,buf2);
			if(strcmp(buf1,Nx->sName)==0 || strcmp(buf2,Nz->sName)==0 ){
				//printf("\n两个都不包含字母 %s %s\n",Nx->sName,Nz->sName);
				return (delta_matrix[Nx->nodeID][Nz->nodeID]);
			}
			else if(strcmp(strlwr(Nx->sName),strlwr(Nz->sName))!=0)
			{
				//sscanf(strlwr(Nx->sName),"%[a-zA-Z0-9\u4e00-\u9fa5\s]",buf1);
			    //sscanf(strlwr(Nz->sName),"%[a-zA-Z0-9\u4e00-\u9fa5\s]",buf2);
				//Replace(yourStr, @"[a-zA-Z0-9\u4e00-\u9fa5\s]", "");
				similarity = getW2V_similarity(strlwr(Nx->sName),strlwr(Nz->sName));
				//printf("\n%s和%s 相似度是%f\n",Nx->sName,Nz->sName,similarity);
				return (delta_matrix[Nx->nodeID][Nz->nodeID]=MU*LAMBDA2*similarity);
			}
			else{
				return (delta_matrix[Nx->nodeID][Nz->nodeID]=MU*LAMBDA2);
			}
		}
		else
			return delta_matrix[Nx->nodeID][Nz->nodeID]; // already there
	}
	
	if(strcmp(Nx->sName,Nz->sName)!=0){
		//printf(" 两个不一样啦！！\n");
		if(Nx->iNoOfChildren==0 && Nz->iNoOfChildren==0){
			char buf1[50] ="null";
			char buf2[50] ="null";
			sscanf(strlwr(Nx->sName),"%[^a-z]",buf1);
			sscanf(strlwr(Nz->sName),"%[^a-z]",buf2);
			//printf("\nbuf1=%s,buf2=%s",buf1,buf2);
			//printf("\njinlaile %s %s %s %s\n",Nx->sName,Nz->sName,buf1,buf2);
			if(strcmp(buf1,Nx->sName)==0 || strcmp(buf2,Nz->sName)==0 ){
				//printf("\n两个都不包含字母 %s %s\n",Nx->sName,Nz->sName);
				if(strcmp(Nx->sName,Nz->sName)==0)
				  return (delta_matrix[Nx->nodeID][Nz->nodeID]=MU*LAMBDA2);
				else
				  return (delta_matrix[Nx->nodeID][Nz->nodeID]=0);
			}
			else if(strcmp(strlwr(Nx->sName),strlwr(Nz->sName))!=0)
			{
				//sscanf(strlwr(Nx->sName),"%[a-zA-Z0-9\u4e00-\u9fa5\s]",buf1);
			    //sscanf(strlwr(Nz->sName),"%[a-zA-Z0-9\u4e00-\u9fa5\s]",buf2);
				//Replace(yourStr, @"[a-zA-Z0-9\u4e00-\u9fa5\s]", "");
				similarity = getW2V_similarity(strlwr(Nx->sName),strlwr(Nz->sName));
				//printf("\n%s和%s 相似度是%f\n",Nx->sName,Nz->sName,similarity);
				return (delta_matrix[Nx->nodeID][Nz->nodeID]=MU*LAMBDA2*similarity);
			}
			else
				return (delta_matrix[Nx->nodeID][Nz->nodeID]=MU*LAMBDA2);
		}
		else
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=0);
	}
	
	else if(Nx->iNoOfChildren==0 || Nz->iNoOfChildren==0){
		if(Nx->iNoOfChildren==0 && Nz->iNoOfChildren==0){
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=MU*LAMBDA2);
		}
		else if(flag){ //PT remove leaves
			//printf("\nxiaci:delta_matrix = %f ,Nx(%s)->ID = %d ,Nz(%s)->ID = %d\n",0,Nx->sName,Nx->nodeID,Nz->sName,Nz->nodeID);
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=0);
		}
		else{  //PT
			//printf("\ndelta_matrix = %f ,Nx(%s)->ID = %d ,Nz(%s)->ID = %d\n",MU*LAMBDA2,Nx->sName,Nx->nodeID,Nz->sName,Nz->nodeID);
			return (delta_matrix[Nx->nodeID][Nz->nodeID]=MU*LAMBDA2);
		}
	}
	else{
		sum = MU*(LAMBDA2+Delta_SK(Nx->pChild, Nz->pChild,Nx->iNoOfChildren, Nz->iNoOfChildren,flag)); 
		return (delta_matrix[Nx->nodeID][Nz->nodeID]=sum);
	}
	return 0;
}
*/
/*
double Delta_PT( TreeNode * Nx, TreeNode * Nz,int flag){
	double sum=0;
	double similarity = 1.0;
	flag = 0;

	//printf("\nwbcs:delta_matrix = %f ,Nx(%s)->ID = %d ,Nz(%s)->ID = %d\n",delta_matrix[Nx->nodeID][Nz->nodeID],Nx->sName,Nx->nodeID,Nz->sName,Nz->nodeID);
	if(delta_matrix[Nx->nodeID][Nz->nodeID]!=1) {
		return delta_matrix[Nx->nodeID][Nz->nodeID]; // already there
	}

	else if(Nx->iNoOfChildren==0 || Nz->iNoOfChildren==0){
		return delta_matrix[Nx->nodeID][Nz->nodeID];
	}
	else{
		sum = MU*(LAMBDA2+Delta_SK(Nx->pChild, Nz->pChild,Nx->iNoOfChildren, Nz->iNoOfChildren,flag)); 
		return (delta_matrix[Nx->nodeID][Nz->nodeID]=sum);
	}
	return 0;
}
*/
double hasdigit(char * Nx_sname,char * Nz_sname){//先对两个单词去掉符号，再比较如果是纯数字，则直接字符串匹配，否则使用词向量

	/*for(int i=0;i<strlen(word);i++){
		if((word[i]>='0')&& (word[i]<='9'))
			return 1;
	}
	return 0;
	*/
	vector<char> sx;
	vector<char> sz;
	vector<char>::iterator it;
	for(int i=0;i<strlen(Nx_sname);i++)
		sx.push_back(Nx_sname[i]);
	for(int i=0;i<strlen(Nz_sname);i++)
		sz.push_back(Nz_sname[i]);
	//先对两个词去除符号
	sx.erase(remove_if(sx.begin(),sx.end(),static_cast<int(*)(int)>(&ispunct) ),sx.end());
	sz.erase(remove_if(sz.begin(),sz.end(),static_cast<int(*)(int)>(&ispunct) ),sz.end()); 
	int i=0;
	int sxsize = sx.size();
	int szsize = sz.size();
	char *s1=new char[sxsize+1];
	char *s2=new char[szsize+1];
	//printf(" sxsize=%d,szsize=%d ",sxsize,szsize);
	//对两个词去除数字
	for(it=sx.begin();it!=sx.end();++it){
		if((*it>='0')&& (*it<='9')){
			return 1;
		}else{
			s1[i]=*it;
			i++;
		}
	}
	s1[i]='\0';
	int j=0;
	for(it=sz.begin();it!=sz.end();++it){
		if((*it>='0')&& (*it<='9')){
			return 1;
		}else{
			s2[j]=*it;
			j++;
		}
	}
	s2[j]='\0';
	
	//printf("\n去掉符号后 Nx:%s Nz:%s",s1,s2);
	if((i<=1)||(j<=1)){//去除符号和数字后有一个词长度小于1.
		return 1;
	}
	return 0;
}

double getW2V_similarity(char * Nx_sname,char *Nz_sname){//先对两个单词去掉符号，再比较如果是纯数字，则直接字符串匹配，否则使用词向量
	vector<char> sx;
	vector<char> sz;
	vector<char>::iterator it;
	for(int i=0;i<strlen(Nx_sname);i++)
		sx.push_back(Nx_sname[i]);
	for(int i=0;i<strlen(Nz_sname);i++)
		sz.push_back(Nz_sname[i]);
	//先对两个词去除符号
	sx.erase(remove_if(sx.begin(),sx.end(),static_cast<int(*)(int)>(&ispunct) ),sx.end());
	sz.erase(remove_if(sz.begin(),sz.end(),static_cast<int(*)(int)>(&ispunct) ),sz.end()); 
	int i=0;
	int sxsize = sx.size();
	int szsize = sz.size();
	char *s1=new char[sxsize+1];
	char *s2=new char[szsize+1];
	//printf(" sxsize=%d,szsize=%d ",sxsize,szsize);
	//对两个词去除数字
	for(it=sx.begin();it!=sx.end();++it){
		if((*it>='0')&& (*it<='9')){
			//return 0.0;
			continue;
		}else{
			s1[i]=*it;
			i++;
		}
	}
	s1[i]='\0';
	int j=0;
	for(it=sz.begin();it!=sz.end();++it){
		if((*it>='0')&& (*it<='9')){
			//return 0.0;
			continue;
		}else{
			s2[j]=*it;
			j++;
		}
	}
	s2[j]='\0';
	
	//printf("\n去掉符号后 Nx:%s Nz:%s",s1,s2);
	if((i<=1)||(j<=1)){//去除符号和数字后有一个词长度小于1.
		if(i==0||j==0){//对比的两个词至少有一个词只包含数字或符号
			if(strcmp(strlwr(Nx_sname),strlwr(Nz_sname))==0)
				return 1.0;
			else
				return 0.0;
		}

		else if(strcmp(strlwr(s1),strlwr(s2))!=0){//对比的两个词至少有一个词只包含一个字母
			delete(s1);
			delete(s2);
			return 0.0;
	    }

		else{
			delete(s1);
			delete(s2);
			return 1.0;
		}
	}

	else if(strcmp(strlwr(s1),strlwr(s2))==0){
		delete(s1);
		delete(s2);
		return 1.0;
	}

	else
	{
		//printf("\n去掉符号后不相同");
		FILE *fl = fopen("new_dict_50_30_test_STS.input.SMTeuroparl.stopwords.txt","r");
		int state1=0;
		int state2=0;
		int size=0;

		if (fl == NULL) {
			printf("Embedding file not found\n");
			return 0.0;
		}

		fscanf(fl,"%d",&size);
		char str[100];
		double *v1 = new double [size];
		double *v2 = new double [size];
		double temp;

		//for(int i=0;i<wordNum;i++){
		while(fscanf(fl,"%s",str)!=EOF){

			if(state1&&state2){
				break;
			}

			if(strcmp(str,strlwr(s1))==0 ){
				//printf("\n找到了%s v = \n",Nx_sname);
				for(int a=0;a<size;a++)
					fscanf(fl,"%lf",&v1[a]);
				state1=1;
			}
			else if(strcmp(str,strlwr(s2))==0){
				//printf("\n找到了第二个词%s v = \n",Nz_sname);
				state2=1;
				for(int a=0;a<size;a++)
					fscanf(fl,"%lf",&v2[a]);
			}
			else{
				for(int a=0;a<size;a++)
					fscanf(fl,"%lf",&temp);
			}
		}
		fclose(fl);
		if(state1&&state2){
			//printf("在小词典找到了！\n");
            double t1 = 0, t2 = 0, t3 = 0;
			double dist=0;
			for (int i = 0; i < size; i++) {
				//printf(" ");
				t1 += v1[i] * v2[i];
				t2 += v1[i] * v1[i];
				t3 += v2[i] * v2[i];
				
				//dist += (v1[i]-v2[i])*(v1[i]-v2[i]);

			}
			double sim = t1 / sqrt(t2 + eps) / sqrt(t3 + eps);
			//double sim = t1 / sqrt(t2) / sqrt(t3);
			//double sim = (double(1.0))/((double(1.0))+sqrt(dist));
		//	sim = double(0.5)+double(0.5)*sim;
			if(sim<=0.0)
				sim=0.0;
			
			printf("\n%s和%s的相似度是%f\n",Nx_sname,Nz_sname,sim);
			delete(v1);
			delete(v2);
			delete(s1);
			delete(s2);
			return sim;
		}
		else
			return 0.0;
	}
}
/*
	double getW2V_similarity(char * Nx_sname,char *Nz_sname){
		//printf("jinlaileyo\n");
		if(strcmp(Nx_sname,Nz_sname)==0)
			return MU*LAMBDA2;
		else
			return 0.0;
    }
*/

/*
double getW2V_similarity(char * Nx_sname,char *Nz_sname){//先对两个单词去掉符号
	vector<char> sx;
	vector<char> sz;
	vector<char>::iterator it;
	for(int i=0;i<strlen(Nx_sname);i++)
		sx.push_back(Nx_sname[i]);
	for(int i=0;i<strlen(Nz_sname);i++)
		sz.push_back(Nz_sname[i]);
	sx.erase(remove_if(sx.begin(),sx.end(),static_cast<int(*)(int)>(&ispunct) ),sx.end());
	sz.erase(remove_if(sz.begin(),sz.end(),static_cast<int(*)(int)>(&ispunct) ),sz.end()); 
	int i=0;
	int sxsize = sx.size();
	int szsize = sz.size();
	char *s1=new char[sxsize+1];
	char *s2=new char[szsize+1];
	//printf(" sxsize=%d,szsize=%d ",sxsize,szsize);
	for(it=sx.begin();it!=sx.end();++it){
		if(!((*it>='0')&& (*it<='9'))){
			s1[i]=*it;
			i++;
		}
	}
	s1[i]='\0';
	int j=0;
	for(it=sz.begin();it!=sz.end();++it){
		if(!((*it>='0')&& (*it<='9'))){
			s2[j]=*it;
			j++;
		}
	}
	s2[j]='\0';
	
	//printf("\n去掉符号后 Nx:%s Nz:%s",s1,s2);
	if((i==1)||(j==1)){
		if(strcmp(s1,s2)!=0){
			delete(s1);
			delete(s2);
			return 0.0;
	    }
		else{
			delete(s1);
			delete(s2);
			return 1.0;
		}
	}

	else if(strcmp(s1,s2)==0){
		delete(s1);
		delete(s2);
		return 1.0;
	}

	else
	{//先从小字典找，找到了就返回，若没找到就从大字典找，并写入小字典。
		//printf("\n去掉符号后不相同");
		FILE *fl = fopen("dict_test_STS.input.MSRpar.txt","rb");
		int state=0;
		int state1=0;
		int state2=0;
		int size=0;

		if (fl == NULL) {
			printf("Embedding file not found\n");
			return 0.0;
		}

		fscanf(fl,"%d",&size);
		char str[100];
		double *v1 = new double [size];
		double *v2 = new double [size];
		double temp;

		//for(int i=0;i<wordNum;i++){
		while(fscanf(fl,"%s",str)!=EOF){

			if(state1&&state2){
				break;
			}

			if(strcmp(str,s1)==0 ){
				//printf("\n找到了%s v = \n",Nx_sname);
				fscanf(fl,"%lf",&temp);
				if(temp==-111){//没有对应词向量
					state=1;
					break;
				}
				else{
					v1[0]=temp;
					for(int a=1;a<size;a++){
						fscanf(fl,"%lf",&v1[a]);
						//printf(" %lf",v1[a]);
					}
				}
				//printf("\n");
				state1=1;
			}
			if(strcmp(str,s2)==0){
				//printf("\n找到了第二个词%s v = \n",Nz_sname);
				
				state2=1;
				fscanf(fl,"%lf",&temp);
				if(temp==-111){//没有对应词向量
					state=1;
					break;
				}
				else{
					v2[0]=temp;
					for(int a=1;a<size;a++){
						fscanf(fl,"%lf",&v2[a]);
						//printf(" %lf",v2[a]);
					}
				}
				//printf("\n");
			}
			else{
				for(int a=0;a<size;a++)
					fscanf(fl,"%lf",&temp);
			}
		}
		fclose(fl);
		if(state==1){
			delete(v1);
			delete(v2);
			delete(s1);
			delete(s2);
			return 0.0;
		}

		else if(state1&&state2){
			//printf("在小词典找到了！\n");
            double t1 = 0, t2 = 0, t3 = 0;
			for (int i = 0; i < size; i++) {
				//printf(" ");
				t1 += v1[i] * v2[i];
				t2 += v1[i] * v1[i];
				t3 += v2[i] * v2[i];
			}
			double sim = t1 / sqrt(t2 + eps) / sqrt(t3 + eps);
			sim = double(0.5)+double(0.5)*sim;
			//printf("\n%s和%s的相似度是%f\n",Nx_sname,Nz_sname,sim);
			delete(v1);
			delete(v2);
			delete(s1);
			delete(s2);
			return sim;
		}
		else{
			FILE *fp = fopen("enwiki_cbow_200_8_50.txt","rb");
			size=0;
			if (fp == NULL) {
					printf("Embedding file not found\n");
					return 0.0;
			}
			int wordNum=0;
			fscanf(fp,"%d",&wordNum);
			fscanf(fp,"%d",&size);
			//printf("大辞典size = %d\n",size);
			if((state1==0) && (state2!=0)){
				//printf("\n第一个词%s在小词典没找到！ \n",Nx_sname);
				char str1[100];
				while(fscanf(fp,"%s",str1)!=EOF){

					if(strcmp(str1,s1)==0 ){
						//printf("\n找到了%s v = \n",Nx_sname);
						for(int a=0;a<size;a++){
							fscanf(fp,"%lf",&v1[a]);
						}
						state1=1;
						//写小文件
						FILE *fout = fopen("dict_test_STS.input.MSRpar.txt","a+");
						if(fout==NULL)
							;
						else{
							fprintf(fout,"%s",s1);
							for(int a=0;a<size;a++)
								fprintf(fout," %lf",v1[a]);
							fprintf(fout,"\n");

							fclose(fout);
						}
						//计算相似度
						double t1 = 0, t2 = 0, t3 = 0;
						for (int i = 0; i < size; i++) {
							t1 += v1[i] * v2[i];
							t2 += v1[i] * v1[i];
							t3 += v2[i] * v2[i];
						}
						double sim = t1 / sqrt(t2 + eps) / sqrt(t3 + eps);
						sim = double(0.5)+double(0.5)*sim;
						//printf("\n大辞典 %s和%s的相似度是%f\n",Nx_sname,Nz_sname,sim);
						delete(v1);
						delete(v2);
						delete(s1);
						delete(s2);
						fclose(fout);
						fclose(fp);
						return sim;
					}
					else{
						for(int a=0;a<size;a++)
							fscanf(fp,"%lf",&temp);
					}
				}
				fclose(fp);
				//printf("%s在大词典也没找到！\n",Nx_sname);
				FILE *fout = fopen("dict_test_STS.input.MSRpar.txt","a+");
				if(fout==NULL)
					return 0.0;
				else{
					fprintf(fout,"%s",s1);
					temp=-111;
					fprintf(fout," %lf\n",temp);
					
					fclose(fout);
					return 0.0;
				}
			}
			else if((state2==0)&&(state1!=0)){
				//printf("\n第二个词%s在小词典没找到！ \n",Nz_sname);
				char str2[100];
				while(fscanf(fp,"%s",str2)!=EOF){
					if(strcmp(str2,s2)==0 ){
						//printf("\n找到了%s v = \n",Nx_sname);
						for(int a=0;a<size;a++){
							fscanf(fp,"%lf",&v2[a]);
						}
						state1=1;
						//写小文件
						FILE *fout = fopen("dict_test_STS.input.MSRpar.txt","a+");
						if(fout==NULL)
							;
						else{
							fprintf(fout,"%s",s2);
							//printf("第二个词%s\n",s2);
							for(int a=0;a<size;a++)
								fprintf(fout," %lf",v2[a]);
							fprintf(fout,"\n");

							fclose(fout);
						}
						//计算相似度
						double t1 = 0, t2 = 0, t3 = 0;
						for (int i = 0; i < size; i++) {
							t1 += v1[i] * v2[i];
							t2 += v1[i] * v1[i];
							t3 += v2[i] * v2[i];
						}
						double sim = t1 / sqrt(t2 + eps) / sqrt(t3 + eps);
						sim = double(0.5)+double(0.5)*sim;
						//printf("\n大辞典 %s和%s的相似度是%f\n",Nx_sname,Nz_sname,sim);
						delete(v1);
						delete(v2);
						delete(s1);
						delete(s2);
						fclose(fout);
						fclose(fp);
						return sim;
					}
					else{
						for(int a=0;a<size;a++)
							fscanf(fp,"%lf",&temp);
					}
				}
				fclose(fp);
				//printf("%s在大词典也没找到！\n",Nz_sname);
				FILE *fout = fopen("dict_test_STS.input.MSRpar.txt","a+");
				if(fout==NULL)
					return 0.0;
				else{
					fprintf(fout,"%s",s2);
					temp=-111;
					fprintf(fout," %lf\n",temp);
					
					fclose(fout);
					return 0.0;
				}
			}
			else{
				//printf("\n两个词%s和%s在小词典没找到！ \n",Nx_sname,Nz_sname);
				if (fp == NULL) {
					printf("Embedding file not found\n");
					//return 0.0;
				}
				fscanf(fp,"%d",&size);
				char str3[100];

				while(fscanf(fp,"%s",str3)!=EOF){
					if(state1&&state2){
						//写小文件
						FILE *fout = fopen("dict_test_STS.input.MSRpar.txt","a+");
						if(fout==NULL)
							;
						else{
							fprintf(fout,"%s",s1);
							for(int a=0;a<size;a++)
								fprintf(fout," %lf",v1[a]);
							fprintf(fout,"\n");

							fprintf(fout,"%s",s2);
							for(int a=0;a<size;a++)
								fprintf(fout," %lf",v2[a]);
							fprintf(fout,"\n");
					
							fclose(fout);
						}
						//计算相似度
						double t1 = 0, t2 = 0, t3 = 0;
						for (int i = 0; i < size; i++) {
							t1 += v1[i] * v2[i];
							t2 += v1[i] * v1[i];
							t3 += v2[i] * v2[i];
						}
						double sim = t1 / sqrt(t2 + eps) / sqrt(t3 + eps);
						sim = double(0.5)+double(0.5)*sim;
						//printf("\n大辞典 %s和%s的相似度是%f\n",Nx_sname,Nz_sname,sim);
						delete(v1);
						delete(v2);
						delete(s1);
						delete(s2);
						fclose(fout);
						fclose(fp);
						return sim;
					}

					if(strcmp(str3,s1)==0 ){
						//printf("\n找到了%s v = \n",Nx_sname);
						for(int a=0;a<size;a++){
							fscanf(fp,"%lf",&v1[a]);
							//printf(" %lf",v1[a]);
						}
						state1=1;
					}
					else if(strcmp(str3,s2)==0 ){
						//printf("\n找到了%s v = \n",Nx_sname);
						for(int a=0;a<size;a++){
							fscanf(fp,"%lf",&v2[a]);
							//printf(" %lf",v1[a]);
						}
						state2=1;
					}
					else{
						for(int a=0;a<size;a++)
							fscanf(fp,"%lf",&temp);
					}
				}
				fclose(fp);
				if(state1==0 && state2 ==1){
					//printf("在大词典也没找到！\n");
					FILE *fout = fopen("dict_test_STS.input.MSRpar.txt","a+");
					if(fout==NULL)
						return 0.0;
					else{
						fprintf(fout,"%s",s1);
						temp=-111;
						fprintf(fout," %lf\n",temp);

						fclose(fout);
						return 0.0;
					}
				}
				else if(state2==0 && state1==1){
					//printf("在大词典也没找到！\n");
					FILE *fout = fopen("dict_test_STS.input.MSRpar.txt","a+");
					if(fout==NULL)
						return 0.0;
					else{
						temp = -111;
						fprintf(fout,"%s",s2);
						fprintf(fout," %lf\n",temp);

						fclose(fout);
						return 0.0;
					}
				}
				else{
				//	printf("在大词典也没找到！\n");
					FILE *fout = fopen("dict_test_STS.input.MSRpar.txt","a+");
					if(fout==NULL)
						return 0.0;
					else{
						temp=-111;
						fprintf(fout,"%s",s1);
						fprintf(fout," %lf\n",temp);

						fprintf(fout,"%s",s2);
						fprintf(fout," %lf\n",temp);

						fclose(fout);
						return 0.0;
					}
				}
			}
		}
		delete(v1);
		delete(v2);
		delete(s1);
		delete(s2);
	}
}
*/
double evaluate_PT(nodePair *pairs, int n,int flag){
	
	int i;
	double sum=0,contr;
	//printf("%d\n",flag);
	for(i=0;i<n;i++){
		
		//contr=Delta_PT(pairs[i].Nx,pairs[i].Nz,"null","null",flag);//-REMOVE_LEAVES; // remove the contribution of leaves
			contr=Delta_PT(pairs[i].Nx,pairs[i].Nz,flag);																								 //printf("Score: %1.15f\n",contr); fflush(stdout);
		//if(pairs[i].Nx->iNoOfChildren!=0 || pairs[i].Nz->iNoOfChildren!=0)
		//	contr-=REMOVE_LEAVES;
		
		// exit(-1);
		sum+=contr;
	}
	return sum;	  
}

/*
void determine_sub_lists(FOREST *a, FOREST *b, nodePair *intersect, int *n){

   int i=0,j=0,j_old,j_final;
   int n_a,n_b;
   short cfr;
   OrderedTreeNode *list_a, *list_b;

   n_a = a->listSize;
   n_b = b->listSize;
   list_a=a->orderedNodeSet; 
   list_b=b->orderedNodeSet;
   *n=0;

   //先计算好delta_matrix矩阵，然后根据矩阵获取相似点对
   //计算delta_matrix矩阵，如果是叶子则调用词向量的方法计算，否则使用硬匹配
   for(i=0;i<n_a;i++){
	   if(list_a[i].node->pre_terminal==-1){
		   for(j=0;j<n_b;j++){
			   if(list_b[j].node->pre_terminal==-1){ //都是叶子单词,则使用词向量计算
				 /*  if(strcmp(strlwr(list_a[i].sName),strlwr(list_b[j].sName))==0){
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=1.0;
					   printf("\n%s和%s的相似度是%f\n",list_a[i].sName,list_b[j].sName,1.0);
						intersect[*n].Nx=list_a[i].node;
						intersect[*n].Nz=list_b[j].node;
						(*n)++;
				   }
				   else{
					   //先判断是否为反义词，然后再使用词向量方法
					 /*  int ia=isantonym(strlwr(list_a[i].sName),strlwr(list_b[j].sName));
					   if(ia==1) //是反义词
						   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=0.0;
					   else{*/
/*						   double w2v_sim=getW2V_similarity(strlwr(list_a[i].sName),strlwr(list_b[j].sName));
						   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=w2v_sim;//*LAMBDA2;
						   printf("%s和%s的相似度=%f\n",list_a[i].sName,list_b[j].sName,w2v_sim);
						   if(w2v_sim>0.0){
							   intersect[*n].Nx=list_a[i].node;
							   intersect[*n].Nz=list_b[j].node;
							   (*n)++;
						   }
					   //}
			    }
		   }
	   }
	   else{
		   for(j=0;j<n_b;j++){
			   if(list_b[j].node->pre_terminal!=-1){
				   if(strcmp(list_a[i].sName,list_b[j].sName)==0){
					   intersect[*n].Nx=list_a[i].node;
					   intersect[*n].Nz=list_b[j].node;
					   (*n)++;
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=1;
				   }
				   else
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=0.0;
			   }
		   }
	  }
   }
}
*/


void determine_sub_lists(FOREST *a, FOREST *b, nodePair *intersect, int *n){

   int i=0,j=0,j_old,j_final;
   int n_a,n_b;
   short cfr;
   OrderedTreeNode *list_a, *list_b;

   n_a = a->listSize;
   n_b = b->listSize;
   list_a=a->orderedNodeSet; 
   list_b=b->orderedNodeSet;
   *n=0;
   

   //compute delta_matrix first,then find the similar pairs.for the leaves, use word embedding to measure similarity
   float delta_sim[MAX_NUMBER_OF_NODES][MAX_NUMBER_OF_NODES];
  // string leaf_id[MAX_NUMBER_OF_NODES];
  // int size=0;
   for(i=0;i<n_a;i++){
	   if(list_a[i].node->pre_terminal==-1){
		   for(j=0;j<n_b;j++){
			   if(list_b[j].node->pre_terminal==-1){ //both leaves,use word vector
				   if(strcmp(strlwr(list_a[i].sName),strlwr(list_b[j].sName))==0){//string hard mathcing first, if not same,to word vector dictionary
					   intersect[*n].Nx=list_a[i].node;
					   intersect[*n].Nz=list_b[j].node;
					   (*n)++;
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=1;
				   }
				   else{
					   if(hasdigit(list_a[i].sName,list_b[j].sName)==1){
						    if(strcmp(strlwr(list_a[i].sName),strlwr(list_b[j].sName))==0){//string hard mathcing first, if not same,to word vector dictionary
							   intersect[*n].Nx=list_a[i].node;
							   intersect[*n].Nz=list_b[j].node;
							   (*n)++;
							   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=1;
							}
					   }
					   else{
						    dict[list_a[i].sName]=NULL;
							dict[list_b[j].sName]=NULL;
					   }
				   }
			    }
		   }
	   }
	   else{
		   for(j=0;j<n_b;j++){
			   if(list_b[j].node->pre_terminal!=-1){
				   if(strcmp(list_a[i].sName,list_b[j].sName)==0){
					   intersect[*n].Nx=list_a[i].node;
					   intersect[*n].Nz=list_b[j].node;
					   (*n)++;
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=1;
				   }
				   else
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=0.0;
			   }
		   }
	  }
   }
   
   int size=ReadEmbedding("new_enwiki_50_30.txt");
 
   if(size>0){
	   for(i=0;i<n_a;i++){
		   if(list_a[i].node->pre_terminal==-1){
			   for(j=0;j<n_b;j++){
				   if(list_b[j].node->pre_terminal==-1){ //both leaves,use word vector
					  // printf("jinlaile\n");
					   //if((isantonym(list_a[i].sName,list_b[j].sName)==0)){//if not antonym，cosine similarity
						   if((dict[list_a[i].sName]!=NULL)&&(dict[list_b[j].sName]!=NULL)){
							 //  printf("jinlaile\n");
						   
							   double t1 = 0, t2 = 0, t3 = 0;
							   double *v1 = dict[list_a[i].sName];
							   double *v2 = dict[list_b[j].sName];
							   double w2v_sim=0.0;
								for (int i = 0; i < size; i++) {
									t1 += v1[i] * v2[i];
									t2 += v1[i] * v1[i];
									t3 += v2[i] * v2[i];
								}
								// printf("jinlaile\n");
								w2v_sim = t1 / sqrt(t2 + eps) / sqrt(t3 + eps);
								// printf("jinlaile\n");
								//printf("%s和%s similarity is：%f\n",list_a[i].sName,list_b[j].sName,w2v_sim);
							
								if(w2v_sim>0.0){
								   intersect[*n].Nx=list_a[i].node;
								   intersect[*n].Nz=list_b[j].node;
								   (*n)++;
								   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=w2v_sim;//*LAMBDA2;
								}
						   }
					 //  }
				   }
			   }
		   }
	   }
   }
}

//read embedding file for word vector
int ReadEmbedding(const char *file_name) {
   FILE *f_dict = fopen(file_name,"r");
   if(f_dict==NULL){
	   printf("Embedding file not found\n");
		return -1;
	}
	int size;
	const int MAX_STRING = 1000;
	fscanf(f_dict, "%d", &size);
	//printf("size=%d\n",size);
	char str[MAX_STRING];
	double *tmp = new double[size];
	while(fscanf(f_dict,"%s",str)!=EOF){
		//printf("str=%s\n",str);
		map<string, double*>::iterator it = dict.find(str);
		double *v = tmp;
		if (it != dict.end()) {
			//printf("find it！！%s\n",str);
			if (it->second == NULL) {
				it->second = new double[size];
				v = it->second;
			}

		}
		for (int a = 0; a < size; a++)
			fscanf(f_dict, "%lf", &v[a]);
		
		
	}
	fclose(f_dict);
	
	return size;
}

//find antonym dictionary，if word1 and word2 are antonyms ，return 1，else return 0
int isantonym(char *word1, char *word2){
	FILE *fi = fopen("surprise.SMTnews_antonym.txt","r");
	//printf("ound\n");
	if (fi == NULL) {
		printf("Antonym file not found\n");
		return 0;
	}

	char str[100];
	int num=0;
	while(fscanf(fi,"%s",str)!=EOF){
		
		fscanf(fi,"%d",&num);//printf("ound\n");
		if(strcmp(word1,str)==0){
			for(int i=0;i<num;i++){
				fscanf(fi,"%s",str);
				if(strcmp(str,word2)==0){
					fclose(fi);
					printf("%s and %s is antonym\n",word1,word2);
					return 1;
				}
			}
			fclose(fi);
			//printf("%s and %s no antonym\n",word1,word2);
			return 0;
		}
		else{
			for(int i=0;i<num;i++){
				fscanf(fi,"%s",str);
			}
		}
	}
	fclose(fi);
	//printf("%s and %s no antonym\n",word1,word2);
	return 0;
}

/*
void determine_sub_lists(FOREST *a, FOREST *b, nodePair *intersect, int *n){
	
	int i=0,j=0,j_old,j_final;
	int n_a,n_b;
	short cfr;
	OrderedTreeNode *list_a, *list_b;
	
	n_a = a->listSize;
	n_b = b->listSize;
	list_a=a->orderedNodeSet; 
	list_b=b->orderedNodeSet;
	*n=0;
	
	 
	while(i<n_a && j<n_b){
		if((cfr=strcmp(list_a[i].sName,list_b[j].sName))>0)j++;
		else if(cfr<0)i++;
		else{
			j_old=j;
			do{
				do{
					intersect[*n].Nx=list_a[i].node;
					intersect[*n].Nz=list_b[j].node;
					(*n)++;
					if (*n>MAX_NUMBER_OF_PAIRS) { 
						printf ("ERROR: The number of identical parse nodes exceed the current capacityn\n"); 
						exit(-1);
					}
					delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=-1.0;		    
					//  TEST           printf("Evaluating-Pair: (%s  ,  %s) i %d,j %d j_old%d \n",list_a[i].sName,list_b[j].sName,i,j,j_old);fflush(stdout);
					j++;
				}
				while(j<n_b && strcmp(list_a[i].sName,list_b[j].sName)==0);
				i++;j_final=j;j=j_old;
			} 		
			while(i<n_a && strcmp(list_a[i].sName,list_b[j].sName)==0);
			j=j_final;
		}
	}
	
	
	//printf ("number of pairs  %d \n",*n); 
	
}
*/

/*
void determine_sub_lists(FOREST *a, FOREST *b, nodePair *intersect, int *n){

   int i=0,j=0,j_old,j_final;
   int n_a,n_b;
   short cfr;
   OrderedTreeNode *list_a, *list_b;

   n_a = a->listSize;
   n_b = b->listSize;
   list_a=a->orderedNodeSet; 
   list_b=b->orderedNodeSet;
   *n=0;

   //先计算好delta_matrix矩阵，然后根据矩阵获取相似点对
   //计算delta_matrix矩阵，如果是叶子则调用词向量的方法计算，否则使用硬匹配
   char s1[100];char s2[100];
   for(i=0;i<n_a;i++){
	    //printf("新的一轮开始了\n");
	   if(list_a[i].node->pre_terminal==-1){
		/*   printf("%s\n",list_a[i].sName);
		   int size= strlen(list_a[i].sName);
		   
		   int m=0;
		   for(;m<size;m++)
			   s1[m]=list_a[i].sName[m];
		   s1[m]='\0';
		   
		   char * word1 = strtok(s1,"::");
		   char * pos_word1 = strtok(NULL,"::");*/
/*		   double maxsim=-1;
		   int max_aid=0;
		   int max_bid=0;
		   int state=0;
		   for(j=0;j<n_b;j++){
			   if(list_b[j].node->pre_terminal==-1){
				  // printf("%s %s\n",list_a[i].sName,list_b[j].sName);
				   //先分离叶子和其词性
				   /*int size2= strlen(list_b[j].sName);
				   
				   int t=0;
				   for(;t<size2;t++)
					   s2[t]=list_b[j].sName[t];
				   s2[t]='\0';
				   
				   char * word2 = strtok(s2,"::");
				   char * pos_word2 = strtok(NULL,"::");
				   printf("%s::%s ,%s::%s\n",word1,pos_word1,word2,pos_word2);
				   if(strcmp(strlwr(word1),strlwr(word2))==0){*/
/*				   if(strcmp(strlwr(list_a[i].sName),strlwr(list_b[j].sName))==0){
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=1.0;
						intersect[*n].Nx=list_a[i].node;
						intersect[*n].Nz=list_b[j].node;
						(*n)++;
						state=1;
						//printf("找到和自己相同的词了：%s \n",list_a[i].sName);
						//free(s1);
						//free(s2);
						//delete(s1);
						//delete(s2);
						printf("找到和自己相同的词了lalala：%s \n",list_a[i].sName);
						break;
				   }
				   else{
					   //先判断是否为反义词，然后再使用词向量方法
					   double w2v_sim=getW2V_similarity(list_a[i].sName,list_b[j].sName);
					   if(w2v_sim>maxsim){
						   maxsim=w2v_sim;
						   max_aid=i;
						   max_bid=j;
					   }
					  // free(s2);
					   //delete(s2);
				   }
			   }
		   }
		   if(state==0){
			   //printf("进来了\n");
			   delta_matrix[max_aid][max_bid]=maxsim;
			   intersect[*n].Nx=list_a[max_aid].node;
			   intersect[*n].Nz=list_b[max_bid].node;
			   (*n)++;
			  // free(s1);
				//delete(s1);
			   printf("%s找到和自己相关的词了：%s ,相关度：%.4f\n",list_a[max_aid].sName,list_b[max_bid].sName,delta_matrix[max_aid][max_bid]);
		   }
	   }
	   else{
		   for(j=0;j<n_b;j++){
			   if(list_b[j].node->pre_terminal!=-1){
				   if(strcmp(list_a[i].sName,list_b[j].sName)==0){
					   intersect[*n].Nx=list_a[i].node;
					   intersect[*n].Nz=list_b[j].node;
					   (*n)++;
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=1;
				   }
				   else
					   delta_matrix[list_a[i].node->nodeID][list_b[j].node->nodeID]=0.0;
			   }
		   }
		}
   }
}
//
*/
//-------------------------------------------------------------------------------------------------------
// STRING/SEQUENCE KERNEL
//-------------------------------------------------------------------------------------------------------

double SK(TreeNode **Sx, TreeNode ** Sz, int n, int m){
  
  double *DPS =(double*) malloc((m+1)*(n+1)*sizeof(double));
  double *DP = (double*) malloc((m+1)*(n+1)*sizeof(double));
  double *kernel_mat = (double*) malloc((n+1)*sizeof(double));
	
	
  int i,j,l,p;
  double K;
  
  p = n; if (m<n) p=m;
  if(p>MAX_SUBSEQUENCE)p=MAX_SUBSEQUENCE;
	
  for (j=0; j<=m; j++)
		for (i=0; i<=n; i++) 
			DPS(i,j) = DP(i,j) =0;
	
  
	//printf("\nDPS(%d,%d)\n",n,m); fflush(stdout);
	//stampa_math(DPS,n,m); fflush(stdout);
	
  kernel_mat[0]=0;
  for (i=1; i<=n; i++)
		for (j=1; j<=m; j++)
			if(strcmp((*(Sx+i-1))->sName,(*(Sz+j-1))->sName)==0) 
			{
				DPS(i,j)=1;
				kernel_mat[0]+=DPS(i,j);
			}
			else 
				DPS(i,j)=0;
	
	
	//  printf("\n\nDPS(%d,%d)\n",n,m); fflush(stdout);
	//  stampa_math(DPS,n,m); fflush(stdout);
	//  printf("\n\nDP(%d,%d)\n",n,m);  fflush(stdout);
	//  stampa_math(DPS,n,m); fflush(stdout);
	//  printf("\n\nKernel: n=%d m=%d, %s %s \n\n",n,m,(*(Sx))->sName,(*(Sz))->sName);fflush(stdout);
	
	for(l=1;l<p;l++){
		kernel_mat[l]=0;
		for (j=0; j<=m; j++)
			DP(l-1,j)=0;
		for (i=0; i<=n; i++)
			DP(i,l-1)=0;
		
		for (i=l; i<=n; i++)
			for (j=l; j<=m; j++){
				DP(i,j) = DPS(i,j)+LAMBDA*DP(i-1,j)
				+ LAMBDA*DP(i,j-1)
				- LAMBDA2*DP(i-1,j-1);
				
				if(strcmp((*(Sx+i-1))->sName,(*(Sz+j-1))->sName)==0){
					DPS(i,j) = DP(i-1,j-1);
					kernel_mat[l] += DPS(i,j);
				}
				// else DPS[i][j] = 0;             
			}
		//      printf("\n----------------------------------\n"); printf("DPS i:%d, j:%d, l:%d\n",n,m,l+1);stampa_math(DPS,n,m);printf("DP\n");stampa_math(DP,n,m); 
	}
	//  K=kernel_mat[p-1];
	K=0;
	for(l=0;l<p;l++){
		K+=kernel_mat[l];
		//printf("String kernel of legnth %d: %1.7f \n\n",l+1,kernel_mat[l]);
	}
	
  
  free(kernel_mat);
  free(DPS);
  free(DP);
  
  return K;
}


double string_kernel_not_norm(KERNEL_PARM * kernel_parm, DOC * a, DOC * b, int i, int j){
	if(a->num_of_trees == 0 || b->num_of_trees == 0) return 0;
	else   
		if(a->num_of_trees <=i || b->num_of_trees<=j){
			printf("\nERROR: attempt to access to a not-defined item of the tree forest");
			if(a->num_of_trees <=i)
				printf("\n     : position %d of the tree forest of the first instance (line %ld)\n\n",i,a->docnum+1);
			if(b->num_of_trees<=j)
				printf("\n     : position %d of the tree forest of the second instance (line %ld)\n\n",j,b->docnum+1);
			fflush(stdout);
			exit(-1);
		}
		else if(a->forest_vec[i]->root==NULL || b->forest_vec[j]->root==NULL) return 0;
		else return SK(a->forest_vec[i]->root->pChild, b->forest_vec[j]->root->pChild, 
									 a->forest_vec[i]->root->iNoOfChildren, b->forest_vec[j]->root->iNoOfChildren);
}


double string_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b, int i, int j){
	double k;
	//printf("\ntree 1: "); writeTreeString(a->forest_vec[i]->root); printf("\ntree 2: "); 
	//writeTreeString(b->forest_vec[j]->root); printf("\n");
	//fflush(stdout);
	k=SK(a->forest_vec[i]->root->pChild, b->forest_vec[j]->root->pChild, 
			 a->forest_vec[i]->root->iNoOfChildren, b->forest_vec[j]->root->iNoOfChildren);
	//printf("STK:%f\n",k);
	
	return k;
}


//-------------------------------------------------------------------------------------------------------
// General Tree Kernel
//-------------------------------------------------------------------------------------------------------

/*
double tree_kernel_not_norm(KERNEL_PARM * kernel_parm, DOC * a, DOC * b, int i, int j){
	
  int n_pairs=0;
  int m = 0;
  nodePair intersect[MAX_NUMBER_OF_PAIRS];
	
	if(a->forest_vec[i]->orderedNodeSet != NULL && b->forest_vec[j]->orderedNodeSet != NULL) 
		determine_sub_lists(a->forest_vec[i],b->forest_vec[j],intersect,&n_pairs);
	else if(kernel_parm->first_kernel!=6 && a->forest_vec[i]->root != NULL && 
					b->forest_vec[i]->root != NULL){ 
		// if trees are not empty, from empty orderedNodeList => they are sequences and should be run wtih kernel 6
		printf("\nERROR: Tree Kernels cannot be used over sequences (positions %d or %d) \n\n",i,j);fflush(stdout);
		exit(-1);
	} 
	
	switch(kernel_parm->first_kernel) {
			
		case  -1: if(TKGENERALITY > SUBSET_TREE_KERNEL){
						printf("\nERROR: SHALLOW SEMANTIC TK kernel (-F -1) cannot be used on trees of Generality higher than 1, i.e. the subset tree kernel \n\n",i,j);
						fflush(stdout); 
						exit(-1);
				  } 
			      SIGMA = 1;
				  return evaluate_SST_ACL07(intersect,n_pairs); // SSTK kernel ACL07
			
		case  0: if(TKGENERALITY > SUBSET_TREE_KERNEL){
					printf("\nERROR: ST kernel (-F 0) cannot be used on trees of Generality higher than 1, i.e. the subset tree kernel \n\n",i,j);
					fflush(stdout); 
					exit(-1);
				  } 
			     SIGMA = 0;
			     return evaluate_SST_ST(intersect,n_pairs); // ST kernel NISP2001 (Wisnathan and Smola)
			
		case  1: if(TKGENERALITY > SUBSET_TREE_KERNEL){
					printf("\nERROR: SST kernel (-F 1) cannot be used on trees of Generality higher than 1, i.e. the subset tree kernel \n\n",i,j);
					fflush(stdout); 
					exit(-1);
				 }  
			     SIGMA = 1;
			     return evaluate_SST_ST(intersect,n_pairs); // SST kernel (Collins and Duffy, 2002)
			
		case  2: return evaluate_GSST(intersect,n_pairs); // SST kenel + bow kernel on leaves, 
																											// i.e. SST until the leaves
		case  3: REMOVE_LEAVES=0;
			///for(;m<n_pairs;m++){
			//	printf("\npair[%d] = %s + %s \n",m,intersect[m].Nx->sName,intersect[m].Nz->sName);
			//}
			     printf("\n\n进行PT计算了！\n\n");

			     return evaluate_PT(intersect,n_pairs); // PT kernel
		case  4: REMOVE_LEAVES=MU*LAMBDA2;
			     return evaluate_PT(intersect,n_pairs); // PT kernel no leaves
		case  6: return string_kernel(kernel_parm, a, b, i, j);
			
		default: printf("\nERROR: Tree Kernel number %ld not available \n\n",kernel_parm->first_kernel);
			     fflush(stdout);
			     exit(-1);
	}
	
	return 0;
}*/

double tree_kernel_not_norm(KERNEL_PARM * kernel_parm, FOREST * a, FOREST * b){
	
  int n_pairs=0;
  int m = 0;
  int i=0;
  nodePair intersect[MAX_NUMBER_OF_PAIRS];
	
	if(a->orderedNodeSet != NULL && b->orderedNodeSet != NULL) {
		determine_sub_lists(a,b,intersect,&n_pairs);
	}
	else if(kernel_parm->first_kernel!=6 && a->root != NULL && 
					b->root != NULL){ 
		// if trees are not empty, from empty orderedNodeList => they are sequences and should be run wtih kernel 6
		printf("\nERROR: Tree Kernels cannot be used over sequences  \n\n");
		fflush(stdout);
		exit(-1);
	} 
	/*
	 printf("\n\nNormal PRINTing\n");
    for(i=0;i<n_pairs;i++)
		printf("\npairs[%d] :%s  %s ,ID: Nx = %d, Nz = %d,delta = %f",i+1,intersect[i].Nx->sName, intersect[i].Nz->sName,intersect[i].Nx->nodeID,intersect[i].Nz->nodeID,delta_matrix[intersect[i].Nx->nodeID][intersect[i].Nz->nodeID]);
    
     printf("\n\nOrdered PRINTing\n");
	 */
	//printf("\n kernel_parm->first_kernel = %d a->root = %s;b->root = %s \n",kernel_parm->first_kernel,a->root->sName,b->root->sName);
	switch(kernel_parm->first_kernel) {
			
		case  -1: if(TKGENERALITY > SUBSET_TREE_KERNEL){
						printf("\nERROR: SHALLOW SEMANTIC TK kernel cannot be used on trees of Generality higher than 1, i.e. the subset tree kernel \n\n");
						fflush(stdout); 
						exit(-1);
				  } 
			      SIGMA = 1;
				  return evaluate_SST_ACL07(intersect,n_pairs); // SSTK kernel ACL07
			
		case  0: if(TKGENERALITY > SUBSET_TREE_KERNEL){
					printf("\nERROR: ST kernel  cannot be used on trees of Generality higher than 1, i.e. the subset tree kernel \n\n");
					fflush(stdout); 
					exit(-1);
				  } 
			     SIGMA = 0;
			     return evaluate_SST_ST(intersect,n_pairs); // ST kernel NISP2001 (Wisnathan and Smola)
			
		case  1: if(TKGENERALITY > SUBSET_TREE_KERNEL){
					printf("\nERROR: SST kernel cannot be used on trees of Generality higher than 1, i.e. the subset tree kernel \n\n");
					fflush(stdout); 
					exit(-1);
				 }  
			     SIGMA = 1;
			     return evaluate_SST_ST(intersect,n_pairs); // SST kernel (Collins and Duffy, 2002)
			
		case  2: return evaluate_GSST(intersect,n_pairs); // SST kenel + bow kernel on leaves, 
																											// i.e. SST until the leaves
		case  3: //REMOVE_LEAVES=0;
			     return evaluate_PT(intersect,n_pairs,0); // PT kernel

		case  4: //REMOVE_LEAVES=MU*LAMBDA2;
			//{
				//printf("进来了\n");
			     return evaluate_PT(intersect,n_pairs,1); // PT kernel no leaves
			//}
	//	case  6: return string_kernel(kernel_parm, a, b, i, j);
			
		default: printf("\nERROR: Tree Kernel number %ld not available \n\n",kernel_parm->first_kernel);
			     fflush(stdout);
			     exit(-1);
	}
	
	return 0;
}

/*
*STK
*/
/*double tree_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b, int i, int j){

  int n_pairs=0;

  double k=0;
  
  nodePair intersect[MAX_NUMBER_OF_PAIRS];
   
  if (b->num_of_trees > j && a->num_of_trees > i){
  
         determine_sub_lists(a->forest_vec[i],b->forest_vec[j],intersect,&n_pairs);
         k =(evaluateParseTreeKernel(intersect,n_pairs));
         //printf("\n\n(i,j)=(%d,%d)= Kernel :%f norm1,norm2 (%f,%f)\n",i,j,k,a->forest_vec[i]->twonorm_PT, b->forest_vec[j]->twonorm_PT);
  } 
  else{
       printf("\nERROR: attempting to access to a tree not defined in the data\n\n");
       exit(-1);
  }
  
  return k;
}*/

/*
double tree_kernel(KERNEL_PARM *kernel_parm, DOC *a, DOC *b, int i, int j){
	if(a->num_of_trees == 0 || b->num_of_trees == 0) return 0;
	else   
		if(a->num_of_trees <=i || b->num_of_trees<=j){
			printf("\nERROR: attempt to access to a not-defined item of the tree forest");
			if(a->num_of_trees <=i)
				printf("\n     : position %d of the tree forest of the first instance (line %ld)\n\n",i,a->docnum+1);
			if(b->num_of_trees<=j)
				printf("\n     : position %d of the tree forest of the second instance (line %ld)\n\n",j,b->docnum+1);
			fflush(stdout);
			exit(-1);
		}       
		else if(a->forest_vec[i]->root==NULL || b->forest_vec[j]->root==NULL)
			return 0;
		else {
			double   k= tree_kernel_not_norm(kernel_parm, a, b, i, j);
			double m = k/sqrt(a->forest_vec[i]->twonorm_PT * b->forest_vec[j]->twonorm_PT);
			printf("\nkernel %f, = %f ,= %f \n",k,a->forest_vec[i]->twonorm_PT,b->forest_vec[j]->twonorm_PT);
			return k;
		}
	
	
}*/
double tree_kernel(KERNEL_PARM *kernel_parm, FOREST *a, FOREST *b){
	if(a->root==NULL||b->root==NULL)
		return 0;
	else {
		double   k= tree_kernel_not_norm(kernel_parm, a, b);
		double m = k/sqrt(a->twonorm_PT * b->twonorm_PT);
		//printf("\nkernel %f, = %f ,= %f \n",k,a->forest_vec[i]->twonorm_PT,b->forest_vec[j]->twonorm_PT);
		return m;
	}
}

/*-----------------------------------------------------------------------------------------------------*/

double basic_kernel_not_norm(KERNEL_PARM *kernel_parm, DOC *a, DOC *b, int i, int j) 
/* calculate the kernel function */
{
  switch(kernel_parm->second_kernel) {
    case 0: /* linear */ 
			return 0;//sprod_ss(a->vectors[i]->words,b->vectors[j]->words); 
    case 1: /* polynomial */
			return (double) pow(((double)kernel_parm->coef_lin)*(double)1.0//sprod_ss(a->vectors[i]->words,b->vectors[j]->words)
													+(double)kernel_parm->coef_const,(double) kernel_parm->poly_degree);
    case 2: /* radial basis function */
			return 0;//sprod_ss(a->vectors[i]->words,b->vectors[j]->words)+b->vectors[i]->twonorm_sq)));
    case 3: /* sigmoid neural net */
			return 0;//sprod_ss(a->vectors[i]->words,b->vectors[j]->words)+kernel_parm->coef_const)); 
    case 4: /* custom-kernel supplied in file kernel.h*/
			return(custom_kernel(kernel_parm,a,b));
			/* string kernel*/
    case 6: return string_kernel(kernel_parm,a,b,i,j);
    default: printf("Error: The kernel function to be combined with the Tree Kernel is unknown\n"); 
			fflush(stdout);
			exit(1);
  }
}

     /* calculate the kernel function */
/*double basic_kernel(KERNEL_PARM *kernel_parm, DOC *a, DOC *b, int i, int j) 

     
{
	 
  switch(kernel_parm->second_kernel) {
	   printf("\n\nkernel_parm->second_kernel\n\n",kernel_parm->second_kernel);
    case 0: // linear 
            return(sprod_i(a, b, i, j)); 
    case 1: // polynomial 
            return (double) pow(((double)kernel_parm->coef_lin)*(double)sprod_i(a, b, i, j)
                   +(double)kernel_parm->coef_const,(double) kernel_parm->poly_degree);
    case 2: //radial basis function 
            return(exp(-kernel_parm->rbf_gamma*(a->vectors[i]->twonorm_sq-2*sprod_i(a, b, i, j)+b->vectors[i]->twonorm_sq)));
    case 3: // sigmoid neural net 
            return(tanh(kernel_parm->coef_lin*sprod_i(a, b, i, j)+kernel_parm->coef_const)); 
    case 4: // custom-kernel supplied in file kernel.h
		  //  printf("\n\n运行用户定义的核了！！\n\n");
            return(custom_kernel(kernel_parm,a,b));
    case 5: // TREE KERNEL 
            return(tree_kernel(kernel_parm,a,b,i,j));
            	     
    default: printf("Error: The kernel function to be combined with the Tree Kernel is unknown\n"); exit(1);
  }
}
*/

double basic_kernel(KERNEL_PARM *kernel_parm, DOC *a, DOC *b, int i, int j){
	if(a->num_of_vectors == 0 || b->num_of_vectors == 0) return 0;
	else   
		if(a->num_of_vectors <=i || b->num_of_vectors <= j){
			printf("\nERROR: attempt to access to a not-defined item of the vector set");
			if(a->num_of_vectors <=i)
				printf("\n     :position %d of the vector set of the first istance (line %ld)\n\n",i,a->docnum+1);
			if(a->num_of_vectors <=j)
				printf("\n     :position %d of the vector set of the second instance (line %ld)\n\n",j,b->docnum+1);
			fflush(stdout);
			exit(-1);
		}
		else if(a->vectors[i]==NULL || b->vectors[j]==NULL) return 0;
		else return basic_kernel_not_norm(kernel_parm, a, b,  i, j)/
			sqrt(a->vectors[i]->twonorm_STD*b->vectors[j]->twonorm_STD);
}


/*-----------------------------------------------------------------------------------------------------*/

/*void evaluateNorma(KERNEL_PARM * kernel_parm, DOC * d){

    int n_pairs=0,
		te,
        i;
        
    double k=0;
    nodePair intersect[MAX_NUMBER_OF_PAIRS];

//printf("doc ID :%d \n",d->docnum);
//printf("num of vectors:%d \n",d->num_of_vectors);
//fflush(stdout);
		
      for(i=0;i < d->num_of_trees;i++){
      //  TESTS 

        //  printf("\n\n\nnode ID: %d \n", d->forest_vec[i]->root->nodeID); fflush(stdout);

        //  printf("node list length: %d\n", d->forest_vec[i]->listSize);

        //  printf("doc ID :%d \n",d->docnum);

         // printf("tree: <"); writeTreeString(d->forest_vec[i]->root);printf(">");
    
          //printf("\n\n"); fflush(stdout);
          //
          

          // this avoids to check for norm == 0
          //printf ("Norm %f\n",k1);

          determine_sub_lists(d->forest_vec[i],d->forest_vec[i],intersect,&n_pairs);

          k =(evaluateParseTreeKernel(intersect,n_pairs));
      //     k = basic_kernel(kernel_parm, d, d, i, i);

          if(k!=0 && (kernel_parm->normalization == 1 || kernel_parm->normalization == 3)) 
               d->forest_vec[i]->twonorm_PT=k; 
          else d->forest_vec[i]->twonorm_PT=1; 
      }

  // SECOND KERNEL NORM EVALUATION 
	  
    
      for(i=0;i < d->num_of_vectors;i++){
//      for(i=0;i < 61 && i<d->num_of_vectors;i+=60){
        
          d->vectors[i]->twonorm_STD=1; // basic-kernel normalizes the standard kernels
                                        // this also avoids to check for norm == 0
                                        
          k = basic_kernel(kernel_parm, d, d, i, i);
  // printf("\n\njinlaile  k= %f !!\n\n",k);                
          if(k!=0 && (kernel_parm->normalization == 2 || kernel_parm->normalization == 3))
               d->vectors[i]->twonorm_STD=k; // if selected normalization is applied
               
          d->vectors[i]->twonorm_sq=sprod_ss(d->vectors[i]->words,d->vectors[i]->words);
       } 
      // maintain the compatibility with svm-light single linear vector 
        if(d->num_of_vectors>0) d->twonorm_sq=sprod_ss(d->vectors[0]->words,d->vectors[0]->words);
        else d->twonorm_sq=0;
 }*/

void evaluateNorma(KERNEL_PARM * kernel_parm, DOC * d){
	int i;
	double k=0;
	double ks[2]={0,0};
	//double testtest;
	//printf("doc ID :%d \n",d->docnum);
	//printf("num of vectors:%d \n",d->num_of_vectors);
	//fflush(stdout);
	
	short kernel_type_tmp=kernel_parm->first_kernel; //save parameters from command line
	double lambda_tmp=LAMBDA; 
	double mu_tmp=MU;
	short TKG_tmp=TKGENERALITY;
	if(d->num_of_trees>=2){  //ensure pair of tree
		if(tree_kernel_params[0].kernel_type==END_OF_TREE_KERNELS)
					PARAM_VECT=0;
		if(PARAM_VECT==1){
					if (tree_kernel_params[0].kernel_type!=NOKERNEL){
						//printf("\n!=NOKERNEL!\n");
						if(tree_kernel_params[0].normalization==1){
							kernel_parm->first_kernel=tree_kernel_params[0].kernel_type;
							LAMBDA = tree_kernel_params[0].lambda; 
							LAMBDA2 = LAMBDA*LAMBDA;
							MU=tree_kernel_params[0].mu;
							TKGENERALITY=tree_kernel_params[0].TKGENERALITY;
							//printf("%f  %f",LAMBDA,MU);
							d->forest_vec[0]->twonorm_PT = tree_kernel_params[0].weight*tree_kernel_not_norm(kernel_parm, d->forest_vec[0], d->forest_vec[0]);
							d->forest_vec[1]->twonorm_PT = tree_kernel_params[1].weight*tree_kernel_not_norm(kernel_parm, d->forest_vec[1],d->forest_vec[1]);
							if(d->forest_vec[0]->twonorm_PT==0) 
								d->forest_vec[0]->twonorm_PT=1;
							if(d->forest_vec[1]->twonorm_PT==1) 
								d->forest_vec[1]->twonorm_PT=1;
						}
						else {
							d->forest_vec[0]->twonorm_PT=1;
							d->forest_vec[1]->twonorm_PT=1;
						}
					}
				}
				else{
					//printf("\nPARAM_VECT!=1  , PARAM_VECT = %d \n",PARAM_VECT);
					
					k = tree_kernel_not_norm(kernel_parm, d->forest_vec[0], d->forest_vec[1]);
					ks[0] = tree_kernel_not_norm(kernel_parm, d->forest_vec[0], d->forest_vec[0]);
					ks[1] = tree_kernel_not_norm(kernel_parm, d->forest_vec[1], d->forest_vec[1]);
					//printf("\nk=%.4f\n",k);
					//testtest= getW2V_similarity("bathing","sing");
					if(k!=0 && (kernel_parm->normalization == 1 || kernel_parm->normalization == 3)) {
						d->forest_vec[0]->twonorm_PT=ks[0]; 
						d->forest_vec[1]->twonorm_PT=ks[1]; 

						k = k / sqrt(d->forest_vec[0]->twonorm_PT * d->forest_vec[1]->twonorm_PT);
					}
					else{
						d->forest_vec[0]->twonorm_PT=1; 
						d->forest_vec[1]->twonorm_PT=1; 
					}
					d->k_score = k;
					//printf("\nk = % .4f\n",d->k_score);
				}
			
				// this avoids to check for norm == 0
				// printf ("Norm %f\n",k);fflush(stdout);
			
			kernel_parm->first_kernel=kernel_type_tmp; //re-set command line parameters 
			LAMBDA=lambda_tmp; 
			MU=mu_tmp;
			TKGENERALITY=TKG_tmp;
     }
  
  
  /* SECOND KERNEL NORM EVALUATION */
	/*
	for(i=0;i < d->num_of_vectors;i++){
		
		if(d->num_of_vectors>0 && d->vectors[i]!=NULL){
			d->vectors[i]->twonorm_STD=1; // basic-kernel normalizes the standard kernels
																		// this also avoids to check for norm == 0
			k = basic_kernel_not_norm(kernel_parm, d, d, i, i);               
			d->vectors[i]->twonorm_sq=sprod_ss(d->vectors[i]->words,d->vectors[i]->words);
			if(k!=0 && (kernel_parm->normalization == 2 || kernel_parm->normalization == 3))
				d->vectors[i]->twonorm_STD=k; // if selected normalization is applied
																			//            printf ("Norm %f\n",k);
			
		}
		
	} 
	// maintain the compatibility with svm-light single linear vector 
	if(d->num_of_vectors>0 && d->vectors[0]!=NULL) 
		d->twonorm_sq=sprod_ss(d->vectors[0]->words,d->vectors[0]->words);
	else 
		d->twonorm_sq=0;
	
	*/
}
/*
void evaluateNorma(KERNEL_PARM * kernel_parm, DOC * d){
	
	int  i;
	double k=0;
	
	//printf("doc ID :%d \n",d->docnum);
	//printf("num of vectors:%d \n",d->num_of_vectors);
	//fflush(stdout);
	
	short kernel_type_tmp=kernel_parm->first_kernel; //save parameters from command line
	double lambda_tmp=LAMBDA; 
	double mu_tmp=MU;
	short TKG_tmp=TKGENERALITY;
	
	for(i=0;i < d->num_of_trees;i++){
		
		k=0;
		
		if (d->num_of_trees > i && d->forest_vec[i]->root!=NULL){
			
			// TESTS
			 
			 //printf("\n\n\ndoc ID :%ld \n",d->docnum);fflush(stdout);
			 
			// printf("node ID: %d \n", d->forest_vec[i]->root->nodeID); fflush(stdout);
			 
			// printf("node list length: %d\n", d->forest_vec[i]->listSize);fflush(stdout);
			// 
			// printf("tree: <"); writeTreeString(d->forest_vec[i]->root);printf(">");fflush(stdout);
			 
			// printf("\n\n"); fflush(stdout);
			 //
			
			if(tree_kernel_params[i].kernel_type==END_OF_TREE_KERNELS)
				PARAM_VECT=0;
			
			if(PARAM_VECT==1){
				if (tree_kernel_params[i].kernel_type!=NOKERNEL){
					if(tree_kernel_params[i].normalization==1){
						kernel_parm->first_kernel=tree_kernel_params[i].kernel_type;
						LAMBDA = tree_kernel_params[i].lambda; 
						LAMBDA2 = LAMBDA*LAMBDA;
						MU=tree_kernel_params[i].mu;
						TKGENERALITY=tree_kernel_params[i].TKGENERALITY;
						d->forest_vec[i]->twonorm_PT = tree_kernel_params[i].weight*tree_kernel_not_norm(kernel_parm, d, d, i, i);
						if(d->forest_vec[i]->twonorm_PT==0) 
							d->forest_vec[i]->twonorm_PT=1;
					}
					else d->forest_vec[i]->twonorm_PT=1;
				}
			}
			else{
				k = tree_kernel_not_norm(kernel_parm, d, d, i, i);
				printf("\nk = % f\n",k);
				if(k!=0 && (kernel_parm->normalization == 1 || kernel_parm->normalization == 3)) 
					d->forest_vec[i]->twonorm_PT=k; 
				else
					d->forest_vec[i]->twonorm_PT=1; 
			}
			
			// this avoids to check for norm == 0
			// printf ("Norm %f\n",k);fflush(stdout);
			
		}
		kernel_parm->first_kernel=kernel_type_tmp; //re-set command line parameters 
		LAMBDA=lambda_tmp; 
		MU=mu_tmp;
		TKGENERALITY=TKG_tmp;
	}
  
  
  // SECOND KERNEL NORM EVALUATION 
	
	for(i=0;i < d->num_of_vectors;i++){
		
		if(d->num_of_vectors>0 && d->vectors[i]!=NULL){
			d->vectors[i]->twonorm_STD=1; // basic-kernel normalizes the standard kernels
																		// this also avoids to check for norm == 0
			k = basic_kernel_not_norm(kernel_parm, d, d, i, i);               
			d->vectors[i]->twonorm_sq=sprod_ss(d->vectors[i]->words,d->vectors[i]->words);
			if(k!=0 && (kernel_parm->normalization == 2 || kernel_parm->normalization == 3))
				d->vectors[i]->twonorm_STD=k; // if selected normalization is applied
																			//            printf ("Norm %f\n",k);
			
		}
		
	} 
	// maintain the compatibility with svm-light single linear vector 
	if(d->num_of_vectors>0 && d->vectors[0]!=NULL) 
		d->twonorm_sq=sprod_ss(d->vectors[0]->words,d->vectors[0]->words);
	else 
		d->twonorm_sq=0;
	
}
 */
/***************************************************************************************/
/*                           KERNELS OVER SET OF KERNELS                               */
/***************************************************************************************/



// sequence summation of trees

double sequence_tree_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){//all_vs_all_tree_kernel

  int i;
  double k;
  int n_pairs=0;
  nodePair intersect[MAX_NUMBER_OF_PAIRS];
  
//   printf("\n\nDocum %ld and %ld, size=(%d,%d)\n",a->docnum,b->docnum,a->num_of_trees,b->num_of_trees);
   k=0;

   for(i=0; i< a->num_of_trees && i< b->num_of_trees; i++){
        /* 
         printf("\n\n\n nodes: %d  %d\n", a->forest_vec[i]->root->nodeID,b->forest_vec[i]->root->nodeID);
         printf("node list lenghts: %d  %d\n", a->forest_vec[i]->listSize,b->forest_vec[i]->listSize);
         */
   //      printf("\ntree 1: <"); writeTreeString(a->forest_vec[i]->root);
  
    //     printf(">\ntree 2: <"); writeTreeString(b->forest_vec[i]->root);printf(">\n"); fflush(stdout);

      if(a->forest_vec[i]!=NULL && b->forest_vec[i]!=NULL){                                
         determine_sub_lists(a->forest_vec[i],b->forest_vec[i],intersect,&n_pairs);
         k+= (evaluateParseTreeKernel(intersect,n_pairs)/
         sqrt(a->forest_vec[i]->twonorm_PT * b->forest_vec[i]->twonorm_PT));
      }

 //        printf("\n\n(i,i)=(%d,%d)= Kernel-Sequence :%f norm1,norm2 (%f,%f)\n",i,i,k,a->forest_vec[i]->twonorm_PT, b->forest_vec[i]->twonorm_PT);
  }
   return k;
}


// all vs all summation of trees


double AVA_tree_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){//all_vs_all_tree_kernel

  int i,
      j,
      n_pairs=0;

  double k=0;
  nodePair intersect[MAX_NUMBER_OF_PAIRS];
   
   //printf("first elements: %s  %s\n", a->orderedNodeSet->sName,b->orderedNodeSet->sName);
   //printf("\n\n---------------------------------------------------------\n\n");fflush(stdout);
   //printf("doc IDs :%d %d",a->docnum,b->docnum);

   if (b->num_of_trees == 0 || a->num_of_trees==0) return 0;
   
   for(i=0; i< a->num_of_trees; i++)
      for(j=0; j< b->num_of_trees;j++){
  
      if(a->forest_vec[i]!=NULL && b->forest_vec[j]!=NULL){

        /* 
         printf("\n\n\n nodes: %d  %d\n", a->forest_vec[i]->root->nodeID,b->forest_vec[i]->root->nodeID);
         printf("node list lenghts: %d  %d\n", a->forest_vec[i]->listSize,b->forest_vec[i]->listSize);
         */
        //printf("\ntree 1: <"); writeTreeString(a->forest_vec[i]->root);
  
        //printf(">\ntree 2: "); writeTreeString(b->forest_vec[j]->root);printf(">\n");
        //fflush(stdout);

         determine_sub_lists(a->forest_vec[i],b->forest_vec[j],intersect,&n_pairs);

         k+= (evaluateParseTreeKernel(intersect,n_pairs)/
         sqrt(a->forest_vec[i]->twonorm_PT * b->forest_vec[j]->twonorm_PT));
         //printf("\n\n(i,j)=(%d,%d)= Kernel :%f norm1,norm2 (%f,%f)\n",i,j,k,a->forest_vec[i]->twonorm_PT, b->forest_vec[j]->twonorm_PT);
      }
     }
    //printf("\n---------------------------------------------------------------\n"); fflush(stdout);
  

   return k;
}



// sequence summation of vectors


double sequence(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){

  int i;
  double k=0;
    
   for(i=0; i< a->num_of_vectors && i< b->num_of_vectors; i++)
     
      if(a->vectors[i]!=NULL && b->vectors[i]!=NULL){
         k+= basic_kernel(kernel_parm, a, b, i, i)/
         sqrt(a->vectors[i]->twonorm_STD * b->vectors[i]->twonorm_STD);
      }
   return k;
}


// all vs all summation of vectors

double AVA(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){

  int i,
      j;
      
  double k=0;
    
   for(i=0; i< a->num_of_vectors; i++)
      for(j=0; j< b->num_of_vectors;j++){
  
         if(a->vectors[i]!=NULL && b->vectors[j]!=NULL){
             k+= basic_kernel(kernel_parm, a, b, i, j)/
             sqrt(a->vectors[i]->twonorm_STD * b->vectors[j]->twonorm_STD);
         }
      }
  return k;
}


// kernel for entailments [Zanzotto and Moschitti, ACL 2005]

double ACL2006_Entailment_kernel(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){

  int i,
      n_pairs=0;

  double k=0,
         max=0;

  nodePair intersect[MAX_NUMBER_OF_PAIRS];
   
//   LAMBDA = kernel_parm->lambda; //faster access for lambda
//   SIGMA = kernel_parm->sigma;


   // printf("first elements: %s  %s\n", a->orderedNodeSet->sName,b->orderedNodeSet->sName);
   // printf("\n\n---------------------------------------------------------\n\n");fflush(stdout);
   // printf("doc IDs :%d %d",a->docnum,b->docnum);

 if (b->num_of_trees > 1 && a->num_of_trees>1){
   //kk=0;
   //for(i=0; i< 2; i++)
   //  for(j=0; j< 2;j++){
   //   if(a->forest_vec[i]!=NULL && b->forest_vec[j]!=NULL){
        /* 
         printf("\n\n\n nodes: %d  %d\n", a->forest_vec[i]->root->nodeID,b->forest_vec[i]->root->nodeID);
         printf("node list lenghts: %d  %d\n", a->forest_vec[i]->listSize,b->forest_vec[i]->listSize);
         */
        //printf("\ntree 1: <"); writeTreeString(a->forest_vec[i]->root);
        //printf(">\ntree 2: "); writeTreeString(b->forest_vec[j]->root);printf(">\n");
        //fflush(stdout);
         //determine_sub_lists(a->forest_vec[i],b->forest_vec[j],intersect,&n_pairs);
         //kk+= (evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[i]->twonorm_PT * b->forest_vec[j]->twonorm_PT));
         //printf("\n\n(i,j)=(%d,%d)= Kernel :%f norm1,norm2 (%f,%f)\n",i,j,k,a->forest_vec[i]->twonorm_PT, b->forest_vec[j]->twonorm_PT);
    //  }
    // }
        // determine_sub_lists(a->forest_vec[0],b->forest_vec[0],intersect,&n_pairs);
        // kk+= (evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[0]->twonorm_PT * b->forest_vec[0]->twonorm_PT));
        // determine_sub_lists(a->forest_vec[1],b->forest_vec[1],intersect,&n_pairs);
        // kk+= (evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[1]->twonorm_PT * b->forest_vec[1]->twonorm_PT));
   }
   //kk = 0; FMZ added to test max contribution
   
   max = 0;
   if (b->num_of_trees > 2  && a->num_of_trees > 2){ 
   if (b->num_of_trees >  a->num_of_trees) {  
     for(i=2;i<b->num_of_trees ;i+=2){
      if(a->forest_vec[2]!=NULL && b->forest_vec[i]!=NULL){
          //  printf("\n\n\n nodes: %d  %d\n", a->forest_vec[i]->root->nodeID,b->forest_vec[i]->root->nodeID);
        //  printf("node list lenghts: %d  %d\n", a->forest_vec[i]->listSize,b->forest_vec[i]->listSize);
          // printf("\ntree 1: "); writeTreeString(a->forest_vec[i]->root);
          //printf("\ntree 2: "); writeTreeString(b->forest_vec[i]->root);
         // fflush(stdout);
         determine_sub_lists(a->forest_vec[2],b->forest_vec[i],intersect,&n_pairs);
	 k= evaluateParseTreeKernel(intersect,n_pairs)/
         sqrt(a->forest_vec[2]->twonorm_PT * b->forest_vec[i]->twonorm_PT);
         determine_sub_lists(a->forest_vec[3],b->forest_vec[i+1],intersect,&n_pairs);
         k+= evaluateParseTreeKernel(intersect,n_pairs)/
         sqrt(a->forest_vec[3]->twonorm_PT * b->forest_vec[i+1]->twonorm_PT);
         if(max<k)max=k;
         //printf("\n\nKernel :%f \n",k);
       } 
      }
     } else {
	     
      for(i=2;i<a->num_of_trees ;i+=2){
       if(a->forest_vec[i]!=NULL && b->forest_vec[2]!=NULL){
          //  printf("\n\n\n nodes: %d  %d\n", a->forest_vec[i]->root->nodeID,b->forest_vec[i]->root->nodeID);
        //  printf("node list lenghts: %d  %d\n", a->forest_vec[i]->listSize,b->forest_vec[i]->listSize);
          // printf("\ntree 1: "); writeTreeString(a->forest_vec[i]->root);
          //printf("\ntree 2: "); writeTreeString(b->forest_vec[i]->root);
         // fflush(stdout);
         determine_sub_lists(a->forest_vec[i],b->forest_vec[2],intersect,&n_pairs);
	     k= evaluateParseTreeKernel(intersect,n_pairs)/
         sqrt(a->forest_vec[i]->twonorm_PT * b->forest_vec[2]->twonorm_PT);
         determine_sub_lists(a->forest_vec[i+1],b->forest_vec[3],intersect,&n_pairs);
         k+= evaluateParseTreeKernel(intersect,n_pairs)/
         sqrt(a->forest_vec[i+1]->twonorm_PT * b->forest_vec[3]->twonorm_PT);
         if(max<k)max=k;
         //printf("\n\nKernel :%f \n",k);
       }
      }
     }
    }
   //printf("\n---------------------------------------------------------------\n");fflush(stdout);
//printf("\n\nKernel :%f \n",max);

  if(kernel_parm->combination_type=='+' && (a->vectors!=NULL && b->vectors!=NULL))
       return basic_kernel(kernel_parm, a, b, 0, 0)/
              sqrt(a->vectors[0]->twonorm_STD * b->vectors[0]->twonorm_STD)+ 
              kernel_parm->tree_constant*max;
  else return max;
}


// Kernel for re-ranking predicate argument structures, [Moschitti, CoNLL 2006]

double SRL_re_ranking_CoNLL2006(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){//all_vs_all_tree_kernel

  int n_pairs=0;

  double k1=0,k2=0;
  nodePair intersect[MAX_NUMBER_OF_PAIRS];
   
   //printf("first elements: %s  %s\n", a->orderedNodeSet->sName,b->orderedNodeSet->sName);
   //printf("\n\n---------------------------------------------------------\n\n");fflush(stdout);
  
 if(kernel_parm->kernel_type==11 || kernel_parm->kernel_type==12){

   if(a->num_of_trees!=0 && b->num_of_trees!= 0){

      if(a->forest_vec[0]==NULL || a->forest_vec[1]==NULL
         || b->forest_vec[0]==NULL || b->forest_vec[1]==NULL){
         printf("ERROR: two trees for each instance are needed");
         exit(-1);
      }

         determine_sub_lists(a->forest_vec[0],b->forest_vec[0],intersect,&n_pairs);
         k1+= evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[0]->twonorm_PT * b->forest_vec[0]->twonorm_PT);

         determine_sub_lists(a->forest_vec[1],b->forest_vec[1],intersect,&n_pairs);
         k1+= evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[1]->twonorm_PT * b->forest_vec[1]->twonorm_PT);

         determine_sub_lists(a->forest_vec[1],b->forest_vec[0],intersect,&n_pairs);
         k1-= evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[1]->twonorm_PT * b->forest_vec[0]->twonorm_PT);

         determine_sub_lists(a->forest_vec[0],b->forest_vec[1],intersect,&n_pairs);
         k1-= evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[0]->twonorm_PT * b->forest_vec[1]->twonorm_PT);

         if(kernel_parm->kernel_type==12)k1*=kernel_parm->tree_constant;

    //printf("\n\n(i,j)=(%d,%d)= Kernel :%f norm1,norm2 (%f,%f)\n",i,j,k,a->forest_vec[i]->twonorm_PT, b->forest_vec[j]->twonorm_PT);
    //printf("\n---------------------------------------------------------------\n"); fflush(stdout);
   }
}

// to use all the local argument classifier features as in [Toutanova ACL2005]

if(kernel_parm->kernel_type==13 || kernel_parm->kernel_type==12){
          k2+=sequence_ranking(kernel_parm, a, b, 0, 0);
          k2+=sequence_ranking(kernel_parm, a, b, 1, 1);
          k2-=sequence_ranking(kernel_parm, a, b, 0, 1);
          k2-=sequence_ranking(kernel_parm, a, b, 1, 0);
          }


// use only 1 vector for each predicate argument structure
          
//     if(kernel_parm->kernel_type==13 || kernel_parm->kernel_type==12){
//       if(a->num_of_vectors>0 && b->num_of_vectors>0){
//         k2+=basic_kernel(kernel_parm, a, b, 0, 0)/
//          sqrt(a->vectors[0]->twonorm_STD * b->vectors[0]->twonorm_STD);
//	     k2+=basic_kernel(kernel_parm, a, b, 1, 1)/
//          sqrt(a->vectors[1]->twonorm_STD * b->vectors[1]->twonorm_STD);
//         k2-=basic_kernel(kernel_parm, a, b, 0, 1)/
//          sqrt(a->vectors[0]->twonorm_STD * b->vectors[1]->twonorm_STD);
//         k2-=basic_kernel(kernel_parm, a, b, 1, 0)/
//          sqrt(a->vectors[1]->twonorm_STD * b->vectors[0]->twonorm_STD);
//          }
//       }  

    
   // printf("kernel: %f\n",k2);

   return  k1+k2;//(k1*k2 + k2 + (k1+1)*(k1+1));
}


// ranking algorithm based on only trees. It can be used for parse-tree re-ranking

double tree_kernel_ranking(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){//all_vs_all_tree_kernel

  int n_pairs=0;

  double k=0;
  nodePair intersect[MAX_NUMBER_OF_PAIRS];
   
      //printf("first elements: %s  %s\n", a->orderedNodeSet->sName,b->orderedNodeSet->sName);
   //printf("\n\n---------------------------------------------------------\n\n");fflush(stdout);
  
   if(a->num_of_trees==0 || b->num_of_trees== 0) return 0;
   
      if(a->forest_vec[0]==NULL || a->forest_vec[1]==NULL
         || b->forest_vec[0]==NULL || b->forest_vec[1]==NULL){
         printf("ERROR: two trees for each instance are needed");
         exit(-1);
      }

         determine_sub_lists(a->forest_vec[0],b->forest_vec[0],intersect,&n_pairs);
         k+= evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[0]->twonorm_PT * b->forest_vec[0]->twonorm_PT);

         determine_sub_lists(a->forest_vec[1],b->forest_vec[1],intersect,&n_pairs);
         k+= evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[1]->twonorm_PT * b->forest_vec[1]->twonorm_PT);

         determine_sub_lists(a->forest_vec[1],b->forest_vec[0],intersect,&n_pairs);
         k-= evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[1]->twonorm_PT * b->forest_vec[0]->twonorm_PT);

         determine_sub_lists(a->forest_vec[0],b->forest_vec[1],intersect,&n_pairs);
         k-= evaluateParseTreeKernel(intersect,n_pairs)/sqrt(a->forest_vec[0]->twonorm_PT * b->forest_vec[1]->twonorm_PT);

         //printf("\n\n(i,j)=(%d,%d)= Kernel :%f norm1,norm2 (%f,%f)\n",i,j,k,a->forest_vec[i]->twonorm_PT, b->forest_vec[j]->twonorm_PT);
    //printf("\n---------------------------------------------------------------\n"); fflush(stdout);
           
   return k;
}


// ranking algorithm based on only vectors. For example, it can be used for ranking documents wrt a query


double vector_ranking(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){

  double k=0;
 
         if(a->num_of_vectors==0 || b->num_of_vectors==0) return 0;
 
          k+=basic_kernel(kernel_parm, a, b, 0, 0);
          sqrt(a->vectors[0]->twonorm_STD * b->vectors[0]->twonorm_STD);
          k+=basic_kernel(kernel_parm, a, b, 1, 1)/
          sqrt(a->vectors[1]->twonorm_STD * b->vectors[1]->twonorm_STD);
          k-=basic_kernel(kernel_parm, a, b, 0, 1)/
          sqrt(a->vectors[0]->twonorm_STD * b->vectors[1]->twonorm_STD);
          k-=basic_kernel(kernel_parm, a, b, 1, 0)/
          sqrt(a->vectors[1]->twonorm_STD * b->vectors[0]->twonorm_STD);
   return k;
}


// ranking algorithm based on tree forests. In this case the ranked objetcs are described by a forest

double vector_sequence_ranking(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){

  double k=0;

  k+=sequence_ranking(kernel_parm, a, b, 0, 0); // ranking with sequences of vectors
  k+=sequence_ranking(kernel_parm, a, b, 1, 1);
  k-=sequence_ranking(kernel_parm, a, b, 0, 1);
  k-=sequence_ranking(kernel_parm, a, b, 1, 0);
  
  return k;
}


/* uses all the vectors in the vector set for ranking */
/* this means that there are n/2 vectors for the first pair and n/2 for the second pair */

double sequence_ranking(KERNEL_PARM * kernel_parm, DOC * a, DOC * b, int memberA, int memberB){//all_vs_all vectorial kernel

  int i;
  int startA, startB;
  
  double k=0;
  
   startA= a->num_of_vectors*memberA/2;
   startB= b->num_of_vectors*memberB/2;
   
   if(a->num_of_vectors==0 || b->num_of_vectors==0) return 0;
   
//   for(i=0; i< a->num_of_vectors/2 && i< b->num_of_vectors/2; i++)
  for(i=0; i<1 && i< a->num_of_vectors/2 && i< b->num_of_vectors/2 ; i++)     
      if(a->vectors[i+startA]!=NULL && b->vectors[startB+i]!=NULL){
         k+= basic_kernel(kernel_parm, a, b, startA+i, startB+i)/
         sqrt(a->vectors[startA+i]->twonorm_STD * b->vectors[startB+i]->twonorm_STD);
      }
   return k;
}

 
/***************************************************************************************/
/*                                  KERNELS COMBINATIONS                               */
/***************************************************************************************/
 
// select the method to combine a forest of trees
// when will be available more kernel types, remeber to define a first_kernel option (e.g. -F)

double choose_tree_kernel(KERNEL_PARM *kernel_parm, DOC *a, DOC *b){
     /* calculate the kernel function */

  switch(kernel_parm->vectorial_approach_tree_kernel) {

    case 'S': /* TREE KERNEL Sequence k11+k22+k33+..+knn*/
            return sequence_tree_kernel(kernel_parm,a,b);; 
    case 'A': /* TREE KERNEL ALL-vs-ALL k11+k12+k13+..+k23+k33+..knn*/
            return(AVA_tree_kernel(kernel_parm,a,b));         	     
    case 'R': /* re-ranking kernel classic SST*/
            return((CFLOAT)tree_kernel_ranking(kernel_parm,a,b));
//    case 7: /* TREE KERNEL MAX of ALL-vs-ALL */
//            return(AVA_MAX_tree_kernel(kernel_parm,a,b));         	     
//    case 8: /* TREE KERNEL MAX of sequence of pairs Zanzotto et all */
//            return(AVA_MAX_tree_kernel_over_pairs(kernel_parm,a,b));         	     
    default: printf("Error: Unknown tree kernel function\n"); 
		     exit(1);
   }
}


// select the method to combine the set of vectors

double choose_second_kernel(KERNEL_PARM *kernel_parm, DOC *a, DOC *b) 
     /* calculate the kernel function */
{
  switch(kernel_parm->vectorial_approach_standard_kernel) {

   case 'S':/* non structured KERNEL Sequence k11+k22+k33+..+knn*/
            return(sequence(kernel_parm,a,b)); 
   case 'A': /* Linear KERNEL ALL-vs-ALL k11+k12+k13+..+k23+k33+..knn*/
            return(AVA(kernel_parm,a,b));
   case 'R': /* re-ranking kernel*/
            return((CFLOAT)vector_ranking(kernel_parm,a,b));
         	     
//    case 13: /* Linear KERNEL MAX of ALL-vs-ALL */
//            return((CFLOAT)AVA_MAX(kernel_parm,a,b));         	     
//    case 14: /* TREE KERNEL MAX of sequence of pairs Zanzotto et all */
//            return((CFLOAT)AVA_MAX_over_pairs(kernel_parm,a,b));         	     
    default: printf("Error: Unknown kernel combination function\n"); 
		     fflush(stdout);
		     exit(1);
   }
}


// select the data to be used in kenrels:
//            vectors, trees, their sum or their product

double advanced_kernels(KERNEL_PARM * kernel_parm, DOC * a, DOC * b){

  double k1,
         k2;
/* TEST
        tmp = (k1*k2);
     	printf("K1 %f and K2= %f NORMA= %f norma.a= %f  norma.b= %f\n",k1,k2,norma,a->twonorm_sq,b->twonorm_sq);
	printf("\nKernel Evaluation: %1.20f\n", tmp);
*/

 switch(kernel_parm->combination_type) {

    case '+': /* sum first and second kernels*/
              k1 = choose_tree_kernel(kernel_parm, a, b);
              k2 = choose_second_kernel(kernel_parm, a, b);
    	      return k2 + kernel_parm->tree_constant*k1;
    case '*': k1 = choose_tree_kernel(kernel_parm, a, b);
              k2 = choose_second_kernel(kernel_parm, a, b);
              return k1*k2;
    case 'T': /* only trees */
              return choose_tree_kernel(kernel_parm, a, b);
    case 'V': /* only vectors*/
              return choose_second_kernel(kernel_parm, a, b); 
              // otherwise evaluate the vectorial kernel on the basic kernels
    default: printf("Error: Unknown kernel combination\n"); 
		     fflush(stdout);
		     exit(1);
   }
}
