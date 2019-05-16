#ifndef _enumerate_significant_itemsets_c_
#define _enumerate_significant_itemsets_c_

/* LIBRARY INCLUDES */
#include <math.h>


//#define min(X,Y) (((X) < (Y)) ? (X) : (Y))
//#define max(X,Y) (((X) > (Y)) ? (X) : (Y))

inline double max(double X , double Y){
	return (X>Y) ? X : Y;
}

inline double min(double X , double Y){
	return (X>Y) ? Y : X;
}

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"var_declare.h"
#include"transaction_keeping.c"
#include"lcm_var.c"

/* CONSTANT DEFINES */
#define READ_BUF_SIZ 524288 //Size of the buffer to read chars from file

/* GLOBAL VARIABLES */
FILE* results_file;
// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
int N, N_over_2;
// Number of observations in positive class
int n;
int N_minus_n1;
// testable region
int xub;
int xlb;
int testable_low;
int testable_high;
// Number of non-empty transaction
int Neff;
// Original vector of labels (dimension # non-empty transactions)
char *labels;

// Corrected significance threshold
double delta;
// frequency accuracy
double epsilon;
// Array with all values of log(n!) in [0,N] pre-computed
double *loggamma;
// Logarithm of 1/binom(N,n). This terms appears for every evaluation of the hypergeometric
// PDF, so it makes sense to precompute
double log_inv_binom_N_n;
double *psi; // montone version
double *psi_exact; // exact non-monotone version
// log of psi
double *psi_log; // montone version
double *psi_log_exact; // exact non-monotone version
double *Tvalues;
double *probabilities_cache;
char *Tflag;
double *pvalues_cache;
char *pvalues_cache_flag;
int use_pvalues_cache;
double *probabilities_thr;

long count_number_of_ci_rejects;
long count_number_of_simpleub_rejects;
long count_number_of_enumerate_rejects;
long count_number_of_enumerate_nonrejects;
long explored_patterns;
long explored_contingency_tables;


double log10_2;
double log_10;

// Output files
FILE *significant_itemsets_output_file, *pvalues_output_file;

/* FUNCTION DECLARATIONS */
void loggamma_init();
void psi_init();
void get_N_n(char *);
void read_labels_file(char *, char*);
extern void LCMFREQ_output_itemset(int *);
// Profiling variables
long long n_significant_patterns;



#include "unconditional.c"


/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the code
 * Input arguments are self-explanatory
 * */
void enum_sig_itemsets_init(double sig_th, char *labels_file , int min_supp , double epsilon_){
	int j; //Loop variable
	char *labels_buffer;

	get_N_n(labels_file);

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	delta = sig_th;

	// Allocate memory for the buffer containing the class labels, giving an error if it fails
	labels_buffer = (char *)malloc(N*sizeof(char));
	if(!labels_buffer){
		fprintf(stderr,"Error in function enum_sig_itemsets_init: couldn't allocate memory for array labels_buffer\n");
		exit(1);
	}

	// Allocate memory for significance threshold looking only at observed contingency tables
	probabilities_thr = (double *)malloc(N*sizeof(double));
	if(!probabilities_thr){
		fprintf(stderr,"Error in function enum_sig_itemsets_init: couldn't allocate memory for array probabilities_thr\n");
		exit(1);
	}
	for(j=0; j<N; j++){
		probabilities_thr[j] = 0.0;
	}

	use_pvalues_cache = 1;
	double order_of_magnitude = log(n+1) + log(N-n+1);
	order_of_magnitude = order_of_magnitude / log(10.0);
	unsigned long long size_of_T_tables = 0;
	if(order_of_magnitude > 9.0){
		use_pvalues_cache = 0;
		fprintf(stderr,"disabled p-values cache \n");
	}
	else{
		size_of_T_tables = (n+1) * (N-n+1);
		fprintf(stderr,"p-values cache on , size_of_T_tables %ld\n",size_of_T_tables);
	}

	if(use_pvalues_cache == 1){
		pvalues_cache = (double *)malloc(size_of_T_tables*sizeof(double));
		if(!pvalues_cache){
			fprintf(stderr,"Error: couldn't allocate %llu memory for array pvalues_cache\n",size_of_T_tables);
			use_pvalues_cache = 0;
			//exit(1);
		}
	}
	if(use_pvalues_cache == 1){
		pvalues_cache_flag = (char *)malloc(size_of_T_tables*sizeof(char));
		if(!pvalues_cache_flag){
			fprintf(stderr,"Error: couldn't allocate  %llu memory for array pvalues_cache_flag\n",size_of_T_tables);
			use_pvalues_cache = 0;
			//exit(1);
		}
	}

	if(use_pvalues_cache == 1){
		unsigned long cache_index = 0;
		for(cache_index=0; cache_index<size_of_T_tables; cache_index++){
			pvalues_cache_flag[cache_index] = 0;
		}
	}

	/* Allocate memory for the vector of class labels, with labels of empty transactions removed */
	Neff = root_trans_list.siz1;
	labels = (char *)malloc(Neff*sizeof(char));
	if(!labels){
		fprintf(stderr,"Error in function enum_sig_itemsets_init: couldn't allocate memory for array labels\n");
		exit(1);
	}

	// Read file containing class labels and store them in array labels, taking care of removing labels
	// associated with empty transactions
	read_labels_file(labels_file,labels_buffer);
	// Ensure class 1 is the minority class
	if(n > (N/2)){
		for(j=0; j<N; j++) labels_buffer[j] = !labels_buffer[j];
		n = N-n;
	}
	N_minus_n1 = N - n;
	fprintf(stderr,"N = %d\n",N);
	fprintf(stderr,"n1 = %d\n",n);
	//fprintf(stderr,"N_minus_n1 = %d\n",N_minus_n1);

	epsilon = epsilon_;

	for(j=0;j<Neff;j++) labels[j] = labels_buffer[non_empty_trans_idx[j]];
	free(labels_buffer);
	// The array containing the indices of all non-empty transactions is no longer needed
	free(non_empty_trans_idx);

	// Initialise cache for log(x!) and psi(x)
	loggamma_init();
	psi_init();

	testable_low = min_supp;
	testable_high = N-1;

	// Initialise profiling variables
	n_significant_patterns = 0;
	count_number_of_ci_rejects = 0;
	count_number_of_simpleub_rejects = 0;
	count_number_of_enumerate_rejects = 0;
	count_number_of_enumerate_nonrejects = 0;
	explored_patterns = 0;
	explored_contingency_tables = 0;

	#ifdef DVERSION4
	// debug binomial functions
	{
		double pi = 0.3;
		int a0 = 10;
		int n0_ = 10;
		double tail_test = -10000;
		double tail_test2 = 0.0;
		double p_table = 0.0;
		int i = 0;
		for(i=a0 ; i<=n0_ ; i++){
			p_table = (i*log(pi)) + ((n0_-i)*log(1-pi)) + logbincoeff(n0_,i);
			tail_test = sumlogs(tail_test , p_table);
		}
		tail_test = exp(tail_test);
		tail_test2 = incbeta((double)a0, (double)(n0_-a0+1), pi);
		fprintf(stderr, "P(X >= %d) = %e", a0 ,tail_test2 );
		fprintf(stderr, "  , tail_test = %e\n", tail_test );
	}
	{
		double pi = 0.5;
		int a0 = 9000;
		int n0_ = 10000;
		double tail_test = -10000;
		double tail_test2 = 0.0;
		double p_table = 0.0;
		int i = 0;
		for(i=a0 ; i<=n0_ ; i++){
			p_table = (i*log(pi)) + ((n0_-i)*log(1-pi)) + logbincoeff(n0_,i);
			tail_test = sumlogs(tail_test , p_table);
		}
		tail_test = exp(tail_test);
		tail_test2 = incbeta((double)a0, (double)(n0_-a0+1), pi);
		fprintf(stderr, "P(X >= %d) = %e", a0 ,tail_test2 );
		fprintf(stderr, "  , tail_test = %e\n", tail_test );
	}
	{
		double pi = 0.5;
		int a0 = 1;
		int n0_ = 100;
		for(a0=0 ; a0<=n0_ ; a0++){
			double tail_test = -10000;
			double tail_test2 = 0.0;
			double p_table = 0.0;
			int i = 0;
			for(i=a0 ; i<=n0_ ; i++){
				p_table = (i*log(pi)) + ((n0_-i)*log(1-pi)) + logbincoeff(n0_,i);
				tail_test = sumlogs(tail_test , p_table);
			}
			tail_test = exp(tail_test);
			tail_test2 = incbeta((double)a0, (double)(n0_-a0+1), pi);
			fprintf(stderr, "P(X >= %d) = %e", a0 ,tail_test2 );
			fprintf(stderr, "  , tail_test = %e\n", tail_test );
		}
	}
	#endif
}

/* Precompute values of log(x!) storing them in the array loggamma */
void loggamma_init(){
	int x;
	// Allocate memory for log-gamma cache, raising error if it fails
	loggamma = (double *)malloc((N+1)*sizeof(double));
	if(!loggamma){
		fprintf(stderr,"Error in function loggamma_init: couldn't allocate memory for array loggamma\n");
		exit(1);
	}
	// Initialise cache with appropriate values
	for(x=0;x<=N;x++) loggamma[x] = lgamma(x+1);//Gamma(x) = (x-1)!
	// Initialise log_inv_binom_N_n
	log_inv_binom_N_n = loggamma[n] + loggamma[N-n] - loggamma[N];
}


/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
void psi_init(){
	double xi1;
	int x, x_init;
	// Allocate memory for psi, raising error if it fails
	psi = (double *)malloc((N+1)*sizeof(double));
  //psi_log = (double *)malloc((N+1)*sizeof(double));
  //psi_exact = (double *)malloc((N+1)*sizeof(double));
  //psi_log_exact = (double *)malloc((N+1)*sizeof(double));
	if(!psi){// || !psi_log){
		fprintf(stderr,"Error in function psi_and_xi1_init: couldn't allocate memory for array psi\n");
		exit(1);
	}

	int aj = 1;
	xub = 0;
	xlb = N;
	psi[0] = 1.0;
	// don't used when enumerating significant patterns since LCM knows the minimum support of candidates
	for(aj = 1; aj < (N+1) ; aj++){
		psi[aj] = 0.0;
	}

}

/* Free all allocated memory and give some output for debugging purposes */
void enum_sig_itemsets_end(){
	// Print results
	fprintf(results_file,"RESULTS\n");
	fprintf(results_file,"\t Corrected significance threshold: %e\n",delta);
	fprintf(results_file,"\t LCM support: %d\n",LCM_th);
	fprintf(results_file,"\t Number of significant patterns found: %lld\n",n_significant_patterns);

	fprintf(stderr, "Tested patterns %ld\n", explored_patterns);
	fprintf(stderr, "explored_contingency_tables %ld\n", explored_contingency_tables);
	fprintf(stderr, "Number of significant patterns found: %lld\n",n_significant_patterns);
	fprintf(stderr, "count_number_of_ci_rejects %ld\n", count_number_of_ci_rejects);
	fprintf(stderr, "count_number_of_simpleub_rejects %ld\n", count_number_of_simpleub_rejects);
	fprintf(stderr, "count_number_of_enumerate_rejects %ld\n", count_number_of_enumerate_rejects);
	fprintf(stderr, "count_number_of_enumerate_nonrejects %ld\n", count_number_of_enumerate_nonrejects);

	// Free allocated memory
	free(loggamma); free(psi);
	free(labels);
	free(Tvalues);
	free(Tflag);
	free(probabilities_cache);

	// Close output files
	fclose(results_file);
	fclose(significant_itemsets_output_file);
	fclose(pvalues_output_file);

}


double fisher_pval(int a, int x){
	int a_min, a_max, k;
	double p_left, p_right, pval;
	double pre_comp_xterms;

	// Compute the contribution of all terms depending on x but not on a
	pre_comp_xterms = loggamma[x] + loggamma[N-x];
	a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n) ? n : x;//min(x,n)

	// The algorithm to compute the p-value proceeds as follows. We inspect simultaneously probability values on the left and right tails of the
	// hypergeometric distribution, "accepting" each time the one that is smaller. When that happens, we move the index in the appropriate direction,
	// that is, a_min++ if we "accept" a value on the left and a_max-- if we "accept" a value on the right. As soon as the "accepted" value is located
	// in index a, we know that we have already explored all values of the hypergeometric probability mass whose probabilities are smaller or equal
	// than the probability of a. The only tricky case occurs when the probability values on both sides happen to be identical. The way to handle
	// that case is by "accepting" both values simultaneously.
	pval = 0; //Accumulate probabilities in this variable
	while(a_min<a_max){
		p_left = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_min] + loggamma[n-a_min] + loggamma[x-a_min] + loggamma[(N-n)-(x-a_min)]));
		p_right = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_max] + loggamma[n-a_max] + loggamma[x-a_max] + loggamma[(N-n)-(x-a_max)]));
		if(p_left == p_right) {
			pval += (p_left+p_right);
			if((a==a_min) || (a==a_max)) return pval;
			a_min++; a_max--;
		}
		else if(p_left < p_right){
			pval += p_left;
			if(a==a_min) return pval;
			a_min++;
		}
		else{
			pval += p_right;
			if(a==a_max) return pval;
			a_max--;
		}
	}
	// If we get to this part of the code, it means is the mode of the distribution and, therefore, its associated p-value is 1
	return 1;
}

/* -------------------FUNCTIONS TO PROCESS A NEWLY FOUND TESTABLE HYPOTHESIS-------------------------------------- */

/* This code contains 3 difference functions to process newly found hypotheses. All of them are virtually identical
 * and the only thing which differs is the way the function receives the list of observations (transactions) for
 * which the hypothesis has X=1.
 * LCM has a complex structure, with frequent itemsets being found at 4 different locations in the source code
 * and under 3 different circumstances. Therefore it was needed to introduce differentiation in the way the transaction
 * list is fed to the "solution processing functions" in order to keep the "transaction keeping" overhead minimal.
 *
 * To reuse this code for problems other than frequent itemset mining, the only thing that needs to be modified
 * is the line which computes the cell counts, for example, the following line in bm_process_solution:
 * 		for(i=0; i<current_trans.siz; i++) a += labels_perm[j][current_trans.list[i]];
 * 	There, current_trans.siz is the number of transactions for which the hypothesis has X=1, i.e. the margin x
 * 	of the 2x2 contingency table (note in this case it is redundant with the input argument x of the function)
 * 	Similarly, current_trans.list[i] with i ranging from 0 to (x-1) has the list of indices of the observations
 * 	for which X=1.
 * 	Simply changing those two parameters accordingly allows code reuse.
 * */

/* Process a solution involving the bitmap represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
void bm_process_solution(int x, int item, int *mask){
	int i,j;//Loop iterators
	int a; //Cell count of current itemset
	double pvalue_lowerbound , pvalue_upperbound;//P-value of current itemset

	// Sanity-check
	if (x != current_trans.siz) printf("Error: x = %d, current_trans.siz=%d\n",x,current_trans.siz);

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	// Compute the cell-count corresponding to the current itemset
	a = 0;
	for(i=0; i<current_trans.siz; i++) a += labels[current_trans.list[i]];
	// Compute the corresponding p-value
	//pval = computeUncondUBPValue(x , a);//computePvalueBarnard(x , a , 1.0); //fisher_pval(a,x);
	#ifdef VERSION2
	pvalue_lowerbound = computeUncondUBPValue_version2(x , a , epsilon);
	pvalue_upperbound = 1.0;
	#endif
	#ifdef VERSION3
	computeUncondUBPValue_version3(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif
	#ifdef VERSIONNAIVE
	compute_naive_ut(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif
	#ifdef VERSION4
	unconditional_fast_test(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif

	#ifdef DVERSION4
	double pvalue_lowerbound3, pvalue_upperbound3 , pvalue_lowerboundn, pvalue_upperboundn;
	compute_naive_ut(x , a , epsilon , &pvalue_lowerboundn , &pvalue_upperboundn);
	computeUncondUBPValue_version3(x , a , epsilon , &pvalue_lowerbound3 , &pvalue_upperbound3);
	if(pvalue_lowerboundn/pvalue_lowerbound > 1.001 || pvalue_lowerboundn/pvalue_lowerbound < 0.999){
		fprintf(stderr,"x = %d, a = %d \n",x,a);
		fprintf(stderr,"  naive %e,%e\n",pvalue_lowerboundn,pvalue_upperboundn);
		fprintf(stderr,"    V3 %e,%e\n",pvalue_lowerbound3,pvalue_upperbound3);
		fprintf(stderr,"    V4 %e,%e\n",pvalue_lowerbound,pvalue_upperbound);
	}
	#endif

	// If p-value is significant, write current itemset and the corresponding p-value to the output files
	if(pvalue_lowerbound <= delta){
		//printf("pval = %f\n",pval);
		n_significant_patterns++;
		fprintf(pvalues_output_file,"%d,%d,%d,%e,%e,%e\n",a,x-a,x,computeP0(x , a),pvalue_lowerbound,pvalue_upperbound);
		fprintf_current_itemset();
	}


}

/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0(int x){
	int i,j;//Loop iterators
	int a; //Cell count of current itemset
	double pvalue_lowerbound , pvalue_upperbound;//P-value of current itemset

	// Sanity-check
	if (x != bm_trans_list[1].siz) printf("Error: x = %d, bm_trans_list[1].siz=%d\n",x,bm_trans_list[1].siz);

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	// Compute the cell-count corresponding to the current itemset
	a = 0;
	for(i=0; i<bm_trans_list[1].siz; i++) a += labels[bm_trans_list[1].list[i]];
	// Compute the corresponding p-value
	//pval = computeUncondUBPValue(x , a);//computePvalueBarnard(x , a , 1.0); //fisher_pval(a,x);
	#ifdef VERSION2
	pvalue_lowerbound = computeUncondUBPValue_version2(x , a , epsilon);
	pvalue_upperbound = 1.0;
	#endif
	#ifdef VERSION3
	computeUncondUBPValue_version3(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif
	#ifdef VERSIONNAIVE
	compute_naive_ut(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif

	#ifdef VERSION4
	unconditional_fast_test(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif


	#ifdef DVERSION4
	double pvalue_lowerbound3, pvalue_upperbound3 , pvalue_lowerboundn, pvalue_upperboundn;
	compute_naive_ut(x , a , epsilon , &pvalue_lowerboundn , &pvalue_upperboundn);
	computeUncondUBPValue_version3(x , a , epsilon , &pvalue_lowerbound3 , &pvalue_upperbound3);
	if(pvalue_lowerboundn/pvalue_lowerbound > 1.001 || pvalue_lowerboundn/pvalue_lowerbound < 0.999){
		fprintf(stderr,"x = %d, a = %d \n",x,a);
		fprintf(stderr,"  naive %e,%e\n",pvalue_lowerboundn,pvalue_upperboundn);
		fprintf(stderr,"    V3 %e,%e\n",pvalue_lowerbound3,pvalue_upperbound3);
		fprintf(stderr,"    V4 %e,%e\n",pvalue_lowerbound,pvalue_upperbound);
	}
	#endif

	// If p-value is significant, write current itemset and the corresponding p-value to the output files
	if(pvalue_lowerbound <= delta){
		//printf("pval = %f\n",pval);
		n_significant_patterns++;
		fprintf(pvalues_output_file,"%d,%d,%d,%e,%e,%e\n",a,x-a,x,computeP0(x , a),pvalue_lowerbound,pvalue_upperbound);
		LCMFREQ_output_itemset(LCM_add.q+LCM_add.t);
	}
}

/* Process a solution involving the array-list represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// L = pointer to TRANS_LIST struct keeping track of merged transactions
// item = current node of the tree
void ary_process_solution(int x, TRANS_LIST *L, int item, int *mask){
	int j;//Loop iterator
	int aux; //Auxiliary counter
	int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list
	int a; //Cell count of current itemset
	double pvalue_lowerbound , pvalue_upperbound;//P-value of current itemset

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	// Compute the cell-count corresponding to the current itemset, plus sanity-check
	a = 0; aux = 0;
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
			a += labels[*ptr];
			aux++;
		}
	}
	if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);
	// Compute the corresponding p-value
	//pval = computeUncondUBPValue(x , a);//computePvalueBarnard(x , a , 1.0); //fisher_pval(a,x);
	//pval = computeUncondUBPValue_version2(x , a , epsilon);
	#ifdef VERSION2
	pvalue_lowerbound = computeUncondUBPValue_version2(x , a , epsilon);
	pvalue_upperbound = 1.0;
	#endif
	#ifdef VERSION3
	computeUncondUBPValue_version3(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif
	#ifdef VERSIONNAIVE
	compute_naive_ut(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif

	#ifdef VERSION4
	unconditional_fast_test(x , a , epsilon , &pvalue_lowerbound , &pvalue_upperbound);
	#endif


	#ifdef DVERSION4
	double pvalue_lowerbound3, pvalue_upperbound3 , pvalue_lowerboundn, pvalue_upperboundn;
	compute_naive_ut(x , a , epsilon , &pvalue_lowerboundn , &pvalue_upperboundn);
	computeUncondUBPValue_version3(x , a , epsilon , &pvalue_lowerbound3 , &pvalue_upperbound3);
	if(pvalue_lowerboundn/pvalue_lowerbound > 1.001 || pvalue_lowerboundn/pvalue_lowerbound < 0.999){
		fprintf(stderr,"x = %d, a = %d \n",x,a);
		fprintf(stderr,"  naive %e,%e\n",pvalue_lowerboundn,pvalue_upperboundn);
		fprintf(stderr,"    V3 %e,%e\n",pvalue_lowerbound3,pvalue_upperbound3);
		fprintf(stderr,"    V4 %e,%e\n",pvalue_lowerbound,pvalue_upperbound);
	}
	#endif

	// If p-value is significant, write current itemset and the corresponding p-value to the output files
	if(pvalue_lowerbound <= delta){
		//printf("pval = %f\n",pval);
		n_significant_patterns++;
		fprintf(pvalues_output_file,"%d,%d,%d,%e,%e,%e\n",a,x-a,x,computeP0(x , a),pvalue_lowerbound,pvalue_upperbound);
		fprintf_current_itemset();
	}

}

/* Do a first scan of the file containing the class labels to compute the total number of observations, N,
 * and the total number of observations in the positive class, n
 * */
void get_N_n(char *labels_file){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops

	// Initialise both counters to 0 (the variables are defined as global variables in wy.c)
	N = 0; n = 0;

	//Try to open file, giving an error message if it fails
	if(!(f_labels = fopen(labels_file,"r"))){
		fprintf(stderr, "Error in function get_N_n when opening file %s\n",labels_file);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function read_labels_file: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
			fprintf(stderr,"Error in function read_labels_file while reading the file %s\n",labels_file);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			N++;
			if(char_to_int[*read_buf_aux]) n++;
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
}

void read_labels_file(char *labels_file, char *labels_buffer){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
	char *labels_aux = labels_buffer;//Auxiliary pointer to array labels for increments

	//Try to open file, giving an error message if it fails
	if(!(f_labels = fopen(labels_file,"r"))){
		fprintf(stderr, "Error in function read_labels_file when opening file %s\n",labels_file);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function read_labels_file: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
			fprintf(stderr,"Error in function read_labels_file while reading the file %s\n",labels_file);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			*labels_aux++ = char_to_int[*read_buf_aux];
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
}

#endif
