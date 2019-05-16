// source code for unconditional test

// this enables version 3
//#define VERSION3 1
//#define USECONFINTERVALS 1
// this enables version 4
#define VERSION4 1
#define VERSION4UBONLY 1
//#define DVERSION4 1
//#define VERSIONNAIVE 1

#include "incbeta.c"

double fast_lower_bound_box(int x_s , int a_s);


inline double sumlogs(double first_log , double second_log){
	return max(first_log , second_log) + (log1p(exp((min(first_log , second_log) - max(first_log , second_log)))));
}

inline double sublogs(double first_log , double second_log){
	return max(first_log , second_log) + (log(1.0 - exp((min(first_log , second_log) - max(first_log , second_log)))));
}


inline double logbincoeff(int n_, int k_){
	return loggamma[n_] - loggamma[n_-k_] - loggamma[k_];
}


int check_bounds(int x_s , int a_s , double epsilon , double *pvalue_lowerbound , double *pvalue_upperbound){

		  if (x_s == N){
		 		*pvalue_lowerbound = 1.0;
		 		*pvalue_upperbound = 1.0;
				return 1;
			}

			int DEBUG_L = 0;
			if(DEBUG_L){fprintf(stderr,"check_bounds: (x_s = %d , a_1 = %d , a_0 = %d)\n",x_s,a_s,x_s-a_s);}

			// check simple upper bound
		  double pi = (double)x_s / (double)N;
			double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi));
		  int n1 = n;
			p0 = p0 + logbincoeff(N_minus_n1,x_s-a_s) + logbincoeff(n1,a_s);

			// check lower bound from p0
			if(exp(p0) > delta){
				if(DEBUG_L){fprintf(stderr,"  not reject due to p0 (%e)\n",exp(p0));}
				*pvalue_lowerbound = exp(p0);
			 	*pvalue_upperbound = 1.0;
				return 1;
			}

			// check confidence intervals
			double lb_ci_pi; double ub_ci_pi;
			if(confidenceIntervalsNHReject(x_s , a_s , epsilon , &lb_ci_pi , &ub_ci_pi)){
				if(DEBUG_L){fprintf(stderr,"  reject using conf intervals\n");}
				count_number_of_ci_rejects++;
			 	*pvalue_lowerbound = 0.0;
			 	*pvalue_upperbound = 0.0;
			 	return 1;
			}

			double test_upper_bound = p0 + log((double)(n1 + 1)) + log((double)(N-n1 + 1));
			test_upper_bound = exp(test_upper_bound);
			 if(test_upper_bound <= delta){
				if(DEBUG_L){fprintf(stderr,"Rejected with UB2 (%e)!\n",test_upper_bound);}
				 count_number_of_simpleub_rejects++;
		 		*pvalue_lowerbound = exp(p0);
		 	 	*pvalue_upperbound = test_upper_bound;
		 		return 1;
			}

			// check probability thresholds
			if(probabilities_thr[x_s] <= p0){
				if(DEBUG_L){fprintf(stderr,"  not reject due to p0 (%e)\n",exp(p0));}
			 	*pvalue_lowerbound = 1.0;
			 	*pvalue_upperbound = 1.0;
				return 1;
			}

			// check if this p-value was aready computed
			if(use_pvalues_cache == 1){
				unsigned long pcache_index = ((unsigned long)a_s * (unsigned long)(N-n)) + (x_s-a_s);
				if(pvalues_cache_flag[pcache_index]){
					double p_value_toreturn = pvalues_cache[pcache_index];
						if(p_value_toreturn <= delta){
							count_number_of_enumerate_rejects++;
						 if(DEBUG_L){fprintf(stderr,"Rejected!\n");}
					 }
					*pvalue_lowerbound = p_value_toreturn;
					*pvalue_upperbound = test_upper_bound;
					return 1;
				}
			}


			// check upper bound on approximation
			/*double dn0 = (double)n0;
			double dn1 = (double)n1;
			double variance1 = dn1 * pi * (1.0 - pi);
			double variance0 = dn0 * pi * (1.0 - pi);
			double test_ub_approx = 4.0 * variance0 + (2.0 * dn0 * pi - (double)());
			if(test_ub_approx <= delta{
				*pvalue_lowerbound = test_ub_approx;
				*pvalue_upperbound = test_upper_bound;
				return 1;
			}*/

			/*double lower_bound_using_squared_box = fast_lower_bound_box(x_s , a_s);
			fprintf(stderr,"lower_bound_using_squared_box = %e , exp(p0) = %e , p0 = %f\n",lower_bound_using_squared_box,exp(p0),p0);
			if(lower_bound_using_squared_box > delta){
			 fprintf(stderr,"Not Rejected with FastUB (%e)!\n",lower_bound_using_squared_box);
			 *pvalue_lowerbound = lower_bound_using_squared_box;
			 *pvalue_upperbound = test_upper_bound;
			 return 0;
		 }*/

			return 0;
}



double computeMinPvalueLBBarnard_version2(int x_s , double epsilon){

	double pi = (double)x_s / (double)N;

	 int n1 = n;
	 double p_value = -10000; // 0 in log scale
	 double p_table = 0.0;
	 double p_table_precomp = 0.0;
	 int a = 1;
	 int maxAS_ = 1;
	 int minAS_ = 1;
	 int epsilon_transactions = (int)(epsilon * N);
	 int a_s = 0;

	 int DEBUG_B = 0;

	  if(DEBUG_B){
	 	 fprintf(stdout,"Starting Improved Barnard Psi: (x_s = %d , N = %d, n1 = %d)\n",x_s,N,n);
	 	 fprintf(stdout,"   pi = %f \n",pi);
	  }
		double p0 = 1.0;

		// find the minimum p0
		maxAS_ = min(x_s,n1);
		minAS_ = max(0,x_s-(N-n1));
		p_table_precomp = (x_s*log(pi)) + ((N-x_s)*log(1-pi));
		for (a = minAS_; a <= maxAS_; a++){
			p_table = p_table_precomp + logbincoeff(N-n1,x_s-a) + logbincoeff(n1,a);
			if(p_table <= p0){
				p0 = p_table;
				a_s = a;
			}
		}

	  if(DEBUG_B){
	 	 fprintf(stdout,"   p0 = %f \n",p0 / log_10);
	  }

	 	 //cout << "p0 = " << p0 << " with pi = " << pi << endl;

	 	 // generate all possible contingency tables with n1 fixed
	 	 //fprintf(stdout,"Starting generating cont tables\n");
		 int xs = 1;
	 	 // confidence intervals on xs using epsilon
	 	 int xlb_ = max(1 , x_s - epsilon_transactions);
	 	 int xub_ = min(N , x_s + epsilon_transactions);
	 	 for (xs = xlb_; xs <= xub_; xs++){
	 		 maxAS_ = min(xs,n1);
	 		 minAS_ = max(0,xs-(N-n1));
	 		 // confidence intervals on as using epsilon
	 		 minAS_ = min(minAS_ , a_s - epsilon_transactions);
	 		 maxAS_ = min(maxAS_ , a_s + epsilon_transactions);
	 		 p_table_precomp = (xs*log(pi)) + ((N-xs)*log(1-pi));
	 		 for (a = minAS_; a <= maxAS_; a++){
	 			 p_table = p_table_precomp + logbincoeff(N-n1,xs-a) + logbincoeff(n1,a);
	 			 if(p_table <= p0){
	 				 p_value = sumlogs(p_value , p_table);
	 				 if(DEBUG_B){
	 					 /*fprintf(stdout,"   x = %d , a = %d \n",xs,a);
	 					 fprintf(stdout,"   ptable = %f \n",p_table);
	 					 fprintf(stdout,"   pvalue increased = %f \n",p_value);*/
	 				 }
	 			 }
	 		 }
	 	 }

	 	 if(DEBUG_B){
	 		 fprintf(stdout,"   final minimum pvalue = %f \n",p_value / log_10);
	 	 }

	  return p_value;

}

double computeMinPvalueBarnard(int x_s , double pi_precision){

	int steps = max((1.0/pi_precision)-1.0 , 1.0);

	double log_e_const = log(exp(1));

	double *pis = (double*)malloc(sizeof(double)*steps);
	pis[0] = (double)x_s / (double)N;
	pis[0] = min(pis[0] , 0.9999999);
	pis[0] = max(pis[0] , 0.0000001);

	double step = 0.02 * pi_precision;

	int i = 1;
	for(i = 1; i < steps; i++){
		pis[i] = pis[0] - 0.1 + (double)(i) * step;
		pis[i] = min(pis[i] , 0.9999999);
		pis[i] = max(pis[i] , 0.0000001);
	}

	// p_values for values of pi
	double *barn_p_values = (double*)malloc(sizeof(double)*steps);

	double p_value_max = -10000000;
	double pi_star = 0;


	int n1 = n;

	int maxAS = min(x_s,n1);
	int minAS = max(0,x_s-(N-n1));


	// for every value of pi
	for(i = 0; i < steps; i++){

		double pi = pis[i];

		double log_pi = log(pi);
		double log_1_pi = log(1-pi);
		double constfactor1 = log_pi/log_e_const;
		double constfactor2 = log_1_pi/log_e_const;

		double p0 = (x_s*constfactor1) + ((N-x_s)*constfactor2);
		p0 = p0 + min(logbincoeff(N-n1,x_s-maxAS) + logbincoeff(n1,maxAS) , logbincoeff(N-n1,x_s-minAS) + logbincoeff(n1,minAS));

		barn_p_values[i] = p0;

		fprintf(stdout,"(x=%d) , p0 = %f , ",x_s,p0);

		// generate all possible contingency tables with n1 fixed

		double p_table = 0.0;
		double p_table_precomp = 0.0;

		int xs = 1;

		for (xs = 1; xs < N+1; xs++){

			int maxAS_ = min(xs,n1);
			int minAS_ = max(0,xs-(N-n1));

			p_table_precomp = (xs*constfactor1) + ((N-xs)*constfactor2);

			int a = minAS_;
			p_table = p_table_precomp + logbincoeff(N-n1,xs-a) + logbincoeff(n1,a);
			while(p_table <= p0 && a <= maxAS_){
				barn_p_values[i] = sumlogs(barn_p_values[i] , p_table);
				a++;
				p_table = p_table_precomp + logbincoeff(N-n1,xs-a) + logbincoeff(n1,a);
			}
			int low_a = a;
			a = maxAS_;
			p_table = p_table_precomp + logbincoeff(N-n1,xs-a) + logbincoeff(n1,a);
			while(p_table <= p0 && a > low_a){
				barn_p_values[i] = sumlogs(barn_p_values[i] , p_table);
				a--;
				p_table = p_table_precomp + logbincoeff(N-n1,xs-a) + logbincoeff(n1,a);
			}

			/*for(int a = minAS_; a <= maxAS_; a++){

				// compute probability of this contingency table
				p_table = p_table_precomp + logbincoeff(N-n1,xs-a) + logbincoeff(n1,a);
				if(p_table <= p0){
					barn_p_values[i] = sumlogs(barn_p_values[i] , p_table);
				}

			}*/

		}

		barn_p_values[i] = sublogs(barn_p_values[i] , p0);
		p_value_max = max (p_value_max , barn_p_values[i]);
		if(p_value_max == barn_p_values[i])
			pi_star = pi;

	}
	//out_file << x_s << ";" << pi_precision << ";" << p_value_max << ";" << pi_star << endl;

	free(barn_p_values);
	free(pis);

	fprintf(stdout,"psi = %f\n",x_s,p_value_max);

	return p_value_max;


}




 double computeUncondUBPValue(int x_s , int a_s){

	//fprintf(stdout,"Starting: (x_s = %d , a_s = %d)\n",x_s,a_s);

 	int n1 = n;
 	double p_value = -10000;

	double T0 = ((double)a_s / (double)n1) - ((double)(x_s - a_s) / (double)(N - n1));
	double Txa = T0;
	//fprintf(stdout,"   T0 = %f\n",T0);

 	double pi = (double)x_s / (double)N;
	//fprintf(stdout,"   pi = %f \n",pi);
 	double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi));
	//fprintf(stdout,"   p0precomp = %f\n",p0);
 	//cout << "N-n1 = " << N-n1 <<  endl;
 	//cout << "x_s-a_s = " << x_s-a_s <<  endl;
 	//cout << "n1 = " << n1 <<  endl;
 	//cout << "a_s = " << a_s <<  endl;
 	p0 = p0 + logbincoeff(N-n1,x_s-a_s) + logbincoeff(n1,a_s);

	//fprintf(stdout,"   p0 = %f \n",p0);

 	if(exp(p0) > delta)
 		return 1.0;

 		//cout << "p0 = " << p0 << " with pi = " << pi << endl;

 		// generate all possible contingency tables with n1 fixed
		//fprintf(stdout,"Starting generating cont tables\n");
 		double p_table = 0.0;
 		double p_table_precomp = 0.0;

 		int continue_loop = 1;

 		int xs = 1;
 		for (xs = xlb; xs < xub+1; xs++){

 			int maxAS_ = min(xs,n1);
 			int minAS_ = max(0,xs-(N-n1));

			//fprintf(stdout,"   x = %d \n",xs);

			pi = (double)xs / (double)N;
			//fprintf(stdout,"   pi = %f \n",pi);

 			p_table_precomp = (xs*log(pi)) + ((N-xs)*log(1-pi));
			//fprintf(stdout,"   p_table_precomp = %f \n",p_table_precomp);

 			int a = maxAS_;
			Txa = ((double)a/(double)n1) - ((double)(xs - a)/(double)(N - n1));
			continue_loop = Txa >= T0;
			//fprintf(stdout,"   Txa = %f\n",Txa);
 			while(continue_loop){
				//fprintf(stdout,"   a = %d \n",a);
	 			p_table = p_table_precomp + logbincoeff(N-n1,xs-a) + logbincoeff(n1,a);
				//fprintf(stdout,"   ptable = %f , p_value = %f \n",p_table,p_value);
 				p_value = sumlogs(p_value , p_table);
 				a--;
				Txa = ((double)a/(double)n1) - ((double)(xs - a)/(double)(N - n1));
				//fprintf(stdout,"   Txa = %f\n",Txa);
				continue_loop = Txa >= T0;
 				continue_loop = continue_loop && a >= minAS_ && exp(p_value) <= delta;
 			}

			if(exp(p_value) > delta)
				return 1.0;

 		}

	//fprintf(stdout,"Finished: (xs = %d , as = %d) , p0 = %f , p = %f , cl = %d \n",x_s,a_s,p0,p_value,continue_loop);

 	return p_value / log_10;


 }




  int confidenceIntervalsNHReject(int x_s , int a_s , double epsilon , double *lb_ci_pi , double *ub_ci_pi){

 	 int DEBUG_B = 0;
 	 	// check that pi is included in intersection of Ha0 and Ha1
 	 		double a0 = x_s - a_s;
 	 		double a1 = a_s;
 	 		double n0 = N-n;
 	 		double n1 = n;
 	 		double f1 = a1 / n1;
 	 		double f0 = a0 / n0;
 			double epsilon1 = epsilon * (n0+n1) / n1;
 			double epsilon0 = epsilon * (n0+n1) / n0;
 	 		double ilb = (f0 < f1) ? f0 + epsilon0 : f1 + epsilon1;
 			double iub = (f0 > f1) ? f0 - epsilon0 : f1 - epsilon1;
 			double pi = (double)x_s / (double)N;
  			if(DEBUG_B){fprintf(stderr,"Improved Barnard 3: Conf int (x_s = %d , a_1 = %d , a_0 = %d)\n",x_s,a_s,(x_s-a_s));
  			fprintf(stderr,"Ha0 = [%f , %f] , Ha1 = [%f , %f] \n",f0-epsilon0,f0+epsilon0,f1-epsilon1,f1+epsilon1);}
 	 		if( ilb < iub ){
 	 			if(DEBUG_B){fprintf(stderr,"REJECTED! it holds %f < %f \n",ilb,iub);}
 	 			return 1; // reject the null hypothesis with p-value = 0
 	 		}
 	 		// now check that pi is in the intersection
 	 		if( !(iub <= pi && pi <= ilb) ){
 	    	 if(DEBUG_B)fprintf(stderr,"pi is not in the intersection! it should hold %f <= %f <= %f \n",iub,pi,ilb);
 	 		 pi = (pi < iub) ? iub : ilb;
 	 		 if(DEBUG_B)fprintf(stderr,"pi set to %f\n",pi);
 	 		}
 			// we have iub <= pi <= ilb
 			*lb_ci_pi = iub;
 			*ub_ci_pi = ilb;
 		return 0;

  }



	   // unconditional test version 2: test statistic is given by p-value of observed contingency table,
	 	// we find contingency tables such that their probability is =< p0

	 double computeUncondUBPValue_version2(int x_s , int a_s , double epsilon){

	 int DEBUG_B = 0;

	 if (x_s == N) return 1.0;

	 int epsilon_transactions = (int)(epsilon * (double)N);
	 #ifndef USECONFINTERVALS
	 epsilon_transactions = N;
	 #endif

	double lb_ci_pi; double ub_ci_pi;
	 if(confidenceIntervalsNHReject(x_s , a_s , epsilon , &lb_ci_pi , &ub_ci_pi)) return 0.0;

	 double pi = (double)x_s / (double)N;

	 int n1 = n;
	 double p_value = -10000; // 0 in log scale

	  if(DEBUG_B){
	 	 fprintf(stdout,"Starting Improved Barnard: (x_s = %d , a_s = %d , N = %d, n1 = %d)\n",x_s,a_s,N,n);
	 	 fprintf(stdout,"   pi = %f \n",pi);
	  }

	  double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi));
	  p0 = p0 + logbincoeff(N_minus_n1,x_s-a_s) + logbincoeff(n1,a_s);

	  if(DEBUG_B){
	 	 fprintf(stdout,"   delta = %f \n",delta);
	 	 fprintf(stdout,"   p0 = %f \n",p0 / log_10);
	  }

	  if(exp(p0) > delta){
	 	 return exp(p0);
	  }

			double log_delta = log(delta);

	 	 //cout << "p0 = " << p0 << " with pi = " << pi << endl;

	 	 // generate all possible contingency tables with n1 fixed
	 	 if(DEBUG_B){fprintf(stderr,"Starting generating cont tables\n");
	   fprintf(stderr,"   delta = %f \n",log_delta);
	   fprintf(stderr,"   p0 = %f \n",p0);}


	 	 double p_table = 0.0;
	 	 double p_table_precomp = 0.0;

	 	 int continue_loop = 1;
		 int cache_index = 0;
		 int explored_tables = 0;

	 	 int xs = 1;
	 	 int a = 1;
	 	 int maxAS_ = 1;
	 	 int minAS_ = 1;
		 if(DEBUG_B){fprintf(stderr," x_s = %d , a_s = %d \n",x_s,a_s);
		 fprintf(stderr," epsilon_transactions = %d \n",epsilon_transactions);}
	 	 // confidence intervals on xs using epsilon
	 	 int xlb_ = max(testable_low , x_s - epsilon_transactions);
	 	 int xub_ = min(testable_high , x_s + epsilon_transactions);
		 if(DEBUG_B){fprintf(stderr," xlb_ = %d \n",xlb_);
		 fprintf(stderr," xub_ = %d \n",xub_);}
	 	 for (xs = xlb_; xs <= xub_; xs++){
	 		 maxAS_ = min(xs,n1);
	 		 minAS_ = max(0,xs-N_minus_n1);
			 if(DEBUG_B){fprintf(stderr," xs = %d \n",xs);
			 fprintf(stderr," N_n1 = %d \n",N_minus_n1);
			 fprintf(stderr," xs-N_n1 = %d \n",xs-N_minus_n1);
			 fprintf(stderr," minAS_ = %d \n",minAS_);
			 fprintf(stderr," maxAS_ = %d \n",maxAS_);}
	 		 // confidence intervals on as using epsilon
	 		 minAS_ = max(minAS_ , a_s - epsilon_transactions);
	 		 maxAS_ = min(maxAS_ , a_s + epsilon_transactions);
			 if(DEBUG_B){fprintf(stderr," minAS_ = %d \n",minAS_);
			 fprintf(stderr," maxAS_ = %d \n",maxAS_);}
	 		 p_table_precomp = (xs*log(pi)) + ((N-xs)*log(1-pi));
	 		 for (a = minAS_; a <= maxAS_; a++){
	 			 p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
				 explored_tables++;
	 			 if(p_table <= p0){
	 				 p_value = sumlogs(p_value , p_table);
	 				probabilities_cache[cache_index] = p_table;
	 				cache_index++;
	 				 /*if(DEBUG_B){
	 					 fprintf(stdout,"   x = %d , a = %d \n",xs,a);
	 					 fprintf(stdout,"   ptable = %f \n",p_table);
	 					 fprintf(stdout,"   pvalue increased = %f \n",p_value);
	 				 }*/
	 				 if(exp(p_value) > delta){
	 					 p_value = exp(p_value);
	 					 return p_value;
	 				 }
	 			 }
				 /*if(last_p_table > 0.0 && p_table < last_p_table && (p_value - p_table) > 50.0 && p_value > -1000.0){
					 a = maxAS_ + 1;
					 fprintf(stderr,"stopped computing ptables %f\n",p_table);
						fprintf(stderr,"last_p_table %f\n",last_p_table);
						fprintf(stderr,"p_value %f\n",p_value);
						fprintf(stderr,"delta %f\n",delta);
				 }*/
				 if(DEBUG_B){fprintf(stderr," x = %d , a = %d \n",xs,a);
					fprintf(stderr,"   p_table %f\n",p_table);
					fprintf(stderr,"   p0 %f\n",p0);
					fprintf(stderr,"   p_value %f\n",p_value);
					fprintf(stderr,"   delta %f\n",delta);}

	 		 }
	 	 }

	 	 //sort_(cache_index);
		 //print_cache(cache_index);
		 //double test = sum_(cache_index);
	 	 //fprintf(stderr,"   sum of probabilities = %f \n",p_value);
	 	 //fprintf(stderr,"   test = %f \n",test);

		 // normalization factor
		 //computeT(x_s , a_s , epsilon);
		 int index = (a_s * (N-n)) + (x_s-a_s);
		 p_value = p_value - Tvalues[index];
	 	 // convert to linear
	 	 p_value = exp(p_value);
		 fprintf(stderr,"   cache_index (v2) %d\n",cache_index);
		 fprintf(stderr,"   explored_tables (v2) %d\n",explored_tables);

	 	 if(DEBUG_B){
	 		 fprintf(stdout,"   final pvalue = %f \n",p_value);
	 	 }

	  return p_value;

	 }


	double min_ptable(int x , int a , double ub_ci_pi , double lb_ci_pi){

		double pi_ub = min(ub_ci_pi , 1.0);
		double pi_lb = max(0.0 , lb_ci_pi);
		double p_table_precomp_ub = min( (x*log(pi_ub)) + ((N-x)*log(1-pi_ub)) , (x*log(pi_lb)) + ((N-x)*log(1-pi_lb)) );
		double p_table_ub = p_table_precomp_ub + logbincoeff(N-n,x-a) + logbincoeff(n,a);
		return p_table_ub;

	}


	  // unconditional test version 3: test statistic is given by p-value of observed contingency table,
		// we find contingency tables such that their probability is > p0

	void computeUncondUBPValue_version3(int x_s , int a_s , double epsilon , double *pvalue_lowerbound , double *pvalue_upperbound){

		int DEBUG_L = 0;
		if(DEBUG_L){fprintf(stderr,"Improved Barnard 3: (x_s = %d , a_1 = %d , a_0 = %d)\n",x_s,a_s,x_s-a_s);}

		explored_patterns++;

		if(check_bounds(x_s , a_s , epsilon , pvalue_lowerbound , pvalue_upperbound) == 1){
			return;
		}

		// check simple upper bound
	  double pi = (double)x_s / (double)N;
		double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi));
	  int n1 = n;
		p0 = p0 + logbincoeff(N_minus_n1,x_s-a_s) + logbincoeff(n1,a_s);


		// if here, we need to explore contingency tables

	  int DEBUG_B = 0;
		int epsilon_transactions = (int)(epsilon * (double)N);
		#ifndef USECONFINTERVALS
		epsilon_transactions = N;
		#endif

		explored_contingency_tables++;

	  double p_value = -10000; // 0 in log scale
		double p_value_ub = -10000; // 0 in log scale

	   if(DEBUG_B){
	  	 fprintf(stderr,"Starting Improved Barnard 3: (x_s = %d , a_s = %d , N = %d, n1 = %d)\n",x_s,a_s,N,n);
	  	 fprintf(stderr,"   pi = %f \n",pi);
	   }

	   if(DEBUG_B){
	  	 fprintf(stderr,"   delta = %f \n",delta);
	  	 fprintf(stderr,"   p0 = %f \n",p0 / log_10);
	   }

	 		double log_delta = log(delta);

	  	 // generate contingency tables with n1 fixed
			 if(DEBUG_B){
	  	fprintf(stderr,"Starting generating cont tables\n");
	    fprintf(stderr,"   delta = %f \n",log_delta);
	    fprintf(stderr,"   p0 = %f \n",p0);
		}


			double p_table = 0.0;
			double p_table_ub = 0.0;
			double p_table_precomp = 0.0;
			double p_table_precomp_ub = 0.0;
			double a0 = x_s - a_s;
 	 		double a1 = a_s;
 	 		double n0 = N-n;
			double n1_ = n;
 	 		double f1 = a1 / n1;
 	 		double f0 = a0 / n0;
 			double epsilon1 = epsilon * (n0+n1_) / n1;
 			double epsilon0 = epsilon * (n0+n1_) / n0;
 	 		double ilb = (f0 < f1) ? f0 + epsilon0 : f1 + epsilon1;
 			double iub = (f0 > f1) ? f0 - epsilon0 : f1 - epsilon1;
 			double lb_ci_pi = iub;
 			double ub_ci_pi = ilb;
			double pi_ub = min(ub_ci_pi , 1.0);
			double pi_lb = max(0.0 , lb_ci_pi);
			if(DEBUG_L){fprintf(stderr,"   pi = %f \n",pi);
			fprintf(stderr,"   pi_ub = %f \n",pi_ub);
			fprintf(stderr,"   p0 = %f \n",p0);
			fprintf(stderr,"   p0min = %f \n",(x_s*log(pi_ub)) + ((N-x_s)*log(1-pi_ub))+logbincoeff(N_minus_n1,x_s-a_s) + logbincoeff(n1,a_s));}
			int maxAS_ = 0;
			int minAS_ = 0;
	  	int xlb_ = max(1 , x_s - epsilon_transactions);
	  	int xub_ = min(N-1 , x_s + epsilon_transactions);
			double n1_over_N = (double)n1 / (double)N;
			int cache_index = 0;
			int explored_tables = 0;

			if(DEBUG_B){fprintf(stderr," Moving x up\n");}
			int xs = x_s;
			int a = floor((double)xs * n1_over_N);
			//lower bound
			p_table_precomp = (xs*log(pi)) + ((N-xs)*log(1-pi));
			p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
			//lower bound
			// upper bound
			//pi_ub = xs/(double)N;
			p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
			// upper bound
			if(DEBUG_B){fprintf(stderr," x = %d\n",xs);
			fprintf(stderr,"  a = %d\n",a);
			fprintf(stderr,"   p0 = %f \n",p0 / log_10);
			fprintf(stderr,"   p_table = %f \n",p_table / log_10);
			fprintf(stderr,"   p_value = %f \n",p_value / log_10);}
			while(xs <= xub_ && p_table > p0){
				maxAS_ = min(xs,n1);
		  	minAS_ = max(0,xs-N_minus_n1);
				// consider confidence intervals
	  		minAS_ = max(minAS_ , a_s - epsilon_transactions);
	  		maxAS_ = min(maxAS_ , a_s + epsilon_transactions);
				if(DEBUG_B){
					fprintf(stderr," x = %d\n",xs);
					fprintf(stderr," Moving a up\n");
				}
				while(a <= maxAS_ && p_table > p0){
					explored_tables++;
					if(p_table_ub > p_table){
						fprintf(stderr,"Problem!  p_table_ub = %f , p_table = %f , x = %d , a = %d\n",p_table_ub,p_table,xs,a);
						p_table_precomp_ub = min( (xs*log(pi_ub)) + ((N-xs)*log(1-pi_ub)) , (xs*log(pi_lb)) + ((N-xs)*log(1-pi_lb)) );
					}
					p_value = sumlogs(p_value , p_table);
					p_value_ub = sumlogs(p_value_ub , p_table_ub);
					if(DEBUG_B){fprintf(stderr,"   sum cont. table with x = %d , a0 = %d , a1 = %d , ptable(pi=xs/n) = %f , ptable(pi=x/n) = %f , ",xs,xs-a,a,p_table,p_table_ub);
					fprintf(stderr,"pvalue_lb = %f , pvalue_ub = %f \n",1.0-exp(p_value),1.0-exp(p_value_ub));}
					if(DEBUG_B){fprintf(stderr,"  a = %d\n",a);
					fprintf(stderr,"   p0 = %f \n",p0 / log_10);
					fprintf(stderr,"   p_table = %f \n",p_table / log_10);
					fprintf(stderr,"   p_value = %f \n",p_value / log_10);}
					a = a + 1;
					p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
					// upper bound
					p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
					// upper bound
				}


				//double test_tail = incbeta(10, 10, 0.3);



				if(DEBUG_B){
				fprintf(stderr," Moving a down\n");}
				a = floor((double)xs * n1_over_N) - 1;
				p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
				// upper bound
				p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
				// upper bound
				while(a >= minAS_ && p_table > p0){
					explored_tables++;
					if(p_table_ub > p_table){
						fprintf(stderr,"Problem!  p_table_ub = %f , p_table = %f , x = %d , a = %d\n",p_table_ub,p_table,xs,a);
						p_table_precomp_ub = min( (xs*log(pi_ub)) + ((N-xs)*log(1-pi_ub)) , (xs*log(pi_lb)) + ((N-xs)*log(1-pi_lb)) );
					}
					p_value = sumlogs(p_value , p_table);
					p_value_ub = sumlogs(p_value_ub , p_table_ub);
					if(DEBUG_B){fprintf(stderr,"   sum cont. table with x = %d , a0 = %d , a1 = %d , ptable(pi=xs/n) = %f , ptable(pi=x/n) = %f , ",xs,xs-a,a,p_table,p_table_ub);
					fprintf(stderr,"pvalue_lb = %f , pvalue_ub = %f \n",1.0-exp(p_value),1.0-exp(p_value_ub));}
					if(DEBUG_B){fprintf(stderr,"  a = %d\n",a);
					fprintf(stderr,"   p0 = %f \n",p0 / log_10);
					fprintf(stderr,"   p_table = %f \n",p_table / log_10);
					fprintf(stderr,"   p_value = %f \n",p_value / log_10);}
					a = a - 1;
					p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
					// upper bound
					p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
					// upper bound
				}
				xs = xs + 1;
				a = floor((double)xs * n1_over_N);
				p_table_precomp = (xs*log(pi)) + ((N-xs)*log(1-pi));
				p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
				// upper bound
				p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
				// upper bound
			}
			if(DEBUG_B){fprintf(stderr," Moving x down\n");}
			xs = x_s - 1;
			a = floor((double)xs * n1_over_N);
			p_table_precomp = (xs*log(pi)) + ((N-xs)*log(1-pi));
			p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
			// upper bound
			//pi_ub = xs/(double)N;
			p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
			// upper bound
			if(DEBUG_B){fprintf(stderr," x = %d\n",xs);
			fprintf(stderr,"  a = %d\n",a);
			fprintf(stderr,"   p0 = %f \n",p0 / log_10);
			fprintf(stderr,"   p_table = %f \n",p_table / log_10);
			fprintf(stderr,"   p_value = %f \n",p_value / log_10);}
			while(xs >= xlb_ && p_table > p0){
				maxAS_ = min(xs,n1);
		  	minAS_ = max(0,xs-N_minus_n1);
				// consider confidence intervals
	  		minAS_ = max(minAS_ , a_s - epsilon_transactions);
	  		maxAS_ = min(maxAS_ , a_s + epsilon_transactions);
				if(DEBUG_B){fprintf(stderr," x = %d\n",xs);
				fprintf(stderr," Moving a up\n");}
				while(a <= maxAS_ && p_table > p0){
					explored_tables++;
					if(p_table_ub > p_table){
						fprintf(stderr,"Problem!  p_table_ub = %f , p_table = %f , x = %d , a = %d\n",p_table_ub,p_table,xs,a);
						p_table_precomp_ub = min( (xs*log(pi_ub)) + ((N-xs)*log(1-pi_ub)) , (xs*log(pi_lb)) + ((N-xs)*log(1-pi_lb)) );
					}
					p_value = sumlogs(p_value , p_table);
					p_value_ub = sumlogs(p_value_ub , p_table_ub);
					if(DEBUG_B){fprintf(stderr,"   sum cont. table with x = %d , a0 = %d , a1 = %d , ptable(pi=xs/n) = %f , ptable(pi=x/n) = %f , ",xs,xs-a,a,p_table,p_table_ub);
					fprintf(stderr,"pvalue_lb = %f , pvalue_ub = %f \n",1.0-exp(p_value),1.0-exp(p_value_ub));}
					if(DEBUG_B){fprintf(stderr,"  a = %d\n",a);
					fprintf(stderr,"   p0 = %f \n",p0 / log_10);
					fprintf(stderr,"   p_table = %f \n",p_table / log_10);
					fprintf(stderr,"   p_value = %f \n",p_value / log_10);}
					a = a + 1;
					p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
					// upper bound
					p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
					// upper bound
				}
				if(DEBUG_B){fprintf(stderr," Moving a down\n");}
				a = floor((double)xs * n1_over_N) - 1;
				p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
				// upper bound
				p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
				// upper bound
				while(a >= minAS_ && p_table > p0){
					explored_tables++;
					if(p_table_ub > p_table){
						fprintf(stderr,"Problem!  p_table_ub = %f , p_table = %f , x = %d , a = %d\n",p_table_ub,p_table,xs,a);
						p_table_precomp_ub = min( (xs*log(pi_ub)) + ((N-xs)*log(1-pi_ub)) , (xs*log(pi_lb)) + ((N-xs)*log(1-pi_lb)) );
					}
					p_value = sumlogs(p_value , p_table);
					p_value_ub = sumlogs(p_value_ub , p_table_ub);
					if(DEBUG_B){fprintf(stderr,"   sum cont. table with x = %d , a0 = %d , a1 = %d , ptable(pi=xs/n) = %f , ptable(pi=x/n) = %f , ",xs,xs-a,a,p_table,p_table_ub);
					fprintf(stderr,"pvalue_lb = %f , pvalue_ub = %f \n",1.0-exp(p_value),1.0-exp(p_value_ub));}
					if(DEBUG_B){fprintf(stderr,"  a = %d\n",a);
					fprintf(stderr,"   p0 = %f \n",p0 / log_10);
					fprintf(stderr,"   p_table = %f \n",p_table / log_10);
					fprintf(stderr,"   p_value = %f \n",p_value / log_10);}
					a = a - 1;
					p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
					// upper bound
					p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
					// upper bound
				}
				xs = xs - 1;
				a = floor((double)xs * n1_over_N);
				p_table_precomp = (xs*log(pi)) + ((N-xs)*log(1-pi));
				p_table = p_table_precomp + logbincoeff(N_minus_n1,xs-a) + logbincoeff(n1,a);
				// upper bound
				//pi_ub = xs/(double)N;
				p_table_ub = max(p0 , min_ptable(xs , a , ub_ci_pi , lb_ci_pi));
				// upper bound
			}


		 //sort_(cache_index);
		 //print_cache(cache_index);
		 //double test = sum_(cache_index);
		 //fprintf(stderr,"   sum of probabilities (v2) = %f \n",p_value);
		 //fprintf(stderr,"   test (v2)= %f \n",test);


		 #ifdef USECONFINTERVALS
			// normalization factor
			computeT(x_s , a_s , epsilon);
	 	 	int index = (a_s * (N-n)) + (x_s-a_s);
			fprintf(stderr,"  negC (log10) = %f \n",p_value / log_10);
			fprintf(stderr,"  negC = %f\n",exp(p_value));
			fprintf(stderr,"   Tvalues[index] (log10) = %f \n",Tvalues[index] / log_10);
			fprintf(stderr,"   Tvalues[index] = %f \n",exp(Tvalues[index]));
			fprintf(stderr,"   norm. negC (*) = %f \n",exp(p_value) / exp(Tvalues[index]));
			 p_value = p_value - Tvalues[index];
	 		fprintf(stderr,"   norm. negC (log10) = %f \n",p_value / log_10);
			fprintf(stderr,"   norm. negC  = %f \n",exp(p_value));
			#endif
			 // convert to linear and complement
			 double p_value_toreturn = 1.0 - max(0.0 , min(1.0 , exp(p_value)));
			 if(DEBUG_L){fprintf(stderr,"   sum of ub probabilities = %e \n",p_value_ub);}
			 double p_value_ub_toreturn = 1.0 - max(0.0 , min(1.0 , exp(p_value_ub)));
			 //if(DEBUG_B){//} && p_value_toreturn <= delta){
			 //fprintf(stderr,"   p_value = %f \n",p_value / log_10);
			 if(DEBUG_L){fprintf(stderr,"   norm. negC = %f \n",exp(p_value));
			 fprintf(stderr,"   p_value = %e \n",p_value_toreturn);
			 fprintf(stderr,"   p_value_ub = %e \n",p_value_ub_toreturn);
		 		//}
				fprintf(stderr,"  cache_index = %d \n",cache_index);
		 	 fprintf(stderr,"   explored_tables %d\n",explored_tables);}

	  	 if(DEBUG_B){
	  		 fprintf(stderr,"   final pvalue = %f \n",p_value_toreturn);
	  	 }


	 		 if(DEBUG_L){if(p_value_ub_toreturn <= delta){
	 		 	fprintf(stderr,"Rejected with UB!\n");}}

			 if(p_value_toreturn <= delta){
				 count_number_of_enumerate_rejects++;
			 	if(DEBUG_L){fprintf(stderr,"Rejected!\n");
				fprintf(stderr,"Rejected with p-value (%e) but not with UB2!\n",p_value_toreturn);}
				/*test_upper_bound = p0 + log((double)n1) + log((double)(N-n1)) - log((double)explored_tables);
		 	 	test_upper_bound = exp(test_upper_bound);
				if(DEBUG_L){fprintf(stderr,"Improved UB2  (%e) gives %d!\n",test_upper_bound,(int)(test_upper_bound<=delta));}
				if(test_upper_bound<=delta){
					return test_upper_bound;
				}*/
			}
			else{
				if(DEBUG_L){fprintf(stderr,"Not Rejected! %e (exp)\n",(p_value_toreturn-delta));}
		 		count_number_of_enumerate_nonrejects++;
				if(probabilities_thr[x_s] > p0){
					probabilities_thr[x_s] = p0;
				}
			}

			if(use_pvalues_cache == 1){
				unsigned long pcache_index = ((unsigned long)a_s * (unsigned long)(N-n)) + (x_s-a_s);
				pvalues_cache[pcache_index] = p_value_toreturn;
				pvalues_cache_flag[pcache_index] = 1;
			}

		*pvalue_lowerbound = p_value_toreturn;
		return;

	  }


		double computeP0(int x_s , int a_s){

			int n1 = n;

			double pi = (double)x_s / (double)N;

			double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi));
			p0 = p0 + logbincoeff(N-n1,x_s-a_s) + logbincoeff(n1,a_s);

			return exp(p0);

		}

		double compute_tails_a1_fixed_naive(int a1 , double p0 , double pi , int *numtables){

			int DEBUG_B = 0;
			int DEBUG_TAILS = 0;


			 if(DEBUG_B){
				fprintf(stderr,"call to compute_tails_a1_fixed_naive: (a1 = %d , p0 = %f )\n",a1,p0);
			 }

			int a0 = 0;
			int n1 = n;
			int n0 = N-n;
			int max_a0 = n0;
			int min_a0 = 0;
			double prob_a1 = (a1*log(pi)) + ((n1-a1)*log(1-pi)) + logbincoeff(n1,a1);
			double bin_tail = -10000;
			double bin_tail2 = -10000;

			int numtables_ = 0;

			if(DEBUG_B){
				fprintf(stderr,"  min_a0 = %d , max_a0 = %d\n",min_a0,max_a0);
			}

			double p_table;
			for(a0 = 0; a0 <= n0; a0++){
				p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
				if(p_table <= p0){
					p_table = p_table - prob_a1;
					bin_tail = sumlogs(bin_tail , p_table);
					numtables_++;
				}
				else{
					p_table = p_table - prob_a1;
					bin_tail2 = sumlogs(bin_tail2 , p_table);
				}
			}

			if(DEBUG_TAILS){
				fprintf(stderr,"  bin_tail2 = %e , bin_tail = %e\n",1.0-exp(bin_tail2),exp(bin_tail));
			}

			*numtables += (n0+1-numtables_);

			return exp(bin_tail);

		}


		double compute_tail_a1_fixed(int a1 , double p0 , int side , double pi , int *hint_a0){

			int DEBUG_B = 0;
			int DEBUG_TAILS = 0;


			 if(DEBUG_B){
				fprintf(stderr,"call to compute_tail_a1_fixed: (a1 = %d , p0 = %f , side = %d)\n",a1,p0,side);
			 }

			int max_a0 = 0;
			int min_a0 = 0;
			int a0 = 0;
			int n1 = n;
			int n0 = N-n;
			double prob_a1 = (a1*log(pi)) + ((n1-a1)*log(1-pi)) + logbincoeff(n1,a1);
			double bin_tail = -10000;
			double p_table;

			if(*hint_a0 < 0){

					if(side < 0){
						max_a0 = (int)(((double)(n0+1) * pi) + 0.5);
						min_a0 = 0;
						a0 = max_a0;
					}
					else{
							max_a0 = n0;
							min_a0 = (int)(((double)(n0+1) * pi) + 0.5) - 1;
							a0 = min_a0;
					}
					if(DEBUG_B){
						fprintf(stderr,"  min_a0 = %d , max_a0 = %d\n",min_a0,max_a0);
					}

					p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);

					if(DEBUG_B){
						fprintf(stderr,"  a0 = %d , p_table = %f , p0 = %f\n",a0,p_table,p0);
						fprintf(stderr,"  p_table > p0 = %d \n",p_table > p0);
					}

					if(p_table > p0){
						// binary search on a0
						if(DEBUG_B){
							fprintf(stderr,"  started binary search with side = %d\n",side);
						}
						while(max_a0 - min_a0 > 3){
							a0 = (int)((double)(max_a0 + min_a0) / 2.0);
							p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
							if(side > 0){
								if(p_table > p0){min_a0 = a0;}
								else{max_a0 = a0;}
							}
							else{
								if(p_table > p0){max_a0 = a0;}
								else{min_a0 = a0;}
							}
							if(DEBUG_B){
								fprintf(stderr,"    min_a0 = %d , max_a0 = %d , a0 = %d , p_table = %f , p0 = %f\n",min_a0,max_a0,a0,p_table,p0);
							}
						}

						double p_table_min = prob_a1 + (min_a0*log(pi)) + ((n0-min_a0)*log(1-pi)) + logbincoeff(n0,min_a0);
						double p_table_max = prob_a1 + (max_a0*log(pi)) + ((n0-max_a0)*log(1-pi)) + logbincoeff(n0,max_a0);
						double p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
						if(DEBUG_B){fprintf(stderr,"  Bin search finished with values min_a0 = %d , a0 = %d , max_a0 = %d , \n",min_a0,a0,max_a0);
						fprintf(stderr,"    prob values: p_table(min_a0) = %f , p_table(a0) = %f , p_table(max_a0) = %f\n",p_table_min,p_table,p_table_max);}

						// base case of binary search
						if(side < 0){
							a0 = max_a0;
							p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
							if(p_table > p0){
								while(a0 >= min_a0 && p_table > p0){
									a0--;
									p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
									if(DEBUG_B){
										fprintf(stderr,"    min_a0 = %d , max_a0 = %d , a0 = %d , p_table = %f , p0 = %f\n",min_a0,max_a0,a0,p_table,p0);
									}
								}
							}
						}
						else{
							a0 = min_a0;
							p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
							while(a0 <= max_a0 && p_table > p0){
								a0++;
								p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
								if(DEBUG_B){
									fprintf(stderr,"    min_a0 = %d , max_a0 = %d , a0 = %d , p_table = %f , p0 = %f\n",min_a0,max_a0,a0,p_table,p0);
								}
							}
						}

					}

					*hint_a0 = a0;
					if(DEBUG_B){fprintf(stderr,"  hint_a0 set to %d , side = %d \n",*hint_a0,side);}

				}
				// we already have an hint on the value of a0 to find
				else{
					a0 = *hint_a0;
					p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
					if(p_table > p0){
						while(a0 >= 0 && a0 <= n0 && p_table > p0){
							a0 += side;
							if(DEBUG_B){fprintf(stderr,"    hint_a0 update to %d , side = %d , p_table = %f , p0 = %f \n",a0,side,p_table,p0);}
							p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
						}
					}
					else{
						while(a0 >= 0 && a0 <= n0 && p_table <= p0){
							a0 -= side;
							if(DEBUG_B){fprintf(stderr,"    hint_a0 update to %d , side = %d , p_table = %f , p0 = %f \n",a0,side,p_table,p0);}
							p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
						}
						if(a0 >= 0 && a0 <= n0){
							a0 += side;
						}
					}
					*hint_a0 = a0;
					if(DEBUG_B){fprintf(stderr,"  FINISHED updating hint_a0 = %d , side = %d , p_table = %f , p0 = %f \n",*hint_a0,side,p_table,p0);}
				}
			// we need to compute the tail of binomial distribution of variables a0 and parameter pi
			// we use the incomplete regularized beta function
			// Pr(X >= a0) when X ~ Bin(n0 , pi)
			if(DEBUG_B){fprintf(stderr,"  min_a0 = %d , max_a0 = %d , a0 of tail = %d \n",min_a0,max_a0,a0);}
			if(side > 0){
				bin_tail = incbeta((double)a0, (double)(n0-a0+1), pi);
			}
			else{
				bin_tail = incbeta((double)(n0-a0+1) , (double)(a0+1), 1.0-pi);
			}
			//bin_tail = max(bin_tail , TINY);
			//bin_tail = min(bin_tail , 1.0-TINY);
			if(DEBUG_B){
				fprintf(stderr,"  bin_tail = %e , a0 = %d \n",bin_tail,a0);
			}

			// compute naively the tail
			if(DEBUG_TAILS){
				int i = 0;
				double tail2 = -10000;
				double p_table;
				if(side < 0){max_a0 = a0;min_a0 = 0;}
				else{max_a0 = n0;min_a0 = a0;}
				for(i=min_a0 ; i<=max_a0 ; i++){
					p_table = (i*log(pi)) + ((n0-i)*log(1-pi)) + logbincoeff(n0,i);
					tail2 = sumlogs(tail2 , p_table);
				}
				fprintf(stderr,"  bin_tail = %e , bin_tail2 debug = %e \n",bin_tail,exp(tail2));
			}

			return bin_tail;

		}




						void unconditional_fast_test2(int x_s , int a_s , double epsilon , double *pvalue_lowerbound , double *pvalue_upperbound){
							int DEBUG_B = 0;
							int DEBUG_L = 0;

							 if(DEBUG_B){
								fprintf(stdout,"Starting test vesion of unconditional: (x_s = %d , a_s = %d , N = %d, n1 = %d)\n",x_s,a_s,N,n);
							 }

							explored_patterns++;


								if(check_bounds(x_s , a_s , epsilon , pvalue_lowerbound , pvalue_upperbound) == 1){
									return;
								}

						  double pi = (double)x_s / (double)N;
						  int n1 = n;
						  int n0 = N-n;
							double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi))+ logbincoeff(n0,x_s-a_s) + logbincoeff(n1,a_s);


								if(DEBUG_B){fprintf(stderr,"Starting p-value computation\n");}

								double p_table = 0.0;
								double p_value = -10000;
								int xs = 1;
								int a0 = 1;
								int a1 = 1;
								int maxAS_ = 1;
								int minAS_ = 1;
								int max_a0 = 1;
								int min_a0 = 1;
								double prob_a1 = 0.0;
								double bin_tail = 0.0;
								double n1_over_N = (double)n1 /(double)N;
								double naive_tail;
								int numtables = 0;
								double max_diff_tails = 0.0;

								a1 = (int)(((double)(n1+1) * pi) + 0.5);

								for(a1=0; a1<=n1; a1++){
									prob_a1 = (a1*log(pi)) + ((n1-a1)*log(1-pi)) + logbincoeff(n1,a1);
									bin_tail = compute_tails_a1_fixed_naive(a1 , p0 , pi , &numtables);
									if(bin_tail < 1.0){
										bin_tail = log(1.0 - bin_tail);
										bin_tail = bin_tail + prob_a1;
										p_value = sumlogs(p_value , bin_tail);
									}
								}
								p_value = 1.0 - exp(p_value);

								if(DEBUG_B){fprintf(stderr,"p_value (fast2) = %e \n",p_value);}
								if(DEBUG_B){fprintf(stderr,"numtables (fast2) = %d \n",numtables);}

								*pvalue_lowerbound = p_value;
								*pvalue_upperbound = p_value;

						}




		void unconditional_fast_test_old(int x_s , int a_s , double epsilon , double *pvalue_lowerbound , double *pvalue_upperbound){
			int DEBUG_B = 0;
			int DEBUG_L = 0;
			int DEBUG_TAILS = 0;

			 if(DEBUG_B){
				fprintf(stdout,"Starting test vesion of unconditional: (x_s = %d , a_s = %d , N = %d, n1 = %d)\n",x_s,a_s,N,n);
			 }

			explored_patterns++;


				if(check_bounds(x_s , a_s , epsilon , pvalue_lowerbound , pvalue_upperbound) == 1){
					return;
				}

		  double pi = (double)x_s / (double)N;
		  int n1 = n;
		  int n0 = N-n;
			double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi))+ logbincoeff(n0,x_s-a_s) + logbincoeff(n1,a_s);


				if(DEBUG_B){fprintf(stderr,"Starting p-value computation\n");}

				double p_table = 0.0;
				double p_value = -10000;
				int xs = 1;
				int a0 = 1;
				int a1 = 1;
				int maxAS_ = 1;
				int minAS_ = 1;
				int max_a0 = 1;
				int min_a0 = 1;
				double prob_a1 = 0.0;
				double bin_tail = 0.0;
				double n1_over_N = (double)n1 /(double)N;
				double naive_tail;
				int numtables = 0;
				double max_diff_tails = 0.0;
				int a0_hint_side1 = -1;
				int a0_hint_sidem1 = -1;

				a1 = (int)(((double)(n1+1) * pi) + 0.5);
				prob_a1 = (a1*log(pi)) + ((n1-a1)*log(1-pi)) + logbincoeff(n1,a1);
				a0 = (int)(((double)(n0+1) * pi) + 0.5);
				xs = a0 + a1;
				p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
				if(DEBUG_B){fprintf(stderr,"a1 = %d , a0 = %d , p_table = %f\n",a1,a0,p_table);}
				while(a1 <= n1 && p_table > p0){
					// compute first tail
					bin_tail = compute_tail_a1_fixed(a1 , p0 , 1 , pi , &a0_hint_side1);
					if(DEBUG_B){fprintf(stderr,"first tail = %e \n",bin_tail);}
					// add second tail
					bin_tail = bin_tail + compute_tail_a1_fixed(a1 , p0 , -1 , pi , &a0_hint_sidem1);
					if(DEBUG_TAILS){
						//debug
						naive_tail = compute_tails_a1_fixed_naive(a1 , p0 , pi , &numtables);
						if(naive_tail/bin_tail >= 1.001 || naive_tail/bin_tail <= 0.999){
							fprintf(stderr,"sum of two tails = %e , naive tail = %e\n",bin_tail,naive_tail);
							max_diff_tails = max(max_diff_tails , abs(naive_tail - bin_tail));
							fprintf(stderr,"ERROR WITH TAILS! sum of two tails = %e , naive tail = %e\n",bin_tail,naive_tail);
						}
					}
					fprintf(stderr,"bin_tail = %e , 1.0 - bin_tail = %e , log(1.0 - bin_tail) = %e\n",bin_tail,1.0 - bin_tail,log(1.0 - bin_tail));
					bin_tail = log(1.0 - bin_tail);
					if(DEBUG_B){fprintf(stderr,"prob_a1 = %f\n",prob_a1);
					fprintf(stderr,"bin_tail = %e , test = %e\n",bin_tail,log(1.0 - exp(bin_tail)));}
					// increase the pvalue by the log-complement of the sum of the tails
					if(DEBUG_B){fprintf(stderr,"p_value before = %f\n",p_value);}
					bin_tail = bin_tail + prob_a1;
					if(DEBUG_B){fprintf(stderr,"summing = %f\n",bin_tail);}
					p_value = sumlogs(p_value , bin_tail);
					if(DEBUG_B){fprintf(stderr,"p_value after = %f\n",p_value);}
					// we increase a1 of 1
					a1++;
					if(DEBUG_B){fprintf(stderr,"a1 increased to a1 = %d\n",a1);}
					if(a1 <= n1){
						prob_a1 = (a1*log(pi)) + ((n1-a1)*log(1-pi)) + logbincoeff(n1,a1);
						a0 = (int)(((double)(n0+1) * pi) + 0.5);
						p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
						if(DEBUG_B){fprintf(stderr,"a1 = %d , a0 = %d , p_table = %f\n",a1,a0,p_table);}
					}
				}
				if(DEBUG_B){fprintf(stderr,"STOPPED increasing a1. a1 = %d , p_table = %f , p0 = %f\n",a1,p_table,p0);
					if(a1 >= 0 && a1 <=n1){
						for(a0=0; a0<=n0; a0++){
							p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
							if(p_table > p0){
								fprintf(stderr,"ERROR! I found someting > p0 at a0 = %d , p_table = %f , p0 = %f\n",a0,p_table,p0);
							}
						}
					}
				}
				// the same decreasing a1
				a1 = (int)(((double)(n1+1) * pi) + 0.5)-1;
				prob_a1 = (a1*log(pi)) + ((n1-a1)*log(1-pi)) + logbincoeff(n1,a1);
				a0 = (int)(((double)(n0+1) * pi) + 0.5);
				p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
				if(DEBUG_B){fprintf(stderr,"a1 = %d , a0 = %d , p_table = %f\n",a1,a0,p_table);}
				while(a1 >= 0 && p_table > p0){
					// compute first tail
					bin_tail = compute_tail_a1_fixed(a1 , p0 , 1 , pi , &a0_hint_side1);
					if(DEBUG_B){fprintf(stderr,"first tail = %e \n",bin_tail);}
					// add second tail
					bin_tail = bin_tail + compute_tail_a1_fixed(a1 , p0 , -1 , pi , &a0_hint_sidem1);
					if(DEBUG_TAILS){
						//debug
						naive_tail = compute_tails_a1_fixed_naive(a1 , p0 , pi , &numtables);
						if(naive_tail/bin_tail >= 1.0001 || naive_tail/bin_tail <= 0.9999){
							fprintf(stderr,"sum of two tails = %e , naive tail = %e\n",bin_tail,naive_tail);
							max_diff_tails = max(max_diff_tails , abs(naive_tail - bin_tail));
							fprintf(stderr,"ERROR WITH TAILS! sum of two tails = %e , naive tail = %e\n",bin_tail,naive_tail);
						}
					}
					fprintf(stderr,"bin_tail = %e , 1.0 - bin_tail = %e , log(1.0 - bin_tail) = %e\n",bin_tail,1.0 - bin_tail,log(1.0 - bin_tail));
					bin_tail = log(1.0 - bin_tail);
					if(DEBUG_B){fprintf(stderr,"prob_a1 = %f\n",prob_a1);
					fprintf(stderr,"bin_tail = %e , test = %e\n",bin_tail,log(1.0 - exp(bin_tail)));}
					// increase the pvalue by the log-complement of the sum of the tails
					if(DEBUG_B){fprintf(stderr,"p_value before = %f\n",p_value);}
					bin_tail = bin_tail + prob_a1;
					if(DEBUG_B){fprintf(stderr,"summing = %f\n",bin_tail);}
					p_value = sumlogs(p_value , bin_tail);
					if(DEBUG_B){fprintf(stderr,"p_value after = %f\n",p_value);}
					// we decrease a1 of 1
					a1--;
					if(DEBUG_B){fprintf(stderr,"a1 decreased to a1 = %d\n",a1);}
					if(a1 >= 0){
						a0 = (int)(((double)(n0+1) * pi) + 0.5);
						prob_a1 = (a1*log(pi)) + ((n1-a1)*log(1-pi)) + logbincoeff(n1,a1);
						p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
						if(DEBUG_B){fprintf(stderr,"a1 = %d , a0 = %d , p_table = %f\n",a1,a0,p_table);}
					}
				}
				if(DEBUG_B){fprintf(stderr,"STOPPED decreasing a1. a1 = %d , p_table = %f , p0 = %f\n",a1,p_table,p0);
					if(a1 >= 0 && a1 <=n1){
						for(a0=0; a0<=n0; a0++){
							p_table = prob_a1 + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
							if(p_table > p0){
								fprintf(stderr,"ERROR! I found someting > p0 at a0 = %d , p_table = %f , p0 = %f\n",a0,p_table,p0);
							}
						}
					}
				}
				p_value = 1.0 - exp(p_value);
				if(DEBUG_B){fprintf(stderr,"FAST numtables = %d \n",numtables);}
				if(DEBUG_TAILS){fprintf(stderr,"FAST max_diff_tails = %e \n",max_diff_tails);}
				if(DEBUG_B){fprintf(stderr,"p_value (fast) = %e \n",p_value);}

				*pvalue_lowerbound = p_value;
				*pvalue_upperbound = p_value;

				//unconditional_fast_test2(x_s , a_s , epsilon , pvalue_lowerbound , pvalue_upperbound);

		}


		void unconditional_fast_test(int x_s , int a_s , double epsilon , double *pvalue_lowerbound , double *pvalue_upperbound){
			int DEBUG_B = 0;
			int DEBUG_L = 0;
			int DEBUG_TAILS = 0;

			explored_patterns++;


				if(check_bounds(x_s , a_s , epsilon , pvalue_lowerbound , pvalue_upperbound) == 1){
					return;
				}

			#ifdef VERSION4UBONLY
			return;
			#endif

			if(DEBUG_B){
				fprintf(stderr,"Starting new vesion of unconditional: (x_s = %d , a_s = %d , N = %d, n1 = %d)\n",x_s,a_s,N,n);
			}

			explored_contingency_tables++;

			double pi = (double)x_s / (double)N;
			double log_pi = log(pi);
			double log_1_pi = log(1.0 - pi);
			int n1 = n;
			int n0 = N-n;
			double p0 = (x_s*log_pi) + ((N-x_s)*log_1_pi)+ logbincoeff(n0,x_s-a_s) + logbincoeff(n1,a_s);


				if(DEBUG_B){fprintf(stderr,"Starting p-value computation\n");}

				double p_table = 0.0;
				double p_value = -10000;
				int xs = 1;
				int a0 = 1;
				int a1 = 1;
				int maxAS_ = 1;
				int minAS_ = 1;
				int max_a0 = 1;
				int min_a0 = 1;
				double prob_a1 = 0.0;
				double bin_tail = 0.0;
				double n1_over_N = (double)n1 /(double)N;
				double naive_tail;
				int numtables = 0;
				double max_diff_tails = 0.0;
				int a0_hint_side1 = -1;
				int a0_hint_sidem1 = -1;

				a1 = (int)(((double)(n1+1) * pi) + 0.5);
				prob_a1 = (a1*log_pi) + ((n1-a1)*log_1_pi) + logbincoeff(n1,a1);
				a0 = (int)(((double)(n0+1) * pi) + 0.5);
				xs = a0 + a1;
				p_table = prob_a1 + (a0*log_pi) + ((n0-a0)*log_1_pi) + logbincoeff(n0,a0);
				if(DEBUG_B){fprintf(stderr,"a1 = %d , a0 = %d , p_table = %f\n",a1,a0,p_table);}
				while(a1 <= n1 && p_table > p0){
					// compute first tail
					bin_tail = compute_tail_a1_fixed(a1 , p0 , 1 , pi , &a0_hint_side1);
					if(DEBUG_B){fprintf(stderr,"first tail = %e \n",bin_tail);}
					// add second tail
					bin_tail = bin_tail + compute_tail_a1_fixed(a1 , p0 , -1 , pi , &a0_hint_sidem1);
					if(DEBUG_TAILS){
						//debug
						naive_tail = compute_tails_a1_fixed_naive(a1 , p0 , pi , &numtables);
						if(naive_tail/bin_tail >= 1.001 || naive_tail/bin_tail <= 0.999){
							fprintf(stderr,"sum of two tails = %e , naive tail = %e\n",bin_tail,naive_tail);
							max_diff_tails = max(max_diff_tails , abs(naive_tail - bin_tail));
							fprintf(stderr,"ERROR WITH TAILS! sum of two tails = %e , naive tail = %e\n",bin_tail,naive_tail);
						}
					}
					//fprintf(stderr,"bin_tail = %e , 1.0 - bin_tail = %e , log(1.0 - bin_tail) = %e\n",bin_tail,1.0 - bin_tail,log(1.0 - bin_tail));
					bin_tail = log(bin_tail);
					if(DEBUG_B){fprintf(stderr,"prob_a1 = %f\n",prob_a1);
					fprintf(stderr,"bin_tail = %e , test = %e\n",bin_tail,log(1.0 - exp(bin_tail)));}
					// increase the pvalue by the log-complement of the sum of the tails
					if(DEBUG_B){fprintf(stderr,"p_value before = %f\n",p_value);}
					bin_tail = bin_tail + prob_a1;
					if(DEBUG_B){fprintf(stderr,"summing = %f\n",bin_tail);}
					p_value = sumlogs(p_value , bin_tail);
					if(DEBUG_B){fprintf(stderr,"p_value after = %f\n",p_value);}
					// we increase a1 of 1
					a1++;
					if(DEBUG_B){fprintf(stderr,"a1 increased to a1 = %d\n",a1);}
					if(a1 <= n1){
						prob_a1 = (a1*log_pi) + ((n1-a1)*log_1_pi) + logbincoeff(n1,a1);
						a0 = (int)(((double)(n0+1) * pi) + 0.5);
						p_table = prob_a1 + (a0*log_pi) + ((n0-a0)*log_1_pi) + logbincoeff(n0,a0);
						if(DEBUG_B){fprintf(stderr,"a1 = %d , a0 = %d , p_table = %f\n",a1,a0,p_table);}
					}
				}
				if(DEBUG_B){fprintf(stderr,"STOPPED increasing a1. a1 = %d , p_table = %f , p0 = %f\n",a1,p_table,p0);
					if(a1 >= 0 && a1 <=n1){
						for(a0=0; a0<=n0; a0++){
							p_table = prob_a1 + (a0*log_pi) + ((n0-a0)*log_1_pi) + logbincoeff(n0,a0);
							if(p_table > p0){
								fprintf(stderr,"ERROR! I found someting > p0 at a0 = %d , p_table = %f , p0 = %f\n",a0,p_table,p0);
							}
						}
					}
				}
				double right_tail_a1 = 0.0;
				if(a1 >= 0 && a1 <=n1){
					right_tail_a1 = incbeta((double)a1, (double)(n1-a1+1), pi);
				}
				// the same decreasing a1
				a1 = (int)(((double)(n1+1) * pi) + 0.5)-1;
				prob_a1 = (a1*log_pi) + ((n1-a1)*log_1_pi) + logbincoeff(n1,a1);
				a0 = (int)(((double)(n0+1) * pi) + 0.5);
				p_table = prob_a1 + (a0*log_pi) + ((n0-a0)*log_1_pi) + logbincoeff(n0,a0);
				if(DEBUG_B){fprintf(stderr,"a1 = %d , a0 = %d , p_table = %f\n",a1,a0,p_table);}
				while(a1 >= 0 && p_table > p0){
					// compute first tail
					bin_tail = compute_tail_a1_fixed(a1 , p0 , 1 , pi , &a0_hint_side1);
					if(DEBUG_B){fprintf(stderr,"first tail = %e \n",bin_tail);}
					// add second tail
					bin_tail = bin_tail + compute_tail_a1_fixed(a1 , p0 , -1 , pi , &a0_hint_sidem1);
					if(DEBUG_TAILS){
						//debug
						naive_tail = compute_tails_a1_fixed_naive(a1 , p0 , pi , &numtables);
						if(naive_tail/bin_tail >= 1.0001 || naive_tail/bin_tail <= 0.9999){
							fprintf(stderr,"sum of two tails = %e , naive tail = %e\n",bin_tail,naive_tail);
							max_diff_tails = max(max_diff_tails , abs(naive_tail - bin_tail));
							fprintf(stderr,"ERROR WITH TAILS! sum of two tails = %e , naive tail = %e\n",bin_tail,naive_tail);
						}
					}
					//fprintf(stderr,"bin_tail = %e , 1.0 - bin_tail = %e , log(1.0 - bin_tail) = %e\n",bin_tail,1.0 - bin_tail,log(1.0 - bin_tail));
					bin_tail = log(bin_tail);
					if(DEBUG_B){fprintf(stderr,"prob_a1 = %f\n",prob_a1);
					fprintf(stderr,"bin_tail = %e , test = %e\n",bin_tail,log(1.0 - exp(bin_tail)));}
					// increase the pvalue by the log-complement of the sum of the tails
					if(DEBUG_B){fprintf(stderr,"p_value before = %f\n",p_value);}
					bin_tail = bin_tail + prob_a1;
					if(DEBUG_B){fprintf(stderr,"summing = %f\n",bin_tail);}
					p_value = sumlogs(p_value , bin_tail);
					if(DEBUG_B){fprintf(stderr,"p_value after = %f\n",p_value);}
					// we decrease a1 of 1
					a1--;
					if(DEBUG_B){fprintf(stderr,"a1 decreased to a1 = %d\n",a1);}
					if(a1 >= 0){
						a0 = (int)(((double)(n0+1) * pi) + 0.5);
						prob_a1 = (a1*log_pi) + ((n1-a1)*log_1_pi) + logbincoeff(n1,a1);
						p_table = prob_a1 + (a0*log_pi) + ((n0-a0)*log_1_pi) + logbincoeff(n0,a0);
						if(DEBUG_B){fprintf(stderr,"a1 = %d , a0 = %d , p_table = %f\n",a1,a0,p_table);}
					}
				}
				if(DEBUG_B){fprintf(stderr,"STOPPED decreasing a1. a1 = %d , p_table = %f , p0 = %f\n",a1,p_table,p0);
					if(a1 >= 0 && a1 <=n1){
						for(a0=0; a0<=n0; a0++){
							p_table = prob_a1 + (a0*log_pi) + ((n0-a0)*log_1_pi) + logbincoeff(n0,a0);
							if(p_table > p0){
								fprintf(stderr,"ERROR! I found someting > p0 at a0 = %d , p_table = %f , p0 = %f\n",a0,p_table,p0);
							}
						}
					}
				}
				double left_tail_a1 = 0.0;
				if(a1 >= 0 && a1 <=n1){
					left_tail_a1 = incbeta((double)(n1-a1+1) , (double)a1, 1.0-pi);
				}
				p_value = exp(p_value) + right_tail_a1 + left_tail_a1;
				if(DEBUG_B){fprintf(stderr,"FAST numtables = %d \n",numtables);}
				if(DEBUG_TAILS){fprintf(stderr,"FAST max_diff_tails = %e \n",max_diff_tails);}
				if(DEBUG_B){fprintf(stderr,"p_value (fast new) = %e \n",p_value);}




 			 if(p_value <= delta){
 				 count_number_of_enumerate_rejects++;
 			 	if(DEBUG_L){fprintf(stderr,"Rejected!\n");
 				fprintf(stderr,"Rejected with p-value (%e) but not with UB2!\n",p_value);}
 				/*test_upper_bound = p0 + log((double)n1) + log((double)(N-n1)) - log((double)explored_tables);
 		 	 	test_upper_bound = exp(test_upper_bound);
 				if(DEBUG_L){fprintf(stderr,"Improved UB2  (%e) gives %d!\n",test_upper_bound,(int)(test_upper_bound<=delta));}
 				if(test_upper_bound<=delta){
 					return test_upper_bound;
 				}*/
 			}
 			else{
 				if(DEBUG_L){fprintf(stderr,"Not Rejected! %e (exp)\n",(p_value-delta));}
 		 		count_number_of_enumerate_nonrejects++;
 				if(probabilities_thr[x_s] > p0){
 					probabilities_thr[x_s] = p0;
 				}
 			}

 			if(use_pvalues_cache == 1){
 				unsigned long pcache_index = ((unsigned long)a_s * (unsigned long)(N-n)) + (x_s-a_s);
 				pvalues_cache[pcache_index] = p_value;
 				pvalues_cache_flag[pcache_index] = 1;
 			}

 			*pvalue_lowerbound = p_value;
			//*pvalue_upperbound = p_value;

		}






void compute_naive_ut(int x_s , int a_s , double epsilon , double *pvalue_lowerbound , double *pvalue_upperbound){

		/*if(check_bounds(x_s , a_s , epsilon , pvalue_lowerbound , pvalue_upperbound) == 1){
			return;
		}*/

		int DEBUG_L = 0;

		explored_patterns++;
		explored_contingency_tables++;

		double p_value = -100000;
		double p_value2 = -100000;
		double p_value3 = 0;
		double p_value4 = 0;

		int a0 , a1 , xs;
		int n1 = n;
		int n0 = N-n1;
	  double pi = (double)x_s / (double)N;
		double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi)) + logbincoeff(n0,x_s-a_s) + logbincoeff(n1,a_s);
		double p_table;
		double a1_prob;
		int max_a1 = 0;
		int min_a1 = n1;
		int numtables = 0;

		for(a1 = 0; a1 <= n1 ; a1++){
			a1_prob = (a1*log(pi)) + ((n1-a1)*log(1-pi)) + logbincoeff(n1,a1);
			for(a0 = 0; a0 <= n0 ; a0++){
				p_table = a1_prob + (a0*log(pi)) + ((n0-a0)*log(1-pi)) + logbincoeff(n0,a0);
				if(p_table > p0){
					p_value = sumlogs(p_value , p_table);
					min_a1 = min(min_a1 , a1);
					max_a1 = max(max_a1 , a1);
					numtables++;
					p_value4 = p_value4 + exp(p_table);
				}
				else{
					p_value2 = sumlogs(p_value2 , p_table);
					p_value3 = p_value3 + exp(p_table);
				}
			}
		}

		p_value = 1.0 - exp(p_value);
		p_value2 = exp(p_value2);
		p_value4 = 1.0 - p_value4;

		if(DEBUG_L){fprintf(stderr,"NAIVE: pvalue = %e , numtables = %d , min_a1 = %d , max_a1 = %d\n",p_value,numtables,min_a1,max_a1);
		fprintf(stderr,"NAIVE: pvalue = %e , p_value2 = %e , p_value3 = %e , p_value4 = %e\n",p_value,p_value2,p_value3,p_value4);}

		if(p_value <= delta){
			count_number_of_enumerate_rejects++;
		}
		else{
			count_number_of_enumerate_nonrejects++;
		}

		*pvalue_lowerbound = p_value;
		*pvalue_upperbound = p_value;
		return;

	  }


		int find_interval_limit(int x_s , int a_s , int side , int left , int *tosum){

			int DEBUG_B = 0;
			int maxa , mina , a , tota;
			double prob_thr, p_table;
			double pi = (double)x_s / (double)N;
			double log_pi = log(pi);
			double log_1_pi = log(1.0-pi);

			// this is to prevent using n0 and n1 again in the rest of the function
			// we will use tota, maxa and mina variables
			{
				int n1 = n;
				int n0 = N-n;
				prob_thr = (x_s*log_pi) + ((N-x_s)*log_1_pi) + logbincoeff(n1,a_s) + logbincoeff(n0,(x_s-a_s));
				if(side > 0){
					a = (int)(((double)(n1+1) * pi) + 0.5);
					prob_thr -= (a*log_pi) + ((n1-a)*log_1_pi) + logbincoeff(n1,a);
					tota = n0;
				}
				else{
					a = (int)(((double)(n0+1) * pi) + 0.5);
					prob_thr -= (a*log_pi) + ((n0-a)*log_1_pi) + logbincoeff(n0,a);
					tota = n1;
				}
			}

			if(left > 0){
				maxa = (int)(((double)(tota+1) * pi) + 0.5)-1;
				mina = 0;
			}
			else{
				maxa = (int)(((double)(tota+1) * pi) + 0.5);
				mina = tota;
			}
			if(DEBUG_B){
				fprintf(stderr,"  mina = %d , maxa = %d\n",mina,maxa);
			}

			a = mina;
			p_table = (a*log_pi) + ((tota-a)*log_1_pi) + logbincoeff(tota,a);

			if(DEBUG_B){
				fprintf(stderr,"  a = %d , p_table = %f , prob_thr = %f\n",a,p_table,prob_thr);
				fprintf(stderr,"  p_table > prob_thr = %d \n",p_table > prob_thr);
			}

			if(p_table < prob_thr){
				// binary search on a0
				if(DEBUG_B){
					fprintf(stderr,"  started binary search with side = %d\n",side);
				}
				while(abs(maxa - mina) > 3){
					a = (int)((double)(maxa + mina) / 2.0);
					p_table = (a*log_pi) + ((tota-a)*log_1_pi) + logbincoeff(tota,a);
					if(p_table > prob_thr){
						maxa = a;
					}
					else{
						mina = a;
					}
					if(DEBUG_B){
						fprintf(stderr,"    mina = %d , maxa = %d , a = %d , p_table = %f , prob_thr = %f\n",mina,maxa,a,p_table,prob_thr);
					}
				}

				// base case of binary search
				a = maxa;
				p_table = (a*log_pi) + ((tota-a)*log_1_pi) + logbincoeff(tota,a);
				int lb_a = min(mina , maxa);
				int ub_a = max(mina , maxa);
				int direction = (mina > maxa) ? +1 : -1;
				if(p_table > prob_thr){
					while(a >= lb_a && a <= ub_a && p_table > prob_thr){
						a+=direction;
						p_table = (a*log_pi) + ((tota-a)*log_1_pi) + logbincoeff(tota,a);
						if(DEBUG_B){
							fprintf(stderr,"    mina = %d , maxa = %d , a = %d , p_table = %f , prob_thr = %f\n",mina,maxa,a,p_table,prob_thr);
						}
					}
				}
			}
			p_table = (a*log_pi) + ((tota-a)*log_1_pi) + logbincoeff(tota,a);
			*tosum = p_table <= prob_thr;
			if(DEBUG_B){
				double p_table_min = (mina*log_pi) + ((tota-mina)*log_1_pi) + logbincoeff(tota,mina);
				double p_table_max = (maxa*log_pi) + ((tota-maxa)*log_1_pi) + logbincoeff(tota,maxa);
				double p_table = (a*log_pi) + ((tota-a)*log_1_pi) + logbincoeff(tota,a);
				fprintf(stderr,"  Bin search finished with values mina = %d , a = %d , maxa = %d , \n",mina,a,maxa);
				fprintf(stderr,"    prob values: p_table(mina) = %f , p_table(a) = %f , p_table(maxa) = %f\n",p_table_min,p_table,p_table_max);
				if(*tosum == 0) fprintf(stderr,"  NOT SUMMING!\n");
				fprintf(stderr,"  final value of a = %d , p_table = %f , prob_thr = %f\n",a,p_table,prob_thr);
			}
			return a;

		}


		double fast_lower_bound_box(int x_s , int a_s){

			int DEBUG_B = 0;
			if(DEBUG_B){
			 fprintf(stderr,"Starting fast_lower_bound_box: (x_s = %d , a_s = %d)\n",x_s,a_s);
			}

			int n1 = n;
			int n0 = N-n;
			int tosum;
			double pi = (double)x_s / (double)N;
			double lower_bound_1 = 0.0;
			double lower_bound_0 = 0.0;
			int a1_left = find_interval_limit(x_s , a_s , 0 , 1 , &tosum);
			if(tosum == 1){lower_bound_1 += incbeta((double)(n1-a1_left+1) , (double)(a1_left+1), 1.0-pi);}
			fprintf(stderr," tail1 %e\n",lower_bound_1);
			int a1_right = find_interval_limit(x_s , a_s , 0 , 0 , &tosum);
			if(tosum == 1){lower_bound_1 += incbeta((double)a1_right, (double)(n1-a1_right+1), pi);}
			fprintf(stderr," tail2 %e\n",lower_bound_1);
			int a0_left = find_interval_limit(x_s , a_s , 1 , 1 , &tosum);
			if(tosum == 1){lower_bound_0 += incbeta((double)(n0-a0_left+1) , (double)(a0_left+1), 1.0-pi);}
			fprintf(stderr," tail3 %e\n",lower_bound_0);
			int a0_right = find_interval_limit(x_s , a_s , 1 , 0 , &tosum);
			if(tosum == 1){lower_bound_0 += incbeta((double)a0_right, (double)(n0-a0_right+1), pi);}
			fprintf(stderr," tail4 %e\n",lower_bound_0);

			return lower_bound_0 * lower_bound_1;

		}
