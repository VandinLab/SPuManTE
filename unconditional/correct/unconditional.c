// source code for unconditional test

// this enables version 3
//#define VERSION3 1
//#define USECONFINTERVALS 1
// this enables version 4
#define VERSION4 1
//#define VERSION4UBONLY 1
//#define DVERSION4 1
//#define VERSIONNAIVE 1

#include "incbeta.c"


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

		  if (x_s == N || x_s == 1){
		 		*pvalue_lowerbound = 1.0;
		 		*pvalue_upperbound = 1.0;
				return 1;
			}

			int DEBUG_L = 0;
			if(DEBUG_L){fprintf(stderr,"check_bounds: (x_s = %d , a_1 = %d , a_0 = %d)\n",x_s,a_s,x_s-a_s);}

			// check simple upper bound
		  double pi = (double)x_s / (double)N;
		  int n1 = n;


			// check confidence intervals
			double lb_ci_pi; double ub_ci_pi;
			if(confidenceIntervalsNHReject(x_s , a_s , epsilon , &lb_ci_pi , &ub_ci_pi)){
			 	*pvalue_lowerbound = 0.0;
			 	*pvalue_upperbound = 0.0;
			 	return 1;
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



		double computeP0(int x_s , int a_s){

			int n1 = n;

			double pi = (double)x_s / (double)N;

			double p0 = (x_s*log(pi)) + ((N-x_s)*log(1-pi));
			p0 = p0 + logbincoeff(N-n1,x_s-a_s) + logbincoeff(n1,a_s);

			return exp(p0);

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



		void unconditional_fast_test(int x_s , double epsilon , double *pvalue_lowerbound , double *pvalue_upperbound){
			int DEBUG_B = 0;
			int DEBUG_L = 0;
			int DEBUG_TAILS = 0;

			int a_s = min(x_s , n);

			if(DEBUG_B){
				fprintf(stderr,"Starting new vesion of unconditional: (x_s = %d , a_s = %d , N = %d, n1 = %d)\n",x_s,a_s,N,n);
			}

				if(check_bounds(x_s , a_s , epsilon , pvalue_lowerbound , pvalue_upperbound) == 1){
					if(DEBUG_B){fprintf(stderr,"bounds 1\n");}
					return;
				}
				if(DEBUG_B){fprintf(stderr,"bounds 0\n");}

			#ifdef VERSION4UBONLY
			return;
			#endif

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
					if(DEBUG_B){fprintf(stderr,"computing right tail with a1 %d \n",a1);}
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
					if(DEBUG_B){fprintf(stderr,"computing left tail with a1 %d \n",a1);}
					left_tail_a1 = incbeta((double)(n1-a1+1) , (double)a1, 1.0-pi);
				}
				p_value = exp(p_value) + right_tail_a1 + left_tail_a1;
				if(DEBUG_B){fprintf(stderr,"FAST numtables = %d \n",numtables);}
				if(DEBUG_TAILS){fprintf(stderr,"FAST max_diff_tails = %e \n",max_diff_tails);}
				if(DEBUG_B){fprintf(stderr,"p_value (fast new) = %e \n",p_value);}

				if(DEBUG_B){fprintf(stderr,"right_tail_a1 = %e \n",right_tail_a1);
				fprintf(stderr,"left_tail_a1 = %e \n",left_tail_a1);
				fprintf(stderr,"p_value (fast new) = %e \n",p_value);}


 			*pvalue_lowerbound = p_value;
			//*pvalue_upperbound = p_value;

		}




void compute_naive_ut(int x_s , double epsilon , double *pvalue_lowerbound , double *pvalue_upperbound){

		/*if(check_bounds(x_s , a_s , epsilon , pvalue_lowerbound , pvalue_upperbound) == 1){
			return;
		}*/

		int a_s = min(x_s , n);

		int DEBUG_L = 0;

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

		*pvalue_lowerbound = p_value;
		*pvalue_upperbound = p_value;
		return;

	  }
