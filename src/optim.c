// Compile with R CMD SHLIB optim.c multimin.c  -lgsl -lgslcblas -o  optim.so

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <assert.h>
#include "multimin.h"

const double sqrt_2pi = 2.506628274631000502415765;
int sample_size = 0, num_parents = 0;

double log1plusexp(double x){
	// Calculates log(1+exp(x)) and avoids overflow
	return (fmax(x,0)+log(1+exp(-fabs(x))));
}

double logNC_common(const double *g, const double *h, const double *k){
	// log normalizing constant
	// g, h, k: const double *, scalar parameters
	return (log1plusexp((*g)+(*h)*(*h)/2/(*k)-log(*k/2/M_PI)/2));
}

double logNC_multi(const int *n, const double *G, const double *H, const double *k){
	// mean log normalizing constant given rows of (G, H, k)
	// n: const int *
	// G, H: const double *, n-vectors of parameters
	// k: const double *, scalar parameter
	double tmp = 0;
	for (int i=0; i<*n; i++){
		tmp += logNC_common(G+i, H+i, k);
	}
	return (tmp/(*n));
}

double meanlogprob_common(const int *n, const double *WY, const double *g, const double *h, const double *k){
	// Mean log prob for a simple 1-d Hurdle, where each sample share the same scalar parameters
	// n: const int *, sample size
	// W: const double *, 2n-vector, W and Y concatenated: W n-vector of indicators, Y n-vector of data
	// g, h, k: const double *, scalar parameters
	double Wsum = 0, Ysum = 0, Y2sum = 0;
	for (int i=0; i<*n; i++){
		Wsum += WY[i];
		Ysum += WY[*n+i];
		Y2sum += WY[*n+i]*WY[*n+i];
	}
	return (((*g)*Wsum+(*h)*Ysum-(*k)*Y2sum/2)/(*n) - logNC_common(g, h, k));
}

double meanlogprob_multi(const int *n, double *WY, const double *G, const double *H, const double *k){
	// Mean log prob for a simple 1-d Hurdle, where each sample has different parameters (depend on the obs of their parents)
	// n: const int *, sample size
	// WY: const double *, 2n-vector, W and Y concatenated: W n-vector of indicators, Y n-vector of data
	// G, H: const double *, n-vectors of parameters
	// k: const double *, scalar parameters
	double tmp = 0, Y2sum = 0;
	for (int i=0; i<*n; i++){
		tmp += G[i]*WY[i]+H[i]*WY[*n+i];
		Y2sum += WY[*n+i]*WY[*n+i];
	}
	return ((tmp-(*k)*Y2sum/2)/(*n) - logNC_multi(n, G, H, k));
}

void sum_G(const int *n, const int *numpa, const double *par, const double *WYpa, double *G){
	// Calculates G[1:n] from parents, where G[i] = g+q'Wpa[i,]+r'Ypa[i,]
	// n: const int *, sample size
	// numpa: const int *, number of PARENTS
	// par: const double *, vector of canonical parameters of the conditional distribution of length 4*numpa+3: g q r h s t k
	// WYpa: const double *, matrix (n*(2*numpa)-vector) of W and Y for parents
	// G: double *, n-vector that stores the results
	for (int i=0; i<*n; i++){
		G[i] = par[0]; // g
		for (int j=0; j<*numpa; j++)
			G[i] += par[j+1]*WYpa[j*(*n)+i] + par[*numpa+j+1]*WYpa[(*numpa+j)*(*n)+i]; // q'Wpa[i,]+r'Ypa[i,]
	}
}

void sum_H(const int *n, const int *numpa, const double *par, const double *WYpa, double *H){
	// Calculates H[1:n] from parents, where H[i] = h+s'Wpa[i,]-t'Ypa[i,]
	// n: const int *, sample size
	// numpa: const int *, number of PARENTS
	// par: const double *, vector of canonical parameters of the conditional distribution of length 4*numpa+3: g q r h s t k
	// WYpa: const double *, matrix (n*(2*numpa)-vector) of W and Y for parents
	// H: double *, n-vector that stores the results
	par += 2*(*numpa)+1; // to skip g q r
	for (int i=0; i<*n; i++){
		H[i] = par[0]; // h
		for (int j=0; j<*numpa; j++)
			H[i] += par[j+1]*WYpa[j*(*n)+i] - par[*numpa+j+1]*WYpa[(*numpa+j)*(*n)+i]; // s'Wpa[i,]-t'Ypa[i,]; note the minus sign
	}
}

void negloglik(const size_t par_size, const double *par, void *fparams, double *fval){
	assert(par_size == 4*num_parents+3);
	double* WY =  (double *)fparams; // [W Y Wpa Ypa] = (sample_size*(2+2*num_parents))
	double* G = (double*)malloc(sample_size*sizeof(double)), *H = (double*)malloc(sample_size*sizeof(double));
	sum_G(&sample_size, &num_parents, par, WY+2*sample_size, G);
	sum_H(&sample_size, &num_parents, par, WY+2*sample_size, H);
	*fval = -(meanlogprob_multi(&sample_size, WY, G, H, par+4*num_parents+2));
	free(G); free(H);
}

void grad_one(const double w, const double y, const double g, const double h, const double k, double *gr3){
	double exp_part = exp(g+h*h/2/k);
	gr3[0] = (w-1)+1/(1+exp_part*sqrt_2pi/sqrt(k));
	gr3[1] = y-h/k+h/(k+exp_part*sqrt_2pi*sqrt(k));
	gr3[2] = (h*h/k+1)*(0.5/k-1/(2*sqrt(k)*(sqrt(k)+exp_part*sqrt_2pi)))-y*y/2;
}

void grad_full(const int *n, const int *numpa, double *WY, const double *G, const double *H, const double *k, double *grad){
	double* grad_tmp = (double*)malloc(3*sizeof(double));
	for (int j=0; j<4*(*numpa)+3; j++){
		grad[j] = 0;
	}
	for (int i=0; i<*n; i++){
		grad_one(WY[i], WY[i+*n], G[i], H[i], *k, grad_tmp);
		grad[0] += grad_tmp[0]; // g
		grad[2*(*numpa)+1] += grad_tmp[1]; // h
		grad[4*(*numpa)+2] += grad_tmp[2]; // k
		for (int j=0; j<*numpa; j++){
			grad[1+j] += grad_tmp[0]*WY[(2+j)*(*n)+i]; // q
			grad[*numpa+1+j] += grad_tmp[0]*WY[(*numpa+2+j)*(*n)+i]; // r
			grad[2*(*numpa)+2+j] += grad_tmp[1]*WY[(2+j)*(*n)+i]; // s
			grad[3*(*numpa)+2+j] -= grad_tmp[1]*WY[(*numpa+2+j)*(*n)+i]; // t; note the minus sign
		}
	}
	for (int j=0; j<4*(*numpa)+3; j++){
		grad[j] /= -(*n);
	}
}

void grad_negloglik(const size_t par_size, const double *par, void *fparams, double *grad){
	assert(par_size == 4*num_parents+3);
	double* WY =  (double *)fparams;
	double* G = (double*)malloc(sample_size*sizeof(double)), *H = (double*)malloc(sample_size*sizeof(double));
	sum_G(&sample_size, &num_parents, par, WY+2*sample_size, G);
	sum_H(&sample_size, &num_parents, par, WY+2*sample_size, H);
	grad_full(&sample_size, &num_parents, WY, G, H, par+4*num_parents+2, grad);
	free(G); free(H);
}

void fdf_func(const size_t par_size, const double *par, void *fparams, double *fval, double *grad){
	assert(par_size == 4*num_parents+3);
	double* WY =  (double *)fparams; // [W Y Wpa Ypa] = (sample_size*(2+2*num_parents))
	double* G = (double*)malloc(sample_size*sizeof(double)), *H = (double*)malloc(sample_size*sizeof(double));
	sum_G(&sample_size, &num_parents, par, WY+2*sample_size, G);
	sum_H(&sample_size, &num_parents, par, WY+2*sample_size, H);
	*fval = -(meanlogprob_multi(&sample_size, WY, G, H, par+4*num_parents+2));
	grad_full(&sample_size, &num_parents, WY, G, H, par+4*num_parents+2, grad);
	free(G); free(H);
}


void optim(const int *sample_size_input, const int *num_parents_input, double *minimum, double *grad, double *par_init, const double *WY, double *step_size, double *tol, unsigned *maxiter, double *epsabs, double *maxsize, unsigned *method, unsigned *verbosity){
	sample_size = *sample_size_input;
	num_parents = *num_parents_input;
	unsigned par_size = 4*num_parents+3;
	unsigned *t = (unsigned*)malloc(par_size*sizeof(unsigned));
	double *xmin = (double*)malloc(par_size*sizeof(double));
	for (int j=0; j<par_size-1; j++){
		t[j] = 0; xmin[j] = -INFINITY;
	}
	t[par_size-1] = 4; xmin[par_size-1] = 0;
	struct multimin_params optim_par = {*step_size, *tol, *maxiter, *epsabs, *maxsize, *method, *verbosity};
	multimin(par_size,par_init,minimum,t,xmin,NULL,&negloglik,&grad_negloglik,&fdf_func,(void *) WY,optim_par);
	grad_negloglik(par_size, par_init, (void *)WY, grad);
	free(t); free(xmin);
}

