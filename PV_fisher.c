/* ===============================================================================*/
/*   Version 1.0.             Cullan Howlett                                      */
/*   Copyright (c) 2017       International Centre for Radio Astronomy Research,  */
/*   The MIT License (MIT)    University of Western Australia                     */
/*                                                                                */
/* Permission is hereby granted, free of charge, to any person obtaining a copy   */
/* of this software and associated documentation files (the "Software"), to deal  */
/* in the Software without restriction, including without limitation the rights   */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      */
/* copies of the Software, and to permit persons to whom the Software is          */
/* furnished to do so, subject to the following conditions:                       */
/*                                                                                */
/* The above copyright notice and this permission notice shall be included in     */
/* all copies or substantial portions of the Software.                            */
/*                                                                                */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN      */
/* THE SOFTWARE.                                                                  */
/* ===============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>

// Fisher matrix calculation for surveys with velocity and density field measurements.
// ASSUMPTIONS:
//  -  Uncorrelated shot-noise between the density and velocity fields
//  -  I use the trapezium rule to integrate over r. There will be some error due to this, but this makes the most sense as 
//     we are binning the number density anyway, and it makes hella difference in the speed of the code.
//  -  The redshift dependence of the non-linear matter and velocity divergence power spectra is captured using linear interpolation.
//  -  The PV error scales as a fixed pecentage of H0*r.
//  -  Flat LCDM cosmology (but not necessarily GR as gammaval can be changed).
//  -  The damping of the velocity and density fields due to non-linear RSD is redshift independent

// The parameters necessary for the calculation
static int nparams = 4;           // The number of free parameters (we can use any of beta, fsigma8, r_g, sigma_g, sigma_u)
static int Data[4] = {0,1,3,4};   // A vector of flags for the parameters we are interested in (0=beta, 1=fsigma8, 2=r_g, 3=sigma_g, 4=sigma_u). MAKE SURE THE LENGTH OF THIS VECTOR, NPARAMS AND THE ENTRIES AGREE/MAKE SENSE, OR YOU MIGHT GET NONSENSE RESULTS!!
static int nziter = 10;           // Now many bins in redshift between zmin and zmax we are considering
static double zmin = 0.0;         // The minimum redshift to consider (You must have power spectra that are within this range or GSL spline will error out)
static double zmax = 0.5;         // The maximum redshift to consider (You must have power spectra that are within this range or GSL spline will error out)
static double Om = 0.3089;        // The matter density at z=0
static double c = 299792.458;     // The speed of light in km/s
static double gammaval = 0.55;    // The value of gammaval to use in the forecasts (where f(z) = Om(z)^gammaval)
static double r_g = 1.0;          // The cross correlation coefficient between the velocity and density fields
static double beta0 = 0.437;      // The value of beta (at z=0, we'll modify this by the redshift dependent value of bias and f as required)
static double sigma80 = 0.8150;   // The value of sigma8 at z=0
static double sigma_u = 13.00;    // The value of the velocity damping parameter in Mpc/h. I use the values from Jun Koda's paper
static double sigma_g = 4.24;     // The value of the density damping parameter in Mpc/h. I use the values from Jun Koda's paper
static double kmax = 0.2;         // The maximum k to evaluate for dd, dv and vv correlations (Typical values are 0.1 - 0.2, on smaller scales the models are likely to break down).
static double survey_area[3] = {0.0, 0.0, 1.745};   // We need to know the survey area for each survey and the overlap area between the surveys (redshift survey only first, then PV survey only, then overlap. 
                                                    // For fully overlapping we would have {0, 0, size_overlap}. For redshift larger than PV, we would have {size_red-size_overlap, 0, size_overlap}). Units are pi steradians, such that full sky is 4.0, half sky is 2.0 etc.
static double error_rand = 300.0;    // The observational error due to random non-linear velocities (I normally use 300km/s as in Jun Koda's paper)
static double error_dist = 0.05;     // The percentage error on the distance indicator (Typically 0.05 - 0.10 for SNe IA, 0.2 or more for Tully-Fisher or Fundamental Plane) 
static double verbosity = 0;         // How much output to give: 0 = only percentage errors on fsigma8, 1 = other useful info and nuisance parameters, 2 = full fisher and covariance matrices

// The number of redshifts and the redshifts themselves of the input matter and velocity divergence power spectra. 
// These numbers are multiplied by 100, converted to ints and written in the form _z0p%02d which is then appended to the filename Pvel_file. See routine read_power. 
static double nzin = 11;
static double zin[11] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
char * Pvel_file = "./example_files/example_pk";                                                  // The file containing the velocity divergence power spectrum. Don't include .dat as we'll append the redshifts on read in

// The files containing the number density of the surveys. First is the PV survey, then the redshift survey. These files MUST have the same binning and redshift range, 
// so that the sum over redshift bins works (would be fine if we used splines), i.e., if one survey is shallower then that file must contain rows with n(z)=0.
char * nbar_file[300] = {"./example_files/example_nbar_vel.dat",
                         "./example_files/example_nbar_red.dat"};      

// Other global parameters and arrays
int NK, * NRED;
double pkkmin;         // The minimum kmin to integrate over, based on the input power spectrum file
double pkkmax;         // The maximum k in the input power spectrum. The maximum k to integrate over is the smallest of this or kmax
double * zarray;
double * rarray;
double * deltararray;
double * growtharray;
double ** nbararray;
double * karray, * deltakarray;
double ** pmmarray, ** pmtarray, ** pttarray;
gsl_spline * growth_spline, * r_spline;
gsl_interp_accel * growth_acc, * r_acc;

// Prototypes
double zeff_integrand(double mu, void * pin);
double mu_integrand(double mu, void * pin);
double ezinv(double x, void *p);
double rz(double red);
double growthfunc(double x, void *p);
double growthz(double red);
void read_nz();
void read_power();

// Calculates the fished matrix for a velocity survey.
int main(int argc, char **argv) {
    
    FILE * fout;
    int i, j;

    // Read in the velocity divergence power spectrum output from the COPTER code (Carlson 2009)
    read_power();

    // Read in the number densities of the surveys
    read_nz();

    // Run some checks
    if (!((survey_area[0] > 0.0) || (survey_area[2] > 0.0))) {
        for (i=0; i<nparams; i++) {
            if (Data[i] == 2) {
                printf("ERROR: r_g is a free parameter, but there is no information in the density field (Fisher matrix will be singular)\n");
                exit(0);
            }
            if (Data[i] == 3) {
                printf("ERROR: sigma_g is a free parameter, but there is no information in the density field (Fisher matrix will be singular)\n");
                exit(0);
            }
        }
    }
    if (!((survey_area[1] > 0.0) || (survey_area[2] > 0.0))) {
        for (i=0; i<nparams; i++) {
            if (Data[i] == 4) {
                printf("ERROR: sigma_u is a free parameter, but there is no information in the velocity field (Fisher matrix will be singular)\n");
                exit(0);
            }
        }
    }
    if ((sizeof(Data)/sizeof(*Data)) != nparams) {
        printf("ERROR: Size of Data vector for parameters of interest must be equal to nparams\n");
        exit(0);
    }

    // Calculate the Fisher matrices for all bins. 
    gsl_matrix * Fisher_Tot = gsl_matrix_alloc(nparams, nparams);
    for (i=0; i<nparams; i++) {
        for (j=0; j<nparams; j++) gsl_matrix_set(Fisher_Tot, i, j, 0.0);
    }

    printf("Evaluating the Fisher Matrix for %d bins between [z_min = %lf, z_max = %lf]\n", nziter, zmin, zmax);

    if (verbosity == 0) printf("#     zmin         zmax         zeff      fsigma8(z_eff)   percentage error(z_eff)\n");

    int ziter;
    for (ziter = 0; ziter<nziter; ziter++) {

        double zbinwidth = (zmax-zmin)/(nziter);
        double zmin_iter = ziter*zbinwidth + zmin;
        double zmax_iter = (ziter+1.0)*zbinwidth + zmin;

        double rzmax = gsl_spline_eval(r_spline, zmax_iter, r_acc);
        double kmin = M_PI/rzmax;

        if (verbosity > 0) printf("Evaluating the Fisher Matrix for [k_min = %lf, k_max = %lf] and [z_min = %lf, z_max = %lf]\n", kmin, kmax, zmin_iter, zmax_iter);

        // Calculate the effective redshift (which I base on the sum of the S/N for the density and velocity fields)
        int numk;
        double k_sum1 = 0.0, k_sum2 = 0.0;
        for (numk=0; numk<NK; numk++) {

            double k = karray[numk]+0.5*deltakarray[numk];
            double deltak = deltakarray[numk];
            if (k < kmin) continue;
            if (k > kmax) continue;

            double result, error;
            double params[4] = {numk, k, zmin_iter, zmax_iter};

            size_t nevals = 1000;
            gsl_function F;
            F.function = &zeff_integrand;
            F.params = &params;
            gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

            gsl_integration_qags(&F, 0.0, 1.0, 0, 5e-3, nevals, w, &result, &error);
            gsl_integration_workspace_free(w);

            k_sum1 += k*k*deltak*result;
            k_sum2 += k*k*deltak;

        }
        double z_eff = k_sum1/k_sum2;
        if (verbosity > 0) printf("Effective redshift z_eff = %lf\n", z_eff);

        double growth_eff = gsl_spline_eval(growth_spline, z_eff, growth_acc);

        // Calculate the fisher matrix, integrating over k, then mu, then r (r is last as it means we are effectively integrating over effective volume).
        // As the input spectra are tabulated we'll just use the trapezium rule to integrate over k
        gsl_matrix * Fisher = gsl_matrix_alloc(nparams, nparams);
        for (i=0; i<nparams; i++) {
            for (j=i; j<nparams; j++) {

                double k_sum = 0.0;
                for (numk=0; numk<NK; numk++) {

                    double k = karray[numk]+0.5*deltakarray[numk];
                    double deltak = deltakarray[numk];
                    if (k < kmin) continue;
                    if (k > kmax) continue;

                    double result, error;
                    double params[6] = {numk, k, Data[i], Data[j], zmin_iter, zmax_iter};

                    size_t nevals = 1000;
                    gsl_function F;
                    F.function = &mu_integrand;
                    F.params = &params;
                    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

                    gsl_integration_qags(&F, 0.0, 1.0, 0, 5e-3, nevals, w, &result, &error);
                    gsl_integration_workspace_free(w);

                    k_sum += k*k*deltak*result;

                }
                //printf("%d, %d, %lf\n", i, j, k_sum/(4.0*M_PI));
                gsl_matrix_set(Fisher, i, j, k_sum/(4.0*M_PI));
                gsl_matrix_set(Fisher, j, i, k_sum/(4.0*M_PI));

            }
        }

        for (i=0; i<nparams; i++) {
            for (j=0; j<nparams; j++) {
                double val = gsl_matrix_get(Fisher_Tot, i, j) + gsl_matrix_get(Fisher, i, j);
                gsl_matrix_set(Fisher_Tot, i, j, val);
            }
        }

        if (verbosity == 2) {
            printf("Fisher Matrix\n======================\n");
            for (i=0; i<nparams; i++) {
                printf("[");
                for (j=0; j<nparams; j++) printf("%15.6lf,\t", gsl_matrix_get(Fisher, i, j));
                printf("],\n");
            }
        }

        // Now invert the Fisher matrix
        int s;
        gsl_permutation * p;
        gsl_matrix * Covariance = gsl_matrix_alloc(nparams, nparams);
        p = gsl_permutation_alloc(nparams);
        gsl_linalg_LU_decomp(Fisher, p, &s);
        gsl_linalg_LU_invert(Fisher, p, Covariance);
        gsl_permutation_free(p);

        double sigma8 = sigma80 * growth_eff;
        double Omz = Om*ezinv(z_eff,NULL)*ezinv(z_eff,NULL)*(1.0+z_eff)*(1.0+z_eff)*(1.0+z_eff);
        double f = pow(Omz, gammaval);
        double beta = f*beta0*growth_eff/pow(Om,0.55);

        if (verbosity == 0) {
            for (i=0; i<nparams; i++) {
                if (Data[i] == 1) printf("%12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf\n", zmin_iter, zmax_iter, z_eff, f*sigma8, 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/(f*sigma8));
            }
        }

        if (verbosity > 0) {
            for (i=0; i<nparams; i++) {
                if (Data[i] == 0) {
                    printf("beta = %12.6lf +/- %12.6lf\n", beta, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on beta\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/beta);
                }
                if (Data[i] == 1) {
                    printf("fsigma8 = %12.6lf +/- %12.6lf\n", f*sigma8, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on fsigma8\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/(f*sigma8));
                }
                if (Data[i] == 2) {
                    printf("r_g = %12.6lf +/- %12.6lf\n", r_g, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on r_g\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/r_g);
                }
                if (Data[i] == 3) {
                    printf("sigma_g = %12.6lf +/- %12.6lf\n", sigma_g, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on sigma_g\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/sigma_g);
                }
                if (Data[i] == 4) {
                    printf("sigma_u = %12.6lf +/- %12.6lf\n", sigma_u, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on sigma_u\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/sigma_u);
                }
            }
        }

        if (verbosity == 2) {
            printf("Covariance Matrix\n======================\n");
            for (i=0; i<nparams; i++) {
                printf("[");
                for (j=0; j<nparams; j++) printf("%15.6lf,\t", gsl_matrix_get(Covariance, i, j));
                printf("],\n");
            }
        }

        gsl_matrix_free(Fisher);
        gsl_matrix_free(Covariance);
    }

    // Now the full Fisher matrix over all redshifts if we had more than 1 redshift bin
    if (nziter > 1) {
        double rzmax = gsl_spline_eval(r_spline, zmax, r_acc);
        double kmin = M_PI/rzmax;

        if (verbosity > 0) printf("Finally, evaluating the Fisher Matrix for [k_min = %lf, k_max = %lf] and [z_min = %lf, z_max = %lf]\n", kmin, kmax, zmin, zmax);

        // Calculate the effective redshift
        int numk;
        double k_sum1 = 0.0, k_sum2 = 0.0;
        for (numk=0; numk<NK; numk++) {

            double k = karray[numk]+0.5*deltakarray[numk];
            double deltak = deltakarray[numk];
            if (k < kmin) continue;
            if (k > kmax) continue;

            double result, error;
            double params[4] = {numk, k, zmin, zmax};

            size_t nevals = 1000;
            gsl_function F;
            F.function = &zeff_integrand;
            F.params = &params;
            gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

            gsl_integration_qags(&F, 0.0, 1.0, 0, 5e-3, nevals, w, &result, &error);
            gsl_integration_workspace_free(w);

            k_sum1 += k*k*deltak*result;
            k_sum2 += k*k*deltak;

        }
        double z_eff = k_sum1/k_sum2;
        if (verbosity > 0) printf("Effective redshift z_eff = %lf\n", z_eff);

        double growth_eff = gsl_spline_eval(growth_spline, z_eff, growth_acc);

        if (verbosity == 2) {
            printf("Fisher Matrix\n======================\n");
            for (i=0; i<nparams; i++) {
                printf("[");
                for (j=0; j<nparams; j++) printf("%15.6lf,\t", gsl_matrix_get(Fisher_Tot, i, j));
                printf("],\n");
            }
        }

        // Now invert the Fisher matrix
        int s;
        gsl_permutation * p;
        gsl_matrix * Covariance = gsl_matrix_alloc(nparams, nparams);
        p = gsl_permutation_alloc(nparams);
        gsl_linalg_LU_decomp(Fisher_Tot, p, &s);
        gsl_linalg_LU_invert(Fisher_Tot, p, Covariance);
        gsl_permutation_free(p);

        double sigma8 = sigma80 * growth_eff;
        double Omz = Om*ezinv(z_eff,NULL)*ezinv(z_eff,NULL)*(1.0+z_eff)*(1.0+z_eff)*(1.0+z_eff);
        double f = pow(Omz, gammaval);
        double beta = f*beta0*growth_eff/pow(Om,0.55);

        if (verbosity == 0) {
            printf("# Full redshift range:\n");
            for (i=0; i<nparams; i++) {
                if (Data[i] == 1) printf("%12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf\n", zmin, zmax, z_eff, f*sigma8, 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/(f*sigma8));
            }
        }

        if (verbosity > 0) {
            for (i=0; i<nparams; i++) {
                if (Data[i] == 0) {
                    printf("beta = %12.6lf +/- %12.6lf\n", beta, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on beta\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/beta);
                }
                if (Data[i] == 1) {
                    printf("fsigma8 = %12.6lf +/- %12.6lf\n", f*sigma8, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on fsigma8\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/(f*sigma8));
                }
                if (Data[i] == 2) {
                    printf("r_g = %12.6lf +/- %12.6lf\n", r_g, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on r_g\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/r_g);
                }
                if (Data[i] == 3) {
                    printf("sigma_g = %12.6lf +/- %12.6lf\n", sigma_g, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on sigma_g\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/sigma_g);
                }
                if (Data[i] == 4) {
                    printf("sigma_u = %12.6lf +/- %12.6lf\n", sigma_u, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on sigma_u\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/sigma_u);
                }
            }
        }

        if (verbosity == 2) {
            printf("Covariance Matrix\n======================\n");
            for (i=0; i<nparams; i++) {
                printf("[");
                for (j=0; j<nparams; j++) printf("%15.6lf,\t", gsl_matrix_get(Covariance, i, j));
                printf("],\n");
            }
        }
    }

    gsl_matrix_free(Fisher_Tot);
    gsl_spline_free(growth_spline);
    gsl_interp_accel_free(growth_acc);
   
    return 0;
}

// The integrand to calculate the effective redshift. I'm not actually sure how this is done in the case of 
// density and velocity field measurements, but it seems logical to base it on the sum of the integral of the density spectra and the velocity power spectra
// weighted by their effective signal to noise. In this way the effective redshift is calculated in the same way as a redshift survey, but there is some dependence
// on the S/N in the velocity power spectrum too. In any case, the S/N of the density field measurement is always much higher (because there are
// no errors and the number density is higher) and so this dominates the effective redshift calculation.
double zeff_integrand(double mu, void * pin) {

    int i, j, m, q, u, surv;
    double * p = (double *)pin;

    int numk = (int)p[0];
    double k = p[1];
    double zminval = p[2];
    double zmaxval = p[3];
    gsl_interp_accel * Pmm_acc, * Pmt_acc, * Ptt_acc;
    gsl_spline * Pmm_spline, * Pmt_spline, * Ptt_spline;
    if (nzin > 1) {
        double * Pmm_array = (double *)malloc(nzin*sizeof(double));
        double * Pmt_array = (double *)malloc(nzin*sizeof(double));
        double * Ptt_array = (double *)malloc(nzin*sizeof(double));
        for (j=0; j<nzin; j++) {
            Pmm_array[j] = pmmarray[j][numk];
            Pmt_array[j] = pmtarray[j][numk];
            Ptt_array[j] = pttarray[j][numk];
        }
        Pmm_acc    = gsl_interp_accel_alloc();
        Pmm_spline = gsl_spline_alloc(gsl_interp_cspline, nzin);
        gsl_spline_init(Pmm_spline, zin, Pmm_array, nzin);
        free(Pmm_array);

        Pmt_acc    = gsl_interp_accel_alloc();
        Pmt_spline = gsl_spline_alloc(gsl_interp_cspline, nzin);
        gsl_spline_init(Pmt_spline, zin, Pmt_array, nzin);
        free(Pmt_array);

        Ptt_acc    = gsl_interp_accel_alloc();
        Ptt_spline = gsl_spline_alloc(gsl_interp_cspline, nzin);
        gsl_spline_init(Ptt_spline, zin, Ptt_array, nzin);
        free(Ptt_array);
    }

    double dendamp = sqrt(1.0/(1.0+0.5*(k*k*mu*mu*sigma_g*sigma_g)));     // This is unitless
    double veldamp = sin(k*sigma_u)/(k*sigma_u);                          // This is unitless

    double dVeff = 0.0, zdVeff = 0.0;
    for (i=0; i<NRED[0]; i++) {

        double zval = zarray[i];
        if (zval < zminval) continue;
        if (zval > zmaxval) break;

        double r_sum = 0.0;
        double r = rarray[i];
        double deltar = deltararray[i];

        double dd_prefac=0.0, vv_prefac=0.0;
        double P_gg=0.0, P_uu=0.0;

        double sigma8 = sigma80 * growtharray[i];

        // First lets calculate the relevant power spectra. Interpolate the power spectra linearly in redshift
        double Pmm, Pmt, Ptt;
        Pmm = gsl_spline_eval(Pmm_spline, zval, Pmm_acc);
        Pmt = gsl_spline_eval(Pmt_spline, zval, Pmm_acc);
        Ptt = gsl_spline_eval(Ptt_spline, zval, Pmm_acc);

        double Omz = Om*ezinv(zval,NULL)*ezinv(zval,NULL)*(1.0+zval)*(1.0+zval)*(1.0+zval);
        double f = pow(Omz, gammaval);
        double beta = f*beta0*growtharray[i]/pow(Om,0.55);

        vv_prefac  = 1.0e2*f*mu*veldamp/k;
        dd_prefac = (1.0/(beta*beta) + 2.0*r_g*mu*mu/beta + mu*mu*mu*mu)*f*f*dendamp*dendamp;
        P_gg = dd_prefac*Pmm;
        P_uu = vv_prefac*vv_prefac*Ptt;

        // We need to do the overlapping and non-overlapping parts of the redshifts and PV surveys separately
        for (surv=0; surv<3; surv++) {
            double surv_sum = 0.0;
            if (survey_area[surv] > 0.0) {
                double error_obs, error_noise, n_g = 0.0, n_u = 0.0;

                // Set the nbar for each section.
                if (surv == 0) {
                    n_g = nbararray[1][i];
                } else if (surv == 1) {
                    error_obs = 100.0*error_dist*r;                              // Percentage error * distance * H0 in km/s (factor of 100.0 comes from hubble parameter)
                    error_noise = error_rand*error_rand + error_obs*error_obs;   // Error_noise is in km^{2}s^{-2}
                    n_u = nbararray[0][i]/error_noise;                   
                } else {
                    error_obs = 100.0*error_dist*r;                              // Percentage error * distance * H0 in km/s (factor of 100.0 comes from hubble parameter)
                    error_noise = error_rand*error_rand + error_obs*error_obs;   // Error_noise is in km^{2}s^{-2}
                    n_u = nbararray[0][i]/error_noise;                   
                    n_g = nbararray[1][i];
                }

                double value1 = n_g/(1.0 + n_g*P_gg);
                double value2 = n_u/(1.0 + n_u*P_uu);
                surv_sum += value1*value1 + value2*value2;

                surv_sum *= survey_area[surv];
                r_sum += surv_sum;
            }
        }

        dVeff += r*r*deltar*r_sum;
        zdVeff += zval*r*r*deltar*r_sum;

    }

    gsl_spline_free(Pmm_spline);
    gsl_spline_free(Pmt_spline);
    gsl_spline_free(Ptt_spline);
    gsl_interp_accel_free(Pmm_acc);
    gsl_interp_accel_free(Pmt_acc);
    gsl_interp_accel_free(Ptt_acc);

    return zdVeff/dVeff;

}

// The integrand for the integral over mu in the Fisher matrix calculation.
// For each mu we need to create a 4x4 matrix of the relevant power spectra derivatives and the inverse of the power spectrum matrix.
// Because there are some regions where the number density goes to zero we have to work directly with the inverse as it is difficult to invert numerically
// but if we deal with the inverse only then we can just set the relevant parts to zero when the number density is zero.
double mu_integrand(double mu, void * pin) {

    int i, j, m, q, u, surv;
    double * p = (double *)pin;
    double result, error;

    int numk = (int)p[0];
    double k = p[1];
    double zminval = p[4];
    double zmaxval = p[5];
    gsl_interp_accel * Pmm_acc, * Pmt_acc, * Ptt_acc;
    gsl_spline * Pmm_spline, * Pmt_spline, * Ptt_spline;
    double * Pmm_array = (double *)malloc(nzin*sizeof(double));
    double * Pmt_array = (double *)malloc(nzin*sizeof(double));
    double * Ptt_array = (double *)malloc(nzin*sizeof(double));
    for (j=0; j<nzin; j++) {
        Pmm_array[j] = pmmarray[j][numk];
        Pmt_array[j] = pmtarray[j][numk];
        Ptt_array[j] = pttarray[j][numk];
    }
    Pmm_acc    = gsl_interp_accel_alloc();
    Pmm_spline = gsl_spline_alloc(gsl_interp_cspline, nzin);
    gsl_spline_init(Pmm_spline, zin, Pmm_array, nzin);
    free(Pmm_array);

    Pmt_acc    = gsl_interp_accel_alloc();
    Pmt_spline = gsl_spline_alloc(gsl_interp_cspline, nzin);
    gsl_spline_init(Pmt_spline, zin, Pmt_array, nzin);
    free(Pmt_array);

    Ptt_acc    = gsl_interp_accel_alloc();
    Ptt_spline = gsl_spline_alloc(gsl_interp_cspline, nzin);
    gsl_spline_init(Ptt_spline, zin, Ptt_array, nzin);
    free(Ptt_array);

    double dendamp = sqrt(1.0/(1.0+0.5*(k*k*mu*mu*sigma_g*sigma_g)));     // This is unitless
    double veldamp = sin(k*sigma_u)/(k*sigma_u);                          // This is unitless

    double result_sum = 0.0;
    for (i=0; i<NRED[0]; i++) {

        double zval = zarray[i];
        double r_sum = 0.0;
        double r = rarray[i];
        double deltar = deltararray[i];

        if (zval < zminval) continue;
        if (zval > zmaxval) break;

        double dd_prefac=0.0, dv_prefac=0.0, vv_prefac=0.0;
        double P_gg=0.0, P_ug=0.0, P_uu=0.0;

        double sigma8 = sigma80 * growtharray[i];

        // First lets calculate the relevant power spectra. Interpolate the power spectra linearly in redshift
        double Pmm, Pmt, Ptt;
        Pmm = gsl_spline_eval(Pmm_spline, zval, Pmm_acc);
        Pmt = gsl_spline_eval(Pmt_spline, zval, Pmm_acc);
        Ptt = gsl_spline_eval(Ptt_spline, zval, Pmm_acc);

        double Omz = Om*ezinv(zval,NULL)*ezinv(zval,NULL)*(1.0+zval)*(1.0+zval)*(1.0+zval);
        double f = pow(Omz, gammaval);
        double beta = f*beta0*growtharray[i]/pow(Om,0.55);

        vv_prefac  = 1.0e2*f*mu*veldamp/k;
        dd_prefac = (1.0/(beta*beta) + 2.0*r_g*mu*mu/beta + mu*mu*mu*mu)*f*f*dendamp*dendamp;
        dv_prefac = (r_g/beta + mu*mu)*f*dendamp;
        P_gg = dd_prefac*Pmm;
        P_ug = vv_prefac*dv_prefac*Pmt;
        P_uu = vv_prefac*vv_prefac*Ptt;

        // And now the derivatives. Need to create a matrix of derivatives for each of the two parameters of interest
        gsl_matrix * dPdt1 = gsl_matrix_calloc(2, 2);
        gsl_matrix * dPdt2 = gsl_matrix_calloc(2, 2);
        double value;
        switch((int)p[2]) {
            // Differential w.r.t betaA
            case 0:
                value = -2.0*(1.0/beta + r_g*mu*mu)*f*f*dendamp*dendamp*Pmm/(beta*beta);
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = -(vv_prefac*f*r_g*dendamp*Pmt)/(beta*beta);
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                break;
            // Differential w.r.t fsigma8
            case 1:
                value = 2.0*(f/(beta*beta) + 2.0*f*r_g*mu*mu/beta + f*mu*mu*mu*mu)*dendamp*dendamp*Pmm/sigma8;
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = 2.0*vv_prefac*(r_g/beta + mu*mu)*dendamp*Pmt/sigma8;
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                value = (2.0*P_uu)/(f*sigma8);       
                gsl_matrix_set(dPdt1, 1, 1, value);
                break;
            // Differential w.r.t r_g
            case 2:
                value = 2.0*(1.0/beta)*mu*mu*f*f*dendamp*dendamp*Pmm;
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = vv_prefac*(1.0/beta)*f*dendamp*Pmt;
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                break;
            // Differential w.r.t sigma_g
            case 3:
                value = -k*k*mu*mu*dendamp*dendamp*sigma_g*P_gg;
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = -0.5*k*k*mu*mu*dendamp*dendamp*sigma_g*P_ug;
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                break;
            // Differential w.r.t sigma_u
            case 4:
                value = P_ug*(k*cos(k*sigma_u)/sin(k*sigma_u) - 1.0/sigma_u);
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                value = 2.0*P_uu*(k*cos(k*sigma_u)/sin(k*sigma_u) - 1.0/sigma_u);
                gsl_matrix_set(dPdt1, 1, 1, value);
                break;
            default:
                break;
        }
        switch((int)p[3]) {
            // Differential w.r.t betaA
            case 0:
                value = -2.0*(1.0/beta + r_g*mu*mu)*f*f*dendamp*dendamp*Pmm/(beta*beta);
                gsl_matrix_set(dPdt2, 0, 0, value);
                value = -(vv_prefac*f*r_g*dendamp*Pmt)/(beta*beta);
                gsl_matrix_set(dPdt2, 0, 1, value);
                gsl_matrix_set(dPdt2, 1, 0, value);
                break;
            // Differential w.r.t fsigma8
            case 1:
                value = 2.0*(f/(beta*beta) + 2.0*f*r_g*mu*mu/beta + f*mu*mu*mu*mu)*dendamp*dendamp*Pmm/sigma8;
                gsl_matrix_set(dPdt2, 0, 0, value);
                value = 2.0*vv_prefac*(r_g/beta + mu*mu)*dendamp*Pmt/sigma8;
                gsl_matrix_set(dPdt2, 0, 1, value);
                gsl_matrix_set(dPdt2, 1, 0, value);
                value = (2.0*P_uu)/(f*sigma8);       
                gsl_matrix_set(dPdt2, 1, 1, value);
                break;
            // Differential w.r.t r_g
            case 2:
                value = 2.0*(1.0/beta)*mu*mu*f*f*dendamp*dendamp*Pmm;
                gsl_matrix_set(dPdt2, 0, 0, value);
                value = vv_prefac*(1.0/beta)*f*dendamp*Pmt;
                gsl_matrix_set(dPdt2, 0, 1, value);
                gsl_matrix_set(dPdt2, 1, 0, value);
                break;
            // Differential w.r.t sigma_g
            case 3:
                value = -k*k*mu*mu*dendamp*dendamp*sigma_g*P_gg;
                gsl_matrix_set(dPdt2, 0, 0, value);
                value = -0.5*k*k*mu*mu*dendamp*dendamp*sigma_g*P_ug;
                gsl_matrix_set(dPdt2, 0, 1, value);
                gsl_matrix_set(dPdt2, 1, 0, value);
                break;
            // Differential w.r.t sigma_u
            case 4:
                value = P_ug*(k*cos(k*sigma_u)/sin(k*sigma_u) - 1.0/sigma_u);
                gsl_matrix_set(dPdt2, 0, 1, value);
                gsl_matrix_set(dPdt2, 1, 0, value);
                value = 2.0*P_uu*(k*cos(k*sigma_u)/sin(k*sigma_u) - 1.0/sigma_u);
                gsl_matrix_set(dPdt2, 1, 1, value);
                break;
            default:
                break;
        }

        // We need to do the overlapping and non-overlapping parts of the surveys separately
        for (surv=0; surv<3; surv++) {
            double surv_sum = 0.0;
            if (survey_area[surv] > 0.0) {
                double error_obs, error_noise, n_g = 0.0, n_u = 0.0;

                // Set the nbar for each section.
                if (surv == 0) {
                    n_g = nbararray[1][i];
                } else if (surv == 1) {
                    error_obs = 100.0*error_dist*r;                              // Percentage error * distance * H0 in km/s (factor of 100.0 comes from hubble parameter)
                    error_noise = error_rand*error_rand + error_obs*error_obs;   // Error_noise is in km^{2}s^{-2}
                    n_u = nbararray[0][i]/error_noise;                   
                } else {
                    error_obs = 100.0*error_dist*r;                              // Percentage error * distance * H0 in km/s (factor of 100.0 comes from hubble parameter)
                    error_noise = error_rand*error_rand + error_obs*error_obs;   // Error_noise is in km^{2}s^{-2}
                    n_u = nbararray[0][i]/error_noise;                   
                    n_g = nbararray[1][i];
                }

                //printf("%lf, %lf, %lf\n", r, n_g, 1.0e6*n_u);

                if (!((n_u > 0.0) || (n_g > 0.0))) continue;

                // First we need the determinant.
                double det = 1.0 + n_u*n_g*(P_gg*P_uu - P_ug*P_ug) + n_u*P_uu + n_g*P_gg;

                // Now the inverse matrix.
                gsl_matrix * iP = gsl_matrix_calloc(2, 2);
                value = n_u*n_g*P_uu + n_g;
                gsl_matrix_set(iP, 0, 0, value);
                value = n_g*n_u*P_gg + n_u;
                gsl_matrix_set(iP, 1, 1, value);
                value = - n_g*n_u*P_ug;
                gsl_matrix_set(iP, 0, 1, value);
                gsl_matrix_set(iP, 1, 0, value);
                
                // Finally we need to compute the Fisher integrand by summing over the inverse and differential matrices
                for (j=0; j<2; j++) {
                    for (m=0; m<2; m++) {
                        for (u=0; u<2; u++) {
                            for (q=0; q<2; q++) {
                                value = gsl_matrix_get(dPdt1, j, q)*gsl_matrix_get(iP, q, u)*gsl_matrix_get(dPdt2, u, m)*gsl_matrix_get(iP, m, j);
                                surv_sum += value;
                            }
                        }
                    }
                }
                surv_sum /= det*det;
                surv_sum *= survey_area[surv];
                r_sum += surv_sum;
                gsl_matrix_free(iP);
                //printf("%d, %lf, %lf, %lf, %lf\n", surv, k, mu, r, r_sum);

            }
        }
        //printf("%lf, %lf, %lf, %lf\n", k, mu, r, r_sum);

        result_sum += r*r*deltar*r_sum;

        gsl_matrix_free(dPdt1);
        gsl_matrix_free(dPdt2);
    }

    gsl_spline_free(Pmm_spline);
    gsl_spline_free(Pmt_spline);
    gsl_spline_free(Ptt_spline);
    gsl_interp_accel_free(Pmm_acc);
    gsl_interp_accel_free(Pmt_acc);
    gsl_interp_accel_free(Ptt_acc);

    return result_sum;
}


// Routine to read in the number density as a function of redshift. We need a file containing the left-most edge of each redshift bin and teh number density in that bin.
// From this we create arrays to store the bin centre, the bin width, the comoving distance and growth factor at the bin centre and the number density.
// The last bin width and bin centre is constructed from the last row of the input and the value of zmax at the top of the code. 
// ITS VERY IMPORTANT THAT THE NUMBER OF ROWS AND THE REDSHIFTS OF BOTH THE DENSITY AND PV NUMBER DENSITIES MATCH AS THE INTEGRATION OVER Z IS DONE USING THE TRAPEZIUM RULE.
void read_nz() {

    FILE * fp;
    char buf[500];
    int i, nsamp;

    NRED = (int *)calloc(2, sizeof(int));
    nbararray = (double **)calloc(2, sizeof(double*));
    double * zinarray;

    for (nsamp = 0; nsamp < 2; nsamp++) {

        if(!(fp = fopen(nbar_file[nsamp], "r"))) {
            printf("\nERROR: Can't open nbar file '%s'.\n\n", nbar_file[nsamp]);
            exit(0);
        }

        NRED[nsamp] = 0;
        while(fgets(buf,500,fp)) {
            if(strncmp(buf,"#",1)!=0) {
                double tz, tnbar;
                if(sscanf(buf, "%lf %lf\n", &tz, &tnbar) != 2) {printf("nbar read error\n"); exit(0);};
                if (tz > zmax) break;
                NRED[nsamp]++;
            }
        }
        fclose(fp);

        if (nsamp == 0) zinarray = (double *)calloc(NRED[nsamp], sizeof(double));
        nbararray[nsamp] = (double *)calloc(NRED[nsamp], sizeof(double));

        NRED[nsamp] = 0;
        fp = fopen(nbar_file[nsamp], "r");
        while(fgets(buf,500,fp)) {
            if(strncmp(buf,"#",1)!=0) {
                double tz, tnbar;
                if(sscanf(buf, "%lf %lf\n", &tz, &tnbar) != 2) {printf("nbar read error\n"); exit(0);};
                if (tz > zmax) break;
                if (nsamp == 0) zinarray[NRED[nsamp]] = tz;
                nbararray[nsamp][NRED[nsamp]] = 1.0e-6*tnbar;
                NRED[nsamp]++;
            }
        }
        fclose(fp);
    }

    if (NRED[1] != NRED[0]) {
        printf("ERROR: The number of redshift bins for each sample must match\n");
        exit(0);
    }   

    zarray = (double *)calloc(NRED[0], sizeof(double));
    rarray = (double *)calloc(NRED[0], sizeof(double));
    deltararray = (double *)calloc(NRED[0], sizeof(double));
    growtharray = (double *)calloc(NRED[0], sizeof(double));

    for (i=0; i<NRED[0]-1; i++) {
        zarray[i] = (zinarray[i+1]+zinarray[i])/2.0;
        rarray[i] = rz(zarray[i]);
        deltararray[i] = rz(zinarray[i+1]) - rz(zinarray[i]);
        growtharray[i] = growthz(zarray[i])/growthz(0.0);
        //printf("%12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n", zarray[i], rarray[i], deltararray[i], growtharray[i], nbararray[0][i], nbararray[1][i]);
    }
    zarray[NRED[0]-1] = (zmax+zinarray[NRED[0]-1])/2.0;
    rarray[NRED[0]-1] = rz(zarray[NRED[0]-1]);
    deltararray[NRED[0]-1] = rz(zmax) - rz(zinarray[NRED[0]-1]);
    growtharray[NRED[0]-1] = growthz(zarray[NRED[0]-1])/growthz(0.0);
    //printf("%12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n", zarray[NRED[0]-1], rarray[NRED[0]-1], deltararray[NRED[0]-1], growtharray[NRED[0]-1], nbararray[0][NRED[0]-1], nbararray[1][NRED[0]-1]);

    growth_acc    = gsl_interp_accel_alloc();
    growth_spline = gsl_spline_alloc(gsl_interp_cspline, NRED[0]);
    gsl_spline_init(growth_spline, zarray, growtharray, NRED[0]);

    free(zinarray);

    // Also create a simple redshift-distance spline
    int nbins = 400;
    double REDMIN = 0.0;
    double REDMAX = 2.0;
    double redbinwidth = (REDMAX-REDMIN)/(double)(nbins-1);
    double RMIN = rz(REDMIN);
    double RMAX = rz(REDMAX);
    double * ztemp = (double *)malloc(nbins*sizeof(double));
    double * rtemp = (double *)malloc(nbins*sizeof(double));
    for (i=0;i<nbins;i++) {
        ztemp[i] = i*redbinwidth+REDMIN;
        rtemp[i] = rz(ztemp[i]);
    }
    r_acc = gsl_interp_accel_alloc();
    r_spline = gsl_spline_alloc(gsl_interp_cspline, nbins);
    gsl_spline_init(r_spline, ztemp, rtemp, nbins);

    free(ztemp);
    free(rtemp);

    return;
}

// Routine to read in the velocity power spectrum.
void read_power() {
    
    FILE * fp;
    char buf[500];
    int i, j;

    pmmarray = (double**)malloc(nzin*sizeof(double*));
    pmtarray = (double**)malloc(nzin*sizeof(double*));
    pttarray = (double**)malloc(nzin*sizeof(double*));

    for (i = 0; i<nzin; i++) {

        char Pvel_file_in[500];
        sprintf(Pvel_file_in, "%s_z0p%02d.dat", Pvel_file, (int)(100.0*zin[i]));

        if(!(fp = fopen(Pvel_file_in, "r"))) {
            printf("\nERROR: Can't open power file '%s'.\n\n", Pvel_file_in);
            exit(0);
        }

        NK = 0;
        while(fgets(buf,500,fp)) {
            if(strncmp(buf,"#",1)!=0) {
                double tk, pkdelta, pkdeltavel, pkvel;
                if(sscanf(buf, "%lf %lf %lf %lf\n", &tk, &pkdelta, &pkdeltavel, &pkvel) != 4) {printf("Pvel read error\n"); exit(0);};
                NK++;
            }
        }
        fclose(fp);

        if (i == 0) {
            karray = (double *)calloc(NK, sizeof(double));
            deltakarray = (double *)calloc(NK-1, sizeof(double));
        }
        pmmarray[i] = (double *)calloc(NK, sizeof(double));
        pmtarray[i] = (double *)calloc(NK, sizeof(double));
        pttarray[i] = (double *)calloc(NK, sizeof(double));

        NK = 0;
        fp = fopen(Pvel_file_in, "r");
        while(fgets(buf,500,fp)) {
            if(strncmp(buf,"#",1)!=0) {
                double tk, pkdelta, pkdeltavel, pkvel;
                if(sscanf(buf, "%lf %lf %lf %lf\n", &tk, &pkdelta, &pkdeltavel, &pkvel) != 4) {printf("Pvel read error\n"); exit(0);};
                if (i == 0) karray[NK] = tk;
                pttarray[i][NK] = pkvel;
                pmmarray[i][NK] = pkdelta;
                pmtarray[i][NK] = pkdeltavel;
                NK++;
            }
        }
        fclose(fp);
    }

    for (i=0; i<NK-1; i++) deltakarray[i] = karray[i+1]-karray[i];

    pkkmin = karray[0];
    pkkmax = karray[NK-1];

    if (pkkmax < kmax) {
        printf("ERROR: The maximum k in the input power spectra id less than k_max\n");
        exit(0);
    }

    return;
}

// Integrand for the comoving distance
double ezinv(double x, void *p) {
  return 1.0/sqrt(Om*(1.0+x)*(1.0+x)*(1.0+x)+(1.0-Om));
}

// Calculates the comoving distance from the redshift
double rz(double red) {
  double result, error;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  F.function = &ezinv;
  gsl_integration_qags(&F, 0.0, red, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return c*result/100.0;
}

// The integrand for the normalised growth factor
double growthfunc(double x, void *p) {
    double red = 1.0/x - 1.0;
    double Omz = Om*ezinv(red,NULL)*ezinv(red,NULL)/(x*x*x);
    double f = pow(Omz, gammaval);
    return f/x;
}

// Calculates the normalised growth factor as a function of redshift given a value of gammaval
double growthz(double red) {
  double result, error;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  F.function = &growthfunc;
  double a = 1.0/(1.0+red);
  gsl_integration_qags(&F, a, 1.0, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return exp(-result);
}



