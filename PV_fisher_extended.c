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

// Fisher matrix calculation for surveys with velocity and density field measurements. This extended version includes the extra parameters
// associated with primordial non-Gaussianity, scale-dependent density and velocity bias, and zero-point offsets, as per Howlett 2017a.
// ASSUMPTIONS:
//  -  Uncorrelated shot-noise between the density and velocity fields
//  -  I use the trapezium rule to integrate over r. There will be some error due to this, but this makes the most sense as 
//     we are binning the number density anyway, and it makes hella difference in the speed of the code.
//  -  The redshift dependence of the non-linear matter and velocity divergence power spectra is captured using linear interpolation.
//  -  The PV error scales as a fixed pecentage of H0*r.
//  -  Flat LCDM cosmology (but not necessarily GR as gammaval can be changed).
//  -  The damping of the velocity and density fields due to non-linear RSD is redshift independent
//  -  The scale dependence of the scale-dependent bias is redshift independent

// The parameters necessary for the calculation
static int nparams = 2;           // The number of free parameters (we can use any of beta, fsigma8, r_g, sigma_g, sigma_u, fnl, Rv, Rd, zp_err)
static int Data[2] = {0,1};     // A vector of flags for the parameters we are interested in (0=beta, 1=fsigma8, 2=r_g, 3=sigma_g, 4=sigma_u, 5=fnl, 6=Rv, 7=Rd, 8=zp_err). MAKE SURE THE LENGTH OF THIS VECTOR, NPARAMS AND THE ENTRIES AGREE/MAKE SENSE, OR YOU MIGHT GET NONSENSE RESULTS!!
static int nziter = 1;            // Now many bins in redshift between zmin and zmax we are considering
static double zmin = 0.0;         // The minimum redshift to consider (You must have power spectra that are within this range or GSL spline will error out)
static double zmax = 0.4;         // The maximum redshift to consider (You must have power spectra that are within this range or GSL spline will error out)
static double Om = 0.3121;        // The matter density at z=0
static double c = 299792.458;     // The speed of light in km/s
static double gammaval = 0.55;    // The value of gammaval to use in the forecasts (where f(z) = Om(z)^gammaval)
static double r_g = 1.0;          // The cross correlation coefficient between the velocity and density fields
static double beta0 = 0.493;      // The value of beta (at z=0, we'll modify this by the redshift dependent value of bias and f as required)
static double sigma80 = 0.8150;   // The value of sigma8 at z=0
static double sigma_u = 13.00;    // The value of the velocity damping parameter in Mpc/h. I use the values from Jun Koda's paper
static double sigma_g = 4.24;     // The value of the density damping parameter in Mpc/h. I use the values from Jun Koda's paper
static double fnl = 0.0;          // The value of fnl contributing to primordial non-Gaussianity
static double Rd = 0.0;           // Scale-dependent bias following the parameterisation in Howlett2017a
static double Rv = 0.0;           // Velocity bias following the parameterisation in Howlett2017a
static double zp_err = 0.05;      // The error in the zeropoint offset if we are interested in the systematic bias caused by such an error (sigma_m in Eq. 34 of Howlett 2017a)
static double do_bias = 0;        // Whether or not to calculate the systematic parameter bias due to neglecting scale-dependent biases.
static double do_zp_bias = 1;     // Whether or not to calculate the systematic parameter bias due to having a zero-point offset.
static double kmax = 0.2;         // The maximum k to evaluate for dd, dv and vv correlations (Typical values are 0.1 - 0.2, on smaller scales the models are likely to break down).
static double survey_area[3] = {0.0, 0.0, 1.65};   // We need to know the survey area for each survey and the overlap area between the surveys (redshift survey only first, then PV survey only, then overlap. 
                                                    // For fully overlapping we would have {0, 0, size_overlap}. For redshift larger than PV, we would have {size_red-size_overlap, 0, size_overlap}). Units are pi steradians, such that full sky is 4.0, half sky is 2.0 etc.
static double error_rand = 300.0;    // The observational error due to random non-linear velocities (I normally use 300km/s as in Jun Koda's paper)
static double error_dist = 0.20;     // The percentage error on the distance indicator (Typically 0.05 - 0.10 for SNe IA, 0.2 or more for Tully-Fisher or Fundamental Plane) 
static double verbosity = 2;         // How much output to give: 0 = only percentage errors on fsigma8, 1 = other useful info and nuisance parameters, 2 = full fisher and covariance matrices

// The number of redshifts and the redshifts themselves of the input matter and velocity divergence power spectra. 
// These numbers are multiplied by 100, converted to ints and written in the form _z0p%02d which is then appended to the filename Pvel_file. See routine read_power. 
static double nzin = 11;
static double zin[11] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
char * Pvel_file = "./example_files/example_pk";          // The file containing the velocity divergence power spectrum. Don't include .dat as we'll append the redshifts on read in
char * Trans_file = "./example_files/example_transfer.dat";   // The file containing the transfer function. Only used for primordial non-Gaussianity and must be normalised to 1 as k->0


// The files containing the number density of the surveys. First is the PV survey, then the redshift survey. These files MUST have the same binning and redshift range, 
// so that the sum over redshift bins works (would be fine if we used splines), i.e., if one survey is shallower then that file must contain rows with n(z)=0.
// I also typically save nbar x 10^6 in the input file to make sure I don't lose precision when outputting small nbar values to files. This is corrected when the nbar file
// is read in, so see the read_nz() routine!
char * nbar_file[300] = {"./example_files/taipan_vel_nbar_final_reformat.dat",
                         "./example_files/taipan_red_nbar_final_reformat.dat"};      

// Other global parameters and arrays
int NK, * NRED;
double pkkmin;         // The minimum kmin to integrate over, based on the input power spectrum file
double pkkmax;         // The maximum k in the input power spectrum. The maximum k to integrate over is the smallest of this or kmax
double * zarray;
double * rarray;
double * ezarray;
double * deltararray;
double * growtharray;
double ** nbararray;
double * karray, * deltakarray;
double ** pmmarray, ** pmtarray, ** pttarray;
gsl_spline * growth_spline, * r_spline;
gsl_interp_accel * growth_acc, * r_acc;
gsl_spline * Trans_spline, * growth_spline, * r_spline;
gsl_interp_accel * Trans_acc, * growth_acc, * r_acc;

// Prototypes
double zeff_integrand(double mu, void * pin);
double mu_integrand(double mu, void * pin);
double ezinv(double x, void *p);
double rz(double red);
double growthfunc(double x, void *p);
double growthz(double red);
void read_nz();
void read_power();
void read_transfer();
double calc_zeff(double kmin, double zmin_iter, double zmax_iter);
void calc_fisher(double kmin, double zmin_iter, double zmax_iter, double z_eff, double bias_flag, gsl_matrix * Covariance, double * Sys);

// Calculates the fished matrix for a velocity survey.
int main(int argc, char **argv) {
    
    int i, j;

    // Read in the velocity divergence power spectrum output from the COPTER code (Carlson 2009)
    read_power();

    // Read in the transfer function
    read_transfer();

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

    printf("Evaluating the Fisher Matrix for %d bins between [z_min = %lf, z_max = %lf]\n", nziter, zmin, zmax);

    if (verbosity == 0) {
        if (do_bias || do_zp_bias) {
            printf("#     zmin         zmax         zeff      fsigma8(z_eff)   TRUE percentage error(z_eff)     percentage bias(z_eff)       INCORRECT percentage error(z_eff)\n");
        } else {
            printf("#     zmin         zmax         zeff      fsigma8(z_eff)   percentage error(z_eff)\n");
        }
    }

    int ziter;
    double * Sys = (double *)malloc(nparams*sizeof(double*));
    gsl_matrix * Covariance = gsl_matrix_alloc(nparams, nparams);
    for (ziter = 0; ziter<nziter; ziter++) {

        double zbinwidth = (zmax-zmin)/(nziter);
        double zmin_iter = ziter*zbinwidth + zmin;
        double zmax_iter = (ziter+1.0)*zbinwidth + zmin;

        double rzmax = gsl_spline_eval(r_spline, zmax_iter, r_acc);
        double kmin = M_PI/rzmax;

        if (verbosity > 0) printf("Evaluating the Fisher Matrix for [k_min = %lf, k_max = %lf] and [z_min = %lf, z_max = %lf]\n", kmin, kmax, zmin_iter, zmax_iter);

        // Calculate the effective redshift (which I base on the sum of the S/N for the density and velocity fields)
        // *********************************************************************************************************
        double z_eff = calc_zeff(kmin, zmin_iter, zmax_iter);

        // Fisher matrix using the TRUE covariance matrix (i.e., with scale-dependent biases). This routine changes based on whether do_zp_bias is true or false.
        // In the former case, this Fisher Matrix DOESN'T include a zero-point offset, as we compute the systematically biased case later. If do_zp_bias is false
        // we are not interested in the bias caused by an offset, and so we treat whatever value is in zp_err as real and the Fisher Matrix DOES include a potential zero-point offset
        // ******************************************************************************************************************
        calc_fisher(kmin, zmin_iter, zmax_iter, z_eff, 0, Covariance, Sys);

        // Bias vector (i.e., Eq. 50 from Howlett 2017a, which uses computes the parameter offset due to neglecting scale-dependent biases or adding a zero-point offset)
        // ******************************************************************************************************************
        if (do_bias || do_zp_bias) calc_fisher(kmin, zmin_iter, zmax_iter, z_eff, 1, Covariance, Sys);

        // Fisher matrix using the INCORRECT covariance matrix (i.e., without scale-dependent biases, and with zero-point offsets). This is used for comparing the 
        // magnitude of any bias with respect to the errors (as the errors could also be biased too small!)
        // ******************************************************************************************************************
        if (do_bias || do_zp_bias) calc_fisher(kmin, zmin_iter, zmax_iter, z_eff, 2, Covariance, Sys);

    }

    // Now the full Fisher matrix over all redshifts if we had more than 1 redshift bin. This is equivalent to assuming a single measurement from all data
    // and NOT the same as summing the individual matrices above (which is equivalent to making separate measurements at each redshift, then marginalising over
    // all redshift-dependent quantities, and so may be different for things like fnl).
    if (nziter > 1) {
        double rzmax = gsl_spline_eval(r_spline, zmax, r_acc);
        double kmin = M_PI/rzmax;

        if (verbosity > 0) printf("Finally, evaluating the Fisher Matrix for [k_min = %lf, k_max = %lf] and [z_min = %lf, z_max = %lf]\n", kmin, kmax, zmin, zmax);

        // Calculate the effective redshift (which I base on the sum of the S/N for the density and velocity fields)
        // *********************************************************************************************************
        double z_eff = calc_zeff(kmin, zmin, zmax);

        // Fisher matrix using the TRUE covariance matrix (i.e., with scale-dependent biases). This routine changes based on whether do_zp_bias is true or false.
        // In the former case, this Fisher Matrix DOESN'T include a zero-point offset, as we compute the systematically biased case later. If do_zp_bias is false
        // we are not interested in the bias caused by an offset, and so we treat whatever value is in zp_err as real and the Fisher Matrix DOES include a potential zero-point offset
        // ******************************************************************************************************************
        calc_fisher(kmin, zmin, zmax, z_eff, 0, Covariance, Sys);

        // Bias vector (i.e., Eq. 50 from Howlett 2017a, which uses computes the parameter offset due to neglecting scale-dependent biases or adding a zero-point offset)
        // ******************************************************************************************************************
        if (do_bias || do_zp_bias) calc_fisher(kmin, zmin, zmax, z_eff, 1, Covariance, Sys);

        // Fisher matrix using the INCORRECT covariance matrix (i.e., without scale-dependent biases, and with zero-point offsets). This is used for comparing the 
        // magnitude of any bias with respect to the errors (as the errors could also be biased too small!)
        // ******************************************************************************************************************
        if (do_bias || do_zp_bias) calc_fisher(kmin, zmin, zmax, z_eff, 2, Covariance, Sys);

    }
        
    free(Sys);
    gsl_matrix_free(Covariance);
    gsl_spline_free(growth_spline);
    gsl_interp_accel_free(growth_acc);
    gsl_spline_free(Trans_spline);
    gsl_interp_accel_free(Trans_acc);

    return 0;
}

double calc_zeff(double kmin, double zmin_iter, double zmax_iter) {

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

    return z_eff;

}

// Calculate the fisher matrix, integrating over k, then mu, then r (r is last as it means we are effectively integrating over effective volume).
// As the input spectra are tabulated we'll just use the trapezium rule to integrate over k
void calc_fisher(double kmin, double zmin_iter, double zmax_iter, double z_eff, double bias_flag, gsl_matrix * Covariance, double * Sys) {

    int i, j, numk;

    double growth_eff = gsl_spline_eval(growth_spline, z_eff, growth_acc);

    double sigma8 = sigma80 * growth_eff;
    double Omz = Om*ezinv(z_eff,NULL)*ezinv(z_eff,NULL)*(1.0+z_eff)*(1.0+z_eff)*(1.0+z_eff);
    double f = pow(Omz, gammaval);
    double beta = f*beta0*growth_eff/pow(Om,0.55);

    if (verbosity > 0) {
        if (do_zp_bias) {
            if (bias_flag == 0) {
                printf("TRUE Fisher matrix without zero-point offset:\n");
            } else if (bias_flag == 1) {
                printf("BIAS vector:\n");
            } else if (bias_flag == 2) {
                printf("INCORRECT Fisher matrix:\n");
            }
        } else {
            if (bias_flag == 0) {
                printf("TRUE Fisher matrix with potential zero-point offset set by zp_err:\n");
            } else if (bias_flag == 1) {
                printf("BIAS vector:\n");
            } else if (bias_flag == 2) {
                printf("INCORRECT Fisher matrix:\n");
            }
        }
    }

    double * Bias = (double *)malloc(nparams*sizeof(double*));
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
                double params[7] = {numk, k, Data[i], Data[j], zmin_iter, zmax_iter, bias_flag};

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
            if (bias_flag == 1) {
                Bias[i] = k_sum/(4.0*M_PI);
            } else {
                gsl_matrix_set(Fisher, i, j, k_sum/(4.0*M_PI));
                gsl_matrix_set(Fisher, j, i, k_sum/(4.0*M_PI));
            }
        }
    }

    if (bias_flag == 1) {

        for (i=0; i<nparams; i++) {
            Sys[i] = 0.0;
            for (j=0; j<nparams; j++) {
                Sys[i] += gsl_matrix_get(Covariance, i, j)*Bias[j];
            }
        }

        if (verbosity > 0) {
            for (i=0; i<nparams; i++) {
                if (Data[i] == 0) {
                    printf("%12.6lf bias on beta = %12.6lf\n", Sys[i], beta);
                }
                if (Data[i] == 1) {
                    printf("%12.6lf bias on fsigma8 = %12.6lf\n", Sys[i], f*sigma8);
                }
                if (Data[i] == 2) {
                    printf("%12.6lf bias on r_g = %12.6lf\n", Sys[i], r_g);
                }
                if (Data[i] == 3) {
                    printf("%12.6lf bias on sigma_g = %12.6lf\n", Sys[i], sigma_g);
                }
                if (Data[i] == 4) {
                    printf("%12.6lf bias on sigma_u = %12.6lf\n", Sys[i], sigma_u);
                }
                if (Data[i] == 5) {
                    printf("%12.6lf bias on fnl = %12.6lf\n", Sys[i], fnl);
                }
            }
        }

    } else {

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
        p = gsl_permutation_alloc(nparams);
        gsl_linalg_LU_decomp(Fisher, p, &s);
        gsl_linalg_LU_invert(Fisher, p, Covariance);
        gsl_permutation_free(p);

        if (verbosity == 0) {
            for (i=0; i<nparams; i++) {
                if (Data[i] == 1) {
                    if (do_bias || do_zp_bias) {
                        if (bias_flag == 0) {
                            printf("%12.6lf  %12.6lf  %12.6lf  %12.6lf        %12.6lf        ", zmin_iter, zmax_iter, z_eff, f*sigma8, 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/(f*sigma8));
                        } else if (bias_flag == 2) {
                            printf("          %12.6lf          %12.6lf\n", 100.0*Sys[i]/sqrt(gsl_matrix_get(Covariance, i, i)), 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/(f*sigma8));
                        }
                    } else {
                        printf("%12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf\n", zmin_iter, zmax_iter, z_eff, f*sigma8, 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/(f*sigma8));
                    }
                }
            }
        }

        if (verbosity > 0) {
            printf("Parameter constraints\n======================================================\n");
            for (i=0; i<nparams; i++) {
                if (Data[i] == 0) {
                    printf("beta = %12.6lf +/- %12.6lf\n", beta, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on beta\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/beta);
                    if (bias_flag == 2) printf("%4.2lf percent bias on beta\n", 100.0*Sys[i]/sqrt(gsl_matrix_get(Covariance, i, i)));
                }
                if (Data[i] == 1) {
                    printf("fsigma8 = %12.6lf +/- %12.6lf\n", f*sigma8, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on fsigma8\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/(f*sigma8));
                    if (bias_flag == 2) printf("%4.2lf percent bias on fsigma8\n", 100.0*Sys[i]/sqrt(gsl_matrix_get(Covariance, i, i)));
                }
                if (Data[i] == 2) {
                    printf("r_g = %12.6lf +/- %12.6lf\n", r_g, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on r_g\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/r_g);
                    if (bias_flag == 2) printf("%4.2lf percent bias on r_g\n", 100.0*Sys[i]/sqrt(gsl_matrix_get(Covariance, i, i)));
                }
                if (Data[i] == 3) {
                    printf("sigma_g = %12.6lf +/- %12.6lf\n", sigma_g, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on sigma_g\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/sigma_g);
                    if (bias_flag == 2) printf("%4.2lf percent bias on sigma_g\n", 100.0*Sys[i]/sqrt(gsl_matrix_get(Covariance, i, i)));
                }
                if (Data[i] == 4) {
                    printf("sigma_u = %12.6lf +/- %12.6lf\n", sigma_u, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on sigma_u\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/sigma_u);
                    if (bias_flag == 2) printf("%4.2lf percent bias on sigma_u\n", 100.0*Sys[i]/sqrt(gsl_matrix_get(Covariance, i, i)));
                }
                if (Data[i] == 5) {
                    printf("fnl = %12.6lf +/- %12.6lf\n", fnl, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf absolute error on fnl\n", sqrt(gsl_matrix_get(Covariance, i, i)));
                    if (bias_flag == 2) printf("%4.2lf percent bias on fnl\n", 100.0*Sys[i]/sqrt(gsl_matrix_get(Covariance, i, i)));
                }
                if (Data[i] == 6) {
                    printf("Rv = %12.6lf +/- %12.6lf\n", Rv, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on Rv\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/Rv);
                }
                if (Data[i] == 7) {
                    printf("Rd = %12.6lf +/- %12.6lf\n", Rd, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on Rd\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/Rd);
                }
                if (Data[i] == 8) {
                    printf("zp_err = %12.6lf +/- %12.6lf\n", zp_err, sqrt(gsl_matrix_get(Covariance, i, i)));
                    printf("%4.2lf percent error on zp_err\n", 100.0*sqrt(gsl_matrix_get(Covariance, i, i))/zp_err);
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

    free(Bias);
    gsl_matrix_free(Fisher);

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

    double Tk = gsl_spline_eval(Trans_spline, k, Trans_acc);
    double Ak0 = 3.0*1.686*Om*1.0e4/(k*k*Tk*c*c);                         // The prefactor for the scale dependent bias induced by primordial non-Gaussianity
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

        double Ak = Ak0/growtharray[i];
        double sigma8 = sigma80 * growtharray[i];

        // First lets calculate the relevant power spectra. Interpolate the power spectra linearly in redshift
        double Pmm, Pmt, Ptt;
        Pmm = gsl_spline_eval(Pmm_spline, zval, Pmm_acc);
        Pmt = gsl_spline_eval(Pmt_spline, zval, Pmm_acc);
        Ptt = gsl_spline_eval(Ptt_spline, zval, Pmm_acc);

        double Omz = Om*ezinv(zval,NULL)*ezinv(zval,NULL)*(1.0+zval)*(1.0+zval)*(1.0+zval);
        double f = pow(Omz, gammaval);
        double beta = f*beta0*growtharray[i]/pow(Om,0.55);
        double bv = 1.0 - Rv*Rv*k*k;
        double bsd = Rd*k*k;
        double beta_sd = 1.0/beta + bsd/f;

        vv_prefac  = 1.0e2*bv*f*mu*veldamp/k;
        dd_prefac = ((1.0+fnl*Ak)*(1.0+fnl*Ak)*(beta_sd*beta_sd) + 2.0*bv*r_g*mu*mu*(1.0+fnl*Ak)*beta_sd + bv*bv*mu*mu*mu*mu - 2.0*fnl*Ak*(1.0+fnl*Ak)*beta_sd/f - 2.0*bv*r_g*mu*mu*fnl*Ak/f + fnl*fnl*Ak*Ak/(f*f))*f*f*dendamp*dendamp;
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
// For each mu we need to create a 2x2 matrix of the relevant power spectra derivatives and the inverse of the power spectrum matrix.
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

    double Tk = gsl_spline_eval(Trans_spline, k, Trans_acc);
    double Ak0 = 3.0*1.686*Om*1.0e4/(k*k*Tk*c*c);                         // The prefactor for the scale dependent bias induced by primordial non-Gaussianity
    double dendamp = sqrt(1.0/(1.0+0.5*(k*k*mu*mu*sigma_g*sigma_g)));     // This is unitless
    double veldamp = sin(k*sigma_u)/(k*sigma_u);                          // This is unitless

    double result_sum = 0.0;
    for (i=0; i<NRED[0]; i++) {

        double zval = zarray[i];
        double r_sum = 0.0;
        double r = rarray[i];
        double ez = ezarray[i];
        double deltar = deltararray[i];

        if (zval < zminval) continue;
        if (zval > zmaxval) break;

        double dd_prefac=0.0, dv_prefac=0.0, vv_prefac=0.0;
        double dd_prefac_bias=0.0, dv_prefac_bias=0.0, vv_prefac_bias=0.0;
        double P_gg=0.0, P_ug=0.0, P_uu=0.0;
        double P_gg_bias=0.0, P_ug_bias=0.0, P_uu_bias=0.0;

        double Ak = Ak0/growtharray[i];
        double sigma8 = sigma80 * growtharray[i];

        // First lets calculate the relevant power spectra. Interpolate the power spectra linearly in redshift
        double Pmm, Pmt, Ptt;
        Pmm = gsl_spline_eval(Pmm_spline, zval, Pmm_acc);
        Pmt = gsl_spline_eval(Pmt_spline, zval, Pmm_acc);
        Ptt = gsl_spline_eval(Ptt_spline, zval, Pmm_acc);

        double Omz = Om*ezinv(zval,NULL)*ezinv(zval,NULL)*(1.0+zval)*(1.0+zval)*(1.0+zval);
        double f = pow(Omz, gammaval);
        double beta = f*beta0*growtharray[i]/pow(Om,0.55);
        double bv = 1.0 - Rv*Rv*k*k;
        double bsd = Rd*k*k;
        double beta_sd = 1.0/beta + bsd/f;

        // If p[6] (bias_flag) == 2 this means we are computing the Fisher matrix INCLUDING systematics, and so we need the covariance matrix without scale-dependent bias
        if ((int)p[6] == 2) {
            bv = 1.0;
            beta_sd = 1.0/beta;
        }

        // In the presence of a zero-point offset we have to calculate the additional velocity recieved. 
        // This is redshift dependent if we assume a fixed offset in the logarithmic distance ratio
        double zp_err_prefac = c*log(10.0)/(1.0 - ((c*(1.0 + zval)*(1.0 + zval))/(100.0*ez*r)))/5.0;
        double zp_err_v = zp_err*zp_err_prefac;

        vv_prefac  = 1.0e2*bv*f*mu*veldamp/k;
        dd_prefac = ((1.0+fnl*Ak)*(1.0+fnl*Ak)*(beta_sd*beta_sd) + 2.0*bv*r_g*mu*mu*(1.0+fnl*Ak)*beta_sd + bv*bv*mu*mu*mu*mu - 2.0*fnl*Ak*(1.0+fnl*Ak)*beta_sd/f - 2.0*bv*r_g*mu*mu*fnl*Ak/f + fnl*fnl*Ak*Ak/(f*f))*f*f*dendamp*dendamp;
        dv_prefac = (r_g*(1.0+fnl*Ak)*beta_sd - r_g*fnl*Ak/f + bv*mu*mu)*f*dendamp;

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
                value = -2.0*((1.0+fnl*Ak)*(1.0+fnl*Ak)*beta_sd + bv*r_g*mu*mu*(1.0+fnl*Ak) - fnl*Ak*(1.0+fnl*Ak)/f)*f*f*dendamp*dendamp*Pmm/(beta*beta);
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = -(vv_prefac*f*r_g*(1.0+fnl*Ak)*dendamp*Pmt)/(beta*beta);
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                break;
            // Differential w.r.t fsigma8
            case 1:
                value = 2.0*(f*(1.0+fnl*Ak)*(1.0+fnl*Ak)*beta_sd/beta + f*r_g*bv*mu*mu*(1.0+fnl*Ak)*(2.0/beta + bsd/f) + f*bv*bv*mu*mu*mu*mu - fnl*Ak*(1.0+fnl*Ak)/beta - bv*r_g*mu*mu*fnl*Ak)*dendamp*dendamp*Pmm/sigma8;
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = vv_prefac*(r_g*(1.0+fnl*Ak)*(2.0/beta + bsd/f) - r_g*fnl*Ak/f + bv*mu*mu)*dendamp*Pmt/sigma8;
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                value = (2.0*P_uu)/(f*sigma8);       
                gsl_matrix_set(dPdt1, 1, 1, value);
                break;
            // Differential w.r.t r_g
            case 2:
                value = 2.0*bv*mu*mu*((1.0+fnl*Ak)*beta_sd - fnl*Ak/f)*f*f*dendamp*dendamp*Pmm;
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = vv_prefac*((1.0+fnl*Ak)*beta_sd - fnl*Ak/f)*f*dendamp*Pmt;
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
            // Differential w.r.t fnl
            case 5:
                value = 2.0*Ak*((1.0+fnl*Ak)*beta_sd*beta_sd + bv*r_g*mu*mu*(beta_sd - 1.0/f) - (1.0+2.0*fnl*Ak)*beta_sd/f + fnl*Ak/(f*f))*f*f*dendamp*dendamp*Pmm;
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = vv_prefac*r_g*Ak*(beta_sd - 1.0/f)*f*dendamp*Pmt;
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                break;    
            // Differential w.r.t Rv
            case 6:
                value = -4.0*Rv*k*k*(r_g*mu*mu*((1.0+fnl*Ak)*beta_sd - fnl*Ak/f) + bv*mu*mu*mu*mu)*f*f*dendamp*dendamp*Pmm;
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = -2.0*Rv*k*k*vv_prefac*(r_g*(1.0+fnl*Ak)*beta_sd/bv - r_g*fnl*Ak/(f*bv) + 2.0*mu*mu)*f*dendamp*Pmt;
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                value = -4.0*Rv*k*k*P_uu/bv;      
                gsl_matrix_set(dPdt1, 1, 1, value);
                break;
            // Differential w.r.t Rd
            case 7:
                value = 2.0*k*k*(f*(1.0+fnl*Ak)*(1.0+fnl*Ak)*beta_sd + bv*r_g*f*mu*mu*(1.0+fnl*Ak) - fnl*Ak*(1.0+fnl*Ak))*dendamp*dendamp*Pmm;
                gsl_matrix_set(dPdt1, 0, 0, value);
                value = vv_prefac*k*k*r_g*(1.0+fnl*Ak)*dendamp*Pmt;
                gsl_matrix_set(dPdt1, 0, 1, value);
                gsl_matrix_set(dPdt1, 1, 0, value);
                break;
            // Differential w.r.t. zp_err
            case 8:
                if (nbararray[0][i] > 0.0) {
                    value = 2.0*zp_err_v*zp_err_prefac/nbararray[0][i];
                    gsl_matrix_set(dPdt1, 1, 1, value);
                }
                break;
            default:
                break;
        }
        // If we want to calculate the Bias vector and look at systematics from scale-dependent bias or zero-point offsets, we simply have to replace dPdt2 with the difference between the models with and without bias. 
        // The subscript "bias" is the systematically biased covariance matrix, which, in the most confusing way possible, is the one which doesn't have the velocity/scale-dependent bias included or does have the zero-point offset
        if ((int)p[6] == 1) {
            // Calculate the model without scale-dependent biases
            if (do_bias) {
                double bv_bias = 1.0;
                double beta_sd_bias = 1.0/beta;
                vv_prefac_bias  = 1.0e2*bv_bias*f*mu*veldamp/k;
                dd_prefac_bias = ((1.0+fnl*Ak)*(1.0+fnl*Ak)*(beta_sd_bias*beta_sd_bias) + 2.0*bv_bias*r_g*mu*mu*(1.0+fnl*Ak)*beta_sd_bias + bv_bias*bv_bias*mu*mu*mu*mu - 2.0*fnl*Ak*(1.0+fnl*Ak)*beta_sd_bias/f - 2.0*bv_bias*r_g*mu*mu*fnl*Ak/f + fnl*fnl*Ak*Ak/(f*f))*f*f*dendamp*dendamp;
                dv_prefac_bias = (r_g*(1.0+fnl*Ak)*beta_sd_bias - r_g*fnl*Ak/f + bv_bias*mu*mu)*f*dendamp;
                P_gg_bias = dd_prefac_bias*Pmm;
                P_ug_bias = vv_prefac_bias*dv_prefac_bias*Pmt;
                P_uu_bias = vv_prefac_bias*vv_prefac_bias*Ptt;
                gsl_matrix_set(dPdt2, 0, 0, P_gg - P_gg_bias);
                gsl_matrix_set(dPdt2, 0, 1, P_ug - P_ug_bias);
                gsl_matrix_set(dPdt2, 1, 0, P_ug - P_ug_bias);
                gsl_matrix_set(dPdt2, 1, 1, P_uu - P_uu_bias);
                //printf("%lf, %lf, %lf, %lf\n",bvA, P_uuAA, P_uuAA_bias, gsl_matrix_get(dPdt2, 1, 1));
            }
            if (do_zp_bias) {
                if (nbararray[0][i] > 0.0) {
                    gsl_matrix_set(dPdt2, 1, 1, zp_err_v*zp_err_v/nbararray[0][i]);
                }
            }
        } else {
            switch((int)p[3]) {
                // Differential w.r.t betaA
                case 0:
                    value = -2.0*((1.0+fnl*Ak)*(1.0+fnl*Ak)*beta_sd + bv*r_g*mu*mu*(1.0+fnl*Ak) - fnl*Ak*(1.0+fnl*Ak)/f)*f*f*dendamp*dendamp*Pmm/(beta*beta);
                    gsl_matrix_set(dPdt2, 0, 0, value);
                    value = -(vv_prefac*f*r_g*(1.0+fnl*Ak)*dendamp*Pmt)/(beta*beta);
                    gsl_matrix_set(dPdt2, 0, 1, value);
                    gsl_matrix_set(dPdt2, 1, 0, value);
                    break;
                // Differential w.r.t fsigma8
                case 1:
                    value = 2.0*(f*(1.0+fnl*Ak)*(1.0+fnl*Ak)*beta_sd/beta + f*r_g*bv*mu*mu*(1.0+fnl*Ak)*(2.0/beta + bsd/f) + f*bv*bv*mu*mu*mu*mu - fnl*Ak*(1.0+fnl*Ak)/beta - bv*r_g*mu*mu*fnl*Ak)*dendamp*dendamp*Pmm/sigma8;
                    gsl_matrix_set(dPdt2, 0, 0, value);
                    value = vv_prefac*(r_g*(1.0+fnl*Ak)*(2.0/beta + bsd/f) - r_g*fnl*Ak/f + bv*mu*mu)*dendamp*Pmt/sigma8;
                    gsl_matrix_set(dPdt2, 0, 1, value);
                    gsl_matrix_set(dPdt2, 1, 0, value);
                    value = (2.0*P_uu)/(f*sigma8);       
                    gsl_matrix_set(dPdt2, 1, 1, value);
                    break;
                // Differential w.r.t r_g
                case 2:
                    value = 2.0*bv*mu*mu*((1.0+fnl*Ak)*beta_sd - fnl*Ak/f)*f*f*dendamp*dendamp*Pmm;
                    gsl_matrix_set(dPdt2, 0, 0, value);
                    value = vv_prefac*((1.0+fnl*Ak)*beta_sd - fnl*Ak/f)*f*dendamp*Pmt;
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
                // Differential w.r.t fnl
                case 5:
                    value = 2.0*Ak*((1.0+fnl*Ak)*beta_sd*beta_sd + bv*r_g*mu*mu*(beta_sd - 1.0/f) - (1.0+2.0*fnl*Ak)*beta_sd/f + fnl*Ak/(f*f))*f*f*dendamp*dendamp*Pmm;
                    gsl_matrix_set(dPdt2, 0, 0, value);
                    value = vv_prefac*r_g*Ak*(beta_sd - 1.0/f)*f*dendamp*Pmt;
                    gsl_matrix_set(dPdt2, 0, 1, value);
                    gsl_matrix_set(dPdt2, 1, 0, value);
                    break;    
                // Differential w.r.t Rv
                case 6:
                    value = -4.0*Rv*k*k*(r_g*mu*mu*((1.0+fnl*Ak)*beta_sd - fnl*Ak/f) + bv*mu*mu*mu*mu)*f*f*dendamp*dendamp*Pmm;
                    gsl_matrix_set(dPdt2, 0, 0, value);
                    value = -2.0*Rv*k*k*vv_prefac*(r_g*(1.0+fnl*Ak)*beta_sd/bv - r_g*fnl*Ak/(f*bv) + 2.0*mu*mu)*f*dendamp*Pmt;
                    gsl_matrix_set(dPdt2, 0, 1, value);
                    gsl_matrix_set(dPdt2, 1, 0, value);
                    value = -4.0*Rv*k*k*P_uu/bv;      
                    gsl_matrix_set(dPdt2, 1, 1, value);
                    break;
                // Differential w.r.t Rd
                case 7:
                    value = 2.0*k*k*(f*(1.0+fnl*Ak)*(1.0+fnl*Ak)*beta_sd + bv*r_g*f*mu*mu*(1.0+fnl*Ak) - fnl*Ak*(1.0+fnl*Ak))*dendamp*dendamp*Pmm;
                    gsl_matrix_set(dPdt2, 0, 0, value);
                    value = vv_prefac*k*k*r_g*(1.0+fnl*Ak)*dendamp*Pmt;
                    gsl_matrix_set(dPdt2, 0, 1, value);
                    gsl_matrix_set(dPdt2, 1, 0, value);
                    break;
                // Differential w.r.t. zp_err
                case 8:
                    if (nbararray[0][i] > 0.0) {
                        value = 2.0*zp_err_v*zp_err_prefac/nbararray[0][i];
                        gsl_matrix_set(dPdt2, 1, 1, value);
                    }
                    break;
                default:
                    break;
            }
        }

        // If we are looking at the bias from the zero-point error, then once we have calculated dPdt2 we want the covariance matrix without
        // systematics and so without a zero-point offset. The same is true if we want the TRUE covariance matrix
        if (do_zp_bias) {
            if (((int)p[6] == 0) || ((int)p[6] == 1)) zp_err_v = 0.0;
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
                    error_obs = 100.0*error_dist*r;                                                  // Percentage error * distance * H0 in km/s (factor of 100.0 comes from hubble parameter)
                    error_noise = error_rand*error_rand + error_obs*error_obs + zp_err_v*zp_err_v;   // Error_noise is in km^{2}s^{-2}
                    n_u = nbararray[0][i]/error_noise;                   
                } else {
                    error_obs = 100.0*error_dist*r;                                                  // Percentage error * distance * H0 in km/s (factor of 100.0 comes from hubble parameter)
                    error_noise = error_rand*error_rand + error_obs*error_obs + zp_err_v*zp_err_v;   // Error_noise is in km^{2}s^{-2}
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
                //printf("%d, %lf, %lf, %lf, %lf, %lf\n", surv, k, mu, r, r_sum, surv_sum);

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
// ALSO MAKE NOTE OF THE FACTOR OF 1.0e-6 ON LINE 827. THIS IS BECAUSE I TYPICALLY SAVE THE VALUE OF NBAR x 10^6 IN THE INPUT FILES< SO THAT I DON'T LOSE PRECISION
// WHEN SMALL VALUES OF THE NUMBER DENSITY ARE WRITTEN TO A FILE!
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
    ezarray = (double *)calloc(NRED[0], sizeof(double));
    deltararray = (double *)calloc(NRED[0], sizeof(double));
    growtharray = (double *)calloc(NRED[0], sizeof(double));

    for (i=0; i<NRED[0]-1; i++) {
        zarray[i] = (zinarray[i+1]+zinarray[i])/2.0;
        rarray[i] = rz(zarray[i]);
        ezarray[i] = 1.0/ezinv(zarray[i], NULL);
        deltararray[i] = rz(zinarray[i+1]) - rz(zinarray[i]);
        growtharray[i] = growthz(zarray[i])/growthz(0.0);
        //printf("%12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n", zarray[i], rarray[i], deltararray[i], growtharray[i], nbararray[0][i], nbararray[1][i]);
    }
    zarray[NRED[0]-1] = (zmax+zinarray[NRED[0]-1])/2.0;
    rarray[NRED[0]-1] = rz(zarray[NRED[0]-1]);
    ezarray[NRED[0]-1] = 1.0/ezinv(zarray[NRED[0]-1], NULL);
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
        printf("ERROR: The maximum k in the input power spectra is less than k_max\n");
        exit(0);
    }

    return;
}

// Routine to read in the transfer function.
void read_transfer() {
    
    FILE * fp;
    char buf[500];
    int i;

    if(!(fp = fopen(Trans_file, "r"))) {
        printf("\nERROR: Can't open transfer function file '%s'.\n\n", Trans_file);
        exit(0);
    }

    int NTK = 0;
    while(fgets(buf,500,fp)) {
        if(strncmp(buf,"#",1)!=0) {
            double tk, ttk;
            if(sscanf(buf, "%lf %lf\n", &tk, &ttk) != 2) {printf("Transfer read error\n"); exit(0);};
            NTK++;
        }
    }
    fclose(fp);

    double * tkarray = (double *)calloc(NTK, sizeof(double));
    double * ttkarray = (double *)calloc(NTK, sizeof(double));

    NTK = 0;
    fp = fopen(Trans_file, "r");
    while(fgets(buf,500,fp)) {
        if(strncmp(buf,"#",1)!=0) {
            double tk, ttk;
            if(sscanf(buf, "%lf %lf\n", &tk, &ttk) != 2) {printf("Transfer read error\n"); exit(0);};
            tkarray[NTK] = tk;
            ttkarray[NTK] = ttk;
            NTK++;
        }
    }
    fclose(fp);

    // Spline the transfer function
    Trans_acc    = gsl_interp_accel_alloc();
    Trans_spline = gsl_spline_alloc(gsl_interp_cspline, NTK);
    gsl_spline_init(Trans_spline, tkarray, ttkarray, NTK);

    free(tkarray);
    free(ttkarray);

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



