# PV_fisher
A short code to forecast errors on the normalised growth rate of structure for separate (partially overlapping) redshift and peculiar velocity surveys. Also works for just one of these, so can be used to forecast your favourite LSS survey. 

This code uses the Fisher Matrix method to forecast the constraints on the normalised growth rate of structure, fsigma8 from combined redshift and peculiar velocity surveys. It does this by considering the information in the anisotropic two-point functions of the density and velocity fields measured from these two surveys, but does not require the two surveys to be fully (or even partially) overlapping.

The required inputs are a bunch of cosmological and survey parameters, the number density of galaxies (x10^6 h^3 Mpc^-3) in the peculiar velocity and redshift surveys and the input matter-matter, matter-velocity and velocity-velocity power spectra at suitable redshifts. These are hardcoded (Ugh, I know!) at the top of the code, and hopefully the exaplanations are sufficient. The input spectra can be calculated using a variety of codes. My method of choice is the implementation of 2nd order Renormalised Perturbation Theory found in COPTER (http://mwhite.berkeley.edu/Copter/), but the simplest method is just to use a linear CAMB (http://camb.info/) matter power spectrum for all three (the three dark matter power spectra are expected to agree on linear scales anyway).

To compile the code you'll need a C compiler (I've only used gcc, but it's a pretty simple code) and the GSL libraries. With these in hand, compilation is as simple as:

    gcc PV_fisher.c -o PV_fisher -lm -lgsl -lgslcblas
    
A set of example files you can use to test the code are given in the /example_files/ directory.
    
This code has been used by me for the papers: Howlett, Staveley-Smith and Blake, 2017 (http://adsabs.harvard.edu/abs/2017MNRAS.464.2517H); Da Cunha et. al., 2017 (http://adsabs.harvard.edu/abs/2017arXiv170601246D); Howlett et. al 2017 (.

The maths is presented mainly in the first of these and if you use the code could you please cite that paper. However, the code given here is not the same as used for the first two (I added the redshift dependence later) and so you won't be able to fully reproduce those results (although they should be reasonably close). However, if you use the files given in the /example_files/ directory and the default parameters in the code, you should have no problem reproducing the last column in Table 1 of the third paper.

Finally, if you read the first paper you'll notice that there are a bunch of things in there that are absent in PV_fisher.c: namely, the ability to handle two distinct redshift and two PV surveys with different parameters (so four datasets in total) and forecasts on primordial non-Gaussianity, scale-dependent galaxy and velocity bias, and zero-point offsets. In have included primordial non-Gaussianity, scale-dependence and zero-points, and a calculation of the systematic bias associated with them, in PV_fisher_extended.c. It makes the code quite a bit more confusing, so it's recommended you only use this if you know what you want to do, or contact me first.


