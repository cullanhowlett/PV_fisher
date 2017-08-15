# PV_fisher
A short code to forecast errors on the normalised growth rate of structure for separate (partially overlapping) redshift and peculiar velocity surveys. Also works for just one of these, so can be used to forecast your favourite LSS survey. 

This code uses the Fisher Matrix method to forecast the constraints on the normalised growth rate of structure, fsigma8 from combined redshift and peculiar velocity surveys. It does this by considering the information in the anisotropic two-point functions of the density and velocity fields measured from these two surveys, but does not require the two surveys to be fully (or even partially) overlapping.

The required inputs are a bunch of cosmological and survey parameters, the number density of galaxies in the peculiar velocity and redshift surveys and the input matter-matter, matter-velocity and velocity-velocity power spectra at suitable redshifts. These spectra can be calculated using a variety of codes. My method of choice is the implementation of 2nd order Renormalised Perturbation Theory found in COPTER (Carlson 2009), but the simplest method is just to use a linear CAMB power spectrum for all three (the three dark matter power spectra are xpected to agree on linear scales anyway).

To compile the code you'll need a C compiler (I've only tested on gcc, but it's a pretty simple code) and the GSL libraries. With these in hand, compilation is as simple as:
    gcc PV_fisher.c -o PV_fisher -lm -lgsl -lgslcblas
    
This code is based on the papers Howlett, Blake and Staveley-Smith, 2017 and Howlett et al., (in preparation). The maths is presented in there and if you use it could you please cite the first of these. If you read this paper you'll notice that there are a bunch of things in there that are absent in this code: namely, the ability to handle two distinct redshift and two PV surveys with different parameters (so four datasets in total) and forecasts on primordial non-Gaussianity, scale-dependent galaxy and velocity bias, and zero-point offsets. I took these out because I wanted the code to be more user-friendly (read "actually usable") and figured a simple version would work for most cases. If you have need of these features, let me know :)


