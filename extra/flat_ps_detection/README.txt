This directory is almost identical to icwt_v0.1 with the following changes:

    1. Different arrays for translations b and scales a
    2. Scales a are now logarithmic, affecting how we integrate for reconstruction
    3. Samples skymap is now over a different region of space
    
    Issues: Reconstruction could be significantly faster if we vectorize the sums
            Slowed down by the starting-up/shutting down of nodes
            
Legend:
    sample_skymap_upper - Sample skymap on upper region of sky
    sample_skymap_upper_2 - Same as above, but computes wavelet coefficients over more
                            translation parameters
    ..._upper_wps - Same as upper_2 but with injected point sources over a grid
                    Intensity varies with x-position
    ..._upper_wps_2 - Same as upper_2 but injected point sources at random positions 
                      with uniform intensity ; locations are saved for efficacy determination