# extlaw_H17
A python code to generate the extinction law between 0.8 -- 2.2 microns as described in Hosek+17. The law is derived using the Westerlund 1 (Wd1) main sequence (A_Ks ~ 0.6 mag) and Arches cluster field Red Clump at the Galactic Center (A_Ks ~ 2.7 mag). To derive the law a Wd1 cluster age of 5 Myr is assumed, though changing the cluster age between 4 Myr -- 7 Myr has no effect on the law.

This extinction law can be applied to highly reddened stellar populations that have similar foreground material as Wd1 and the Arches RC, namely dust from the spiral arms of the Milky Way in the Galactic Plane. 


#### Dependencies
* python (either 2.7 or higher)
* numpy
* scipy
* pylab

#### Functions 
extinction(AKs, wavelength): returns the total extinction in magnitudes (A_lambda) and corresponding error at the specified wavelength(s) and overall A_Ks. The wavelengths must be in microns and A_Ks in magnitudes. Input wavelength(s) can either be a float or array of floats

plot_extinction_law(): Plot the extinction law (with errors) between 0.8 -- 2.2 microns. Saves output plot as extlaw_H17.png in working directory.

###### Example
To get the total extinction at 1.25 microns and A_Ks = 0.5 mags (in an ipython session):

> import extlaw_H17 # Import code

> result, err = extlaw_H17.extinction(0.5, 1.25)

> result[0] # extinction at 1.25 microns

> result[0] + err[0] # extinction + 1-sigma err 

> result[0] - err[0] # extinction - 1-sigma err





