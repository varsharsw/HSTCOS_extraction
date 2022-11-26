# HSTCOS_extraction
Extract target spectra from COS aperture


Crowding in a COS aperture is common when observing targets at greater distances. The Fortran code cosextract.f extracts target spectra from the contaminated 2D flat-fielded detector images. Here we simultaneously fit for both our target and neighboring cluster  by multiplying them with the (wavelength-dependent) cross-dispersion instrumental profile and subsequently subtract from the combined observed spectra. This has been done along the cross-dispersion axis for all wavelengths. Finally, we extracted two spectral components such that the residual is minimum. For more details see Ramachandran et al. 2022 
