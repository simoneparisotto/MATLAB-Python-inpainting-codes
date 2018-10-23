# MATLAB codes for the Image Inpainting Problem
This is a detailed Matlab implementation of classic inpainting methods.

All the scripts provided are used in

```
Schoenlieb, Carola-Bibiane
Partial Differential Equation Methods for Image Inpainting.
Cambridge Monographs on Applied and Computational Mathematics,
Cambridge University Press, 2015
doi:10.1017/CBO9780511734304
```

Please use the following entry to cite this code:

```
@Misc{MATLABinpainting2016,
 author       = {Parisotto, Simone and Sch\"{o}nlieb, Carola},
 title        = {MATLAB Codes for the {Image} {Inpainting} {Problem}},
 howpublished = {GitHub repository, {MATLAB} Central File Exchange},
 month        = {September},
 year         = {2016}
  }
```

License: BSD-3-Clause (https://opensource.org/licenses/BSD-3-Clause)

<h4>1) Absolute Minimizing Lipschitz Extension Inpainting (AMLE)</h4>

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/dataset/amle_input.png" width="150px"> 
<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/results/amle_output.png" width="150px"> 

https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/u_CLA_65.png

**Path**: ./amle/inpainting_amle.m

**Note**: Used to reproduce Figure 4.10 in the book. (Only grayscale images).
      
**Bibliography**:
- Caselles, V., Morel, J. M., & Sbert, C. (1998). An axiomatic approach to image interpolation. Image Processing, IEEE Transactions on, 7(3), 376-386.
- Almansa, A. (2002). Echantillonnage, interpolation et detection: applications en imagerie satellitaire (Doctoral dissertation, Cachan, Ecole normale superieure).

<h4>2) Harmonic Inpainting</h4>

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/dataset/harmonic_input.png" width="150px"> 
<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/results/harmonic_output.png" width="150px"> 

**Path**: ./harmonic/inpainting_harmonic.m

**Note** Used to reproduce Figure 2.2 in the book. This program solves harmonic inpainting via a discrete heat flow. (Both Grayscale / Color Images).

**Bibliography**:
- Shen, J., & Chan, T. F. (2002). Mathematical models for local nontexture inpaintings. SIAM Journal on Applied Mathematics, 62(3), 1019-1043.

<h4>3) Mumford-Shah Inpainting with Ambrosio-Tortorelli approximation</h4>

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/dataset/mumford_shah_input.png" width="150px"> 
<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/results/mumford_shah_output.png" width="150px"> 

**Path**: ./mumford-shah/inpainting_mumford-shah.m

**Note**: Used to reproduce Figure 7.3 in the book.  (Both Grayscale / Color Images).

**Bibliography**: 
- Esedoglu, S., & Shen, J. (2002). Digital inpainting based on the Mumford-Shah-Euler image model. European Journal of Applied Mathematics, 13(04), 353-370.


<h4>4) Cahn-Hilliard Inpainting</h4>

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/dataset/cahn-hilliard_input.png" width="150px"> 
<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/results/cahn-hilliard_output.png" width="150px"> 

**Path**: ./cahn-hilliard/inpainting_cahn_hilliard.m

**Note**: Used to reproduce Figure 5.9 in the book.  (Both Grayscale / Color Images).

**Bibliography**: 
- Bertozzi, A., Esedoglu, S. & Gillette, A. (2007). Inpainting of binary images using the Cahn-Hilliard equation, IEEE Transactions on image processing 16.1 pp. 285-291 (2007).
- Schoenlieb, C.-B. & Bertozzi, A. (2011). Unconditionally stable schemes for higher order inpainting, Communications in Mathematical Sciences, Volume 9, Issue 2, pp. 413-457 (2011).


<h4>5) Transport Inpainting</h4>

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/dataset/transport_input.png" width="150px"> 
<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Codes-for-the-Image-Inpainting-Problem/blob/master/results/transport_output.png" width="150px"> 

**Path**: ./transport/inpainting_transport.m

**Note**: Refer to Section 6.1 in the book. (Both Grayscale / Color Images).

**Bibliography**:
- Bertalmio, M. (2001). Processing of flat and non-flat image information on arbitrary manifolds using partial differential equations.PhD Thesis, 2001.
