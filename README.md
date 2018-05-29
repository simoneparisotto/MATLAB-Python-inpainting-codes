# MATLAB codes for the Image Inpainting Problem
This is a detailed Matlab implementation of classic inpainting methods.

All the scripts provided are used in Partial Differential Equation Methods for Image Inpainting (Carola-Bibiane Schoenlieb, Cambridge University Press, 2015):

```
@book{Schonlieb:2015ux,
 author    = {Sch\"{o}nlieb, Carola-Bibiane},
 title     = {{Partial Differential Equation Methods for Image Inpainting}},
 publisher = {Cambridge University Press},
 month     = {November}
 year      = {2015},
}
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

Copyright (c) 2016, Simone Parisotto and Carola-Bibiane Schoenlieb
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

<h4>1) Absolute Minimizing Lipschitz Extension Inpainting (AMLE)</h4>

**Path**: ./amle/inpainting_amle.m

**Note**: Used to reproduce Figure 4.10 in \cite{Schonlieb:2015ux}. (Only grayscale images).
      
**Bibliography**:
- Caselles, V., Morel, J. M., & Sbert, C. (1998). An axiomatic approach to image interpolation. Image Processing, IEEE Transactions on, 7(3), 376-386.
- Almansa, A. (2002). Echantillonnage, interpolation et détection: applications en imagerie satellitaire (Doctoral dissertation, Cachan, Ecole normale supérieure).

<h4>2) Harmonic Inpainting</h4>

**Path**: ./harmonic/inpainting_harmonic.m

**Note** Used to reproduce Figure 2.2 in \cite{Schonlieb:2015ux}. This program solves harmonic inpainting via a discrete heat flow. (Both Grayscale / Color Images).

**Bibliography**:
- Shen, J., & Chan, T. F. (2002). Mathematical models for local nontexture inpaintings. SIAM Journal on Applied Mathematics, 62(3), 1019-1043.

<h4>3) Mumford-Shah Inpainting with Ambrosio-Tortorelli approximation</h4>

**Path**: ./mumford-shah/inpainting_mumford-shah.m

**Note**: Used to reproduce Figure 7.3 in \cite{Schonlieb:2015ux}.  (Both Grayscale / Color Images).

**Bibliography**: 
- Esedoglu, S., & Shen, J. (2002). Digital inpainting based on the Mumford–Shah–Euler image model. European Journal of Applied Mathematics, 13(04), 353-370.


<h4>4) Cahn-Hilliard Inpainting</h4>

**Path**: ./cahn-hilliard/inpainting_cahn_hilliard.m

**Note**: Used to reproduce Figure 5.9 in \cite{Schonlieb:2015ux}.  (Both Grayscale / Color Images).

**Bibliography**: 
- Bertozzi, A., Esedoglu, S. & Gillette, A. (2007). Inpainting of binary images using the Cahn-Hilliard equation, IEEE Transactions on image processing 16.1 pp. 285-291 (2007).
- Schoenlieb, C.-B. & Bertozzi, A. (2011). Unconditionally stable schemes for higher order inpainting, Communications in Mathematical Sciences, Volume 9, Issue 2, pp. 413-457 (2011).


<h4>5) Transport Inpainting</h4>

**Path**: ./transport/inpainting_transport.m

**Note**: Refer to Section 6.1 in \cite{Schonlieb:2015ux}. (Both Grayscale / Color Images).

**Bibliography**:
- Bertalmio, M. (2001). Processing of flat and non-flat image information on arbitrary manifolds using partial differential equations.PhD Thesis, 2001.
