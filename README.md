# MATLAB/Python Codes for the Image Inpainting Problem
[![DOI](https://zenodo.org/badge/73518110.svg)](https://zenodo.org/badge/latestdoi/73518110) [![View MATLAB/Python Codes for the Image Inpainting Problem on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/55326-matlab-python-codes-for-the-image-inpainting-problem)

This is a detailed MATLAB/Python implementation of classic inpainting methods.

All the scripts provided are used in the book "[Partial Differential Equation Methods for Image Inpainting](https://www.cambridge.org/core/books/partial-differential-equation-methods-for-image-inpainting/F4750CEA96C35354A97E2161130E91DC)" by Carola-Bibiane Schönlieb, Cambridge Monographs on Applied and Computational Mathematics, Cambridge University Press, 2015:

Please use the following entry to cite this code:

```
 @software{ParSch2016,
  author       = {Simone Parisotto and
                  Carola-Bibiane Sch\"{o}nlieb},
  title        = {{MATLAB/Python Codes for the Image Inpainting 
                   Problem}},
  month        = dec,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {3.0.0},
  doi          = {10.5281/zenodo.4315173},
  url          = {https://doi.org/10.5281/zenodo.4315173}
}
```

<h4>1) Absolute Minimizing Lipschitz Extension Inpainting (AMLE)</h4>

Used to reproduce Figure 4.10 in the book.

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/dataset/amle_clean.png"> <img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/results/amle_output.png"> 
      
**Bibliography**:
- Caselles, V., Morel, J. M., & Sbert, C. (1998). An axiomatic approach to image interpolation. Image Processing, IEEE Transactions on, 7(3), 376-386.
- Almansa, A. (2002). Echantillonnage, interpolation et detection: applications en imagerie satellitaire (Doctoral dissertation, Cachan, Ecole normale superieure).

<h4>2) Harmonic Inpainting</h4>

Used to reproduce Figure 2.2 in the book. This program solves harmonic inpainting via a discrete heat flow.

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/dataset/harmonic_input.png"> <img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/results/harmonic_output.png"> 

**Bibliography**:
- Shen, J., & Chan, T. F. (2002). Mathematical models for local nontexture inpaintings. SIAM Journal on Applied Mathematics, 62(3), 1019-1043.

<h4>3) Mumford-Shah Inpainting with Ambrosio-Tortorelli approximation</h4>

**Note**: Used to reproduce Figure 7.3 in the book.

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/dataset/mumford_shah_input.png"  width=45%> <img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/results/mumford_shah_output.png"  width=45%> 

**Bibliography**: 
- Esedoglu, S., & Shen, J. (2002). Digital inpainting based on the Mumford-Shah-Euler image model. European Journal of Applied Mathematics, 13(04), 353-370.

<h4>4) Cahn-Hilliard Inpainting</h4>

Used to reproduce Figure 5.9 in the book.

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/dataset/cahn_hilliard_input.png"> <img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/results/cahn_hilliard_output.png"> 

**Bibliography**: 
- Bertozzi, A., Esedoglu, S. & Gillette, A. (2007). Inpainting of binary images using the Cahn-Hilliard equation, IEEE Transactions on image processing 16.1 pp. 285-291 (2007).
- Schoenlieb, C.-B. & Bertozzi, A. (2011). Unconditionally stable schemes for higher order inpainting, Communications in Mathematical Sciences, Volume 9, Issue 2, pp. 413-457 (2011).

<h4>5) Transport Inpainting</h4>

Refer to Section 6.1 in the book. (Both Grayscale / Color Images).

<img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/dataset/transport_input.png"  width=45%> <img src="https://raw.githubusercontent.com/simoneparisotto/MATLAB-Python-inpainting-codes/master/matlab/results/transport_output.png" width=45%> 

**Bibliography**:
- Bertalmio, M. (2001). Processing of flat and non-flat image information on arbitrary manifolds using partial differential equations.PhD Thesis, 2001.

### License
[BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause)
