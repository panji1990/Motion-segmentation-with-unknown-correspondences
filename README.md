# Motion segmentation with unknown point correspondences

Code for the following paper:

Pan Ji, Hongdong Li, Mathieu Salzmann, and Yuchao Dai. Robust motion segmentation with unknown correspondences. In Europian Conference on Computer Vision (ECCV), 2014. Poster.



1. Install ncuts from ~~https://www.cis.upenn.edu/~jshi/software/~~ https://github.com/panji1990/Ncut9.
2. Install vlfeat from http://www.vlfeat.org/.
3. Download Hopkins 155 outdoor sequences from https://www.dropbox.com/s/qzsw958uhk8r9j1/HopkinsOutdoor.zip?dl=0.
4. Change file path where necessary.
5. Compile the Hungarian algorithm (C code) with your Matlab, e.g. mex assignmentoptimal.c.
6. Install mmread from http://au.mathworks.com/matlabcentral/fileexchange/8028-mmread.
