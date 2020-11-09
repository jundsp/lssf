# Laplace State Space Filter

Matlab code that implements a state space filter for univariate Laplace-distributed data sequences, as detailed in the following paper:

J. Neri, P. Depalle, R. Badeau, "<a href="https://ieeexplore.ieee.org/document/9053185" target="_blank">**Laplace State Space Filter with Exact Inference and Moment Matching**</a>," *IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP)*, pp. 5880-5884, Barcelona, Spain, 2020.

The script `demo.m` runs a demonstration of the Bayesian filter. Three different test data sequences are available to choose from: Laplace noise, outliers, and noise switch.

It uses exact inference to infer the latent state sequence from temporal data. It successfully filters outliers and heavy-tailed noise, in addition to Laplace noise, Gaussian noise, and Cauchy noise. It is as fast as the Kalman filter.

<img src="https://www.music.mcgill.ca/~julian/wp-content/uploads/2020/11/lssf_outlier.png" width="40" height="40" />

