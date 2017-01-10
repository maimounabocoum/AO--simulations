%====================================================================
%                                                           Content of this directory
%====================================================================
% 
% ----------------------------------------------------------------------------------------------------------------------------         
% 1. Grid operations:
% ----------------------------------------------------------------------------------------------------------------------------         
% Create_Grid - This function creates a set of points in the plane froming either Cartesian, 
% Polar, or Recto-Polar grids (several variations). A set of parameters dictate the span of 
% these points, their number, and more. 
% 
% Grid_Evolution - describes the changes we perform when moving from the pseudo-polar 
% grid to the equi-spaced angles pseudo-polar grid, and from there to the Polar grid.
% 
% Distinct - This function counts the number of distict and different ROWS in Array. It is 
% good in order to find how many repetitions we have in a created grid 
%              
% Create_Oversampled_Grid - creates the 'D', 'S' and the 'X' grids such that they are 
% oversampled. This function is similar to 'Create_Grid.m' but uses different 
% parametrization. It is used as part of the effort to perform the interpolation along rays 
% more accurate, and debug the developed algorithm.
%
% ----------------------------------------------------------------------------------------------------------------------------         
% 2. Exact and Approximate Fourier Transforms:
% ----------------------------------------------------------------------------------------------------------------------------         
% Brute_Force_Transform - for a 2D discrete signal given on a Cartesian grid, it computes 
% the Discrete Fourier Transform on an arbitrary required grid defined by [GridX,GridY].
% 
% Transform_Matrix  - obtains a 2D discrete signal given on a cartesian grid, and computes 
% for it the Discrete Fourier Transform matrix on a required grid defined by [GridX,GridY].
% 
% RectoPolar_Trans_New - performs a Recto (Pseudo) Polar transform on a 2D signal X  
% given on a Cartesian grid. If X is N*N, the output will have 2N by 2N output values. The 
% algorithm applied here uses the Fast Fractional Fourier Transfrom. The transform is 
% applied in two pairs of quadrants 1&3 and then 2&4.
% 
% Pseudo_Polar_Evolution - developing an evolution of stages from the brute-force grid 
% based Transform and to the fast Pseudo-Polar FT. The target is getting a fast 
% PseudoPolar transform that will match perfectly the results obtained by brute force. The 
% tested methods are: 
%    1. Brute force transform  with no vector multiplications - (like in C) 
%    2. Brute force transform but with vector multiplications (exploiting Matlab capabilities)
%    3. First stage using regular 1D FFT and second stage using brute-force summations
%    4. First stage using regular 1D FFT and second stage using FR-FFT
% 
% RectoPolar_Trans_Adjoint - performs a Recto (Pseudo) Polar ADJOINT transform on 
% a 2D signal X given on the PP-coordinate system. If X is 2N*2N, the output will have N 
%  by N values. The algorithm applied here uses the Fast Fractional Fourier Transform. 
% This program proposes a set of algorithms (evolution) from the direct approach to the 
% fastest approach. 
%
% PPFFT - Pseudo-polar fast Fourier transform - very similar to 'RectoPolar_Trans_New'
% but cleaner, with ability to obtain oversampling on both axes, and ability to transform 
% arbitrary sized input.
% 
% APPFFT - Adjoint Pseudo-polar fast Fourier transform - very similar to 
% 'RectoPolar_Trans_Adjoint ' but cleaner, with ability to handle oversampling on both 
% axes, and ability to transform arbitrary sized input.
% 
% IPPFFT - Inverse pseudo-polar fast Fourier transform - using a pre-conditioned iterated 
% LS solver that is based on using the PPFFT and the APPFFT. 
%
% AFTUSF_NGP_1 - approximate Fourier Transform at Unequally Spaced Frequencies. 
% The approximation is built on bilinear interpolation applied on an oversampled grid. This is 
% taken from Dave Donoho's code.
%             
% FTUSF_Spline_1 - approximate Fourier Transform at Unequally Spaced Frequencies.  
% The approximation is built on spline interpolation with x- and y-derivatives applied on an 
% oversampled grid. This is taken from Dave Donoho's code.
%                        
% Cartesian_2_Polar - for a 2D discrete signal on a Cartesian grid, it applies Polar-FT using 
% two methods: 
%    1. Direct transform that may take long, and 
%    2. Computing an oversampled 2D Cartesian FFT and bilinear interpolating.         
%
% SPolar_Transform - performs a Fourier transform on a pseudo-polar grid with equaly 
% spaced angles (and not distances as regular Pseudo-Polar). several methods 
% are tested:
%    1. exact
%    2. oversampling and bilinear interpolation
%    3. oversampling and spline interpolation
%    4. oversampling, computing derivatives, and spline interpolation
%           
% XPolar_Transform - performs a Fourier transform on a polar grid by using the S-Polar 
% transform first, and then  equalizing distances along rays. 
%  
% Polar_Transform_New - performs a Fourier Transfrom over the Polar grid, using either 
% spline or Hermite interpolation for the two stages. As opposed to the XPolar_Transform 
% routine, this code starts by regular Pseudo-Polar transform with oversampling by zero 
% padding. Then the interpolations are done by working on the rows and the columns of the 
% produced array. A shortcoming of this algorithm is that the oversampling on the two 
% axes (r,theta) is forced to be the same.
%       
% PFFT - Fast polar Fourier Transform  - very much similar to
% XPolar_Transform, but removes the choice of some parameters. 
%       
% IPFFT - Fast inverse polar Fourier Transform ??????? NOT DONE YET
%
% ----------------------------------------------------------------------------------------------------------------------------         
% 3. Tools for the above Transforms:
% ----------------------------------------------------------------------------------------------------------------------------         
% My_FRFT - My Fractional Fourier Transform, computing the transform in a fast manner.
%            
% My_FRFT_Centered - My Fractional Fourier Transform, computing the transform as 
% before but for centered indices.
% 
% Test_FRFT - analyse the proper way to perform the Fractional Fourier Transform. 
% This part assumes that the output is required for n=0,1,2, .... ,N-1.
%             
% Bilinear_Interp_Matrix - builds a matrix that takes an over-sampled uniform grid (N 
% mulitplied by R), and performs the bilinear interpolation of the values of the FFT for 
% the destination grid [Xd,Yd].
%             
% Compare_Bilinear_Interpolations - compares between mine and Donoho bilinear 
% interpolation, where Donoho applies the interpolation on a signal, and we construct 
% a matrix.
%  
% interp_SPLINE - 1-D interpolation using function values and its derivatives. This 
% function interpolates to find YY, the values of the underlying function Y at the points 
% in the vector XX. The vector X specifies the points at which the data Y is given. 
%
% Debug_interp_SPLINE - Debugging the function interp_SPLINE
%
% Debug_XPolar_FFT - Tracing the origins of the errors in the polar FFT - this is a 
% debugging code. 
% 
% ----------------------------------------------------------------------------------------------------------------------------         
% 4. Worst case Error Analysis as eigenvalue problems:
% ----------------------------------------------------------------------------------------------------------------------------         
% Error_Analysis_1 - performs an error analysis for the Unequaly Spaced Fourier 
% Transform. It finds the worst signal that maximizes the error, under the assumption that 
% the destination grid is cartesian. 
% 
% Error_Analysis_2 - same as Error_Analysis_1 but the destination grid is Polar. 
% 
% Error_Analysis_3 - performs an error analysis for the Polar Fourier transform, computed 
% based on Recto-Polar coordinates. The difference from Error_Analysis_2.m is that 
% the RectoPolar FT are computed exactly (assuming we have a fast algorithm for 
% their computation), and then we interpolate to get the Polar points. The interpolation 
% used is the simplest possible - Nearest-Neighbor.
%
% Ray_Behavior - This function displays one ray of the FFT function in order to show that 
% it is indeed band limitted. This builds the justification for the polynomial interpolation. 
% The program prompts the user to choose a ray and then shows the sampling of its 
% values with growing sampling rate, while showing the accuracy obtained using a 
% polynomial approximation (spline).
%
% ----------------------------------------------------------------------------------------------------------------------------         
% 5. General tools:
% ----------------------------------------------------------------------------------------------------------------------------         
% Disk_Relative_Energy - by the 2D-FFT for the given image it computes the relative  
% energy outside the radious pi disk. Note that for white noise we expect to get that the  
% ratio is 1-pi/4=0.2146 because the disk area is pi^3 while the entire area is 4pi^2. For 
% LPF signals, this ratio is far smaller.
% 
% Distance_from_Grid - This function computes the distance between the polar and its
% approximating grids (caresian and pseudo-polar). The distance is computed as a  
% weighted Euclidean distance, matching for every polar coordinate its nearest neighbor  
% in the approximating grid. The weights are chosen based on the relative energy that  
% exists in each polar coordinate. 
%
% ----------------------------------------------------------------------------------------------------------------------------         
% 6. Scripts to produce the paper/slides figures and results:
% ----------------------------------------------------------------------------------------------------------------------------         
% Comparison_1_USFT_vs_Polar - comparing both the NGP and the spline USFT
% methods to the direct polar FFT using the pseudo-polar approach. The
% comparison is done for specific signals, and by sweeping the oversampling
% parameter.
%
% Comparison_2_USFT_vs_Polar - comparing both the NGP and the spline USFT
% methods to the direct polar FFT using the pseudo-polar approach. The
% comparison is done by building the transform matrices and then analysing
% these matrices from the eigenspace point of view.
%
% Comparison_3_USFT_vs_Polar - comparing both the NGP and the spline USFT
% methods to the direct polar FFT using the pseudo-polar approach. The
% comparison's objective is to show that the error in the frequency domain is 
% spatially dependent. We show this by first implimenting the
% involved transforms on example images, and later by
% worst-case scenario analysis.
%
% Comparison_4_Polar_Corner - computing the condition number of the polar
% FFT matrix to the same one with augmentation of zeros in the corners. The idea is to 
% show that the condition number is very stable and low for such an additional constraint. 
%
% Create_Figures_For_paper - this script generates the figure that appears
% in the paper (in its longer version).
%====================================================================
