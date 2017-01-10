function [ratio]=Disk_Relative_Energy(Im)

%====================================================================
% This function computes the 2DFFT for the given image, and compute the relative energy 
% outside the radious pi disk. Note that for white noise we expect to get that the ratio is 
% 1-pi/4=0.2146 because the disk area is pi^3 while the entire area is 4pi^2. For LPF signals,
% this ratio is far smaller.
% 
% Synopsis: [ratio]=Disk_Relative_Energy(Im)
% 
% Input:    Im - an input image
% Output: ratio - the obtained ratio
%
% Example:
%    A=imread('Example_Image1.jpg ');
%    r=Disk_Relative_Energy(A);
%    disp(r);
%
% Written by Miki Elad on March 20th, 2005.
%====================================================================

[N,M]=size(Im);
Im=double(Im);
Im=fft2(Im);
Im=abs(Im);
Im=fftshift(Im);

[xx,yy]=meshgrid(1:1:M,1:1:N);
Disk=(xx*2*pi/max(xx(:))-pi).^2+(yy*2*pi/max(yy(:))-pi).^2>pi^2;

ratio=sum(Im(Disk).^2)/sum(Im(:).^2);

return;
