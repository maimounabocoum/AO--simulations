function [ Ifft, Iout ] = GetTagged( Fstruct , Iin , Hfilter )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Ifft = Fstruct.fourier(Iin);
Iout = Fstruct.ifourier( Ifft.*Hfilter );



end

