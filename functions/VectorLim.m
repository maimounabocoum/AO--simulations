% this function takes a vector nand return the same vector wihtin the limit
% imposed by the user : 

function [Vout Iout] = VectorLim(Vmin,Vmax,Vin)


I1 = find(Vin >= Vmin);
MinV = Vin(I1);

I2 = find(MinV <= Vmax);
Vout = MinV(I2);


Iout = I1(I2);







end