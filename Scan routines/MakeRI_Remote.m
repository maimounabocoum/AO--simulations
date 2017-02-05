function RI = MakeRI_Remote(f0,fs,N)

LO = round(fs/f0);
RI = zeros(1,N*LO);
NbRI = length(RI);

for i=1:NbRI
   RI(i) = sin(2*pi*f0/fs*i)*exp(-(i-(NbRI-1)/2)^2/((NbRI-1)/0.5)^2);
   %RI(i) = cos(2*pi*f0/fs*(i-(NbRI-1)/2))*exp(-(i-(NbRI-1)/2)^2/((NbRI-1)/4)^2);
end
%figure; title('pulse d''emission');
%plot(RI);