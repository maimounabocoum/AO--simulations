% scrSignal = [];


 Signal = []; 
% Noise = [];

for i=1:100
    SimulateCameraHolography;
    Signal = [Signal SIG]
    % Noise = [Noise NOISE]
end

%  snr = mean(Signal - Noise)/std(Signal) 

% 1024px : 792.0663
% 512px :  310.4991
% no signal 512 px : -0.1820
