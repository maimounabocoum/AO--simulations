clear all
close all
clc
field_init(0);

% Set initial parameters (cf PhantomImage.m)
f0=5.5e6;           % Transducer center frequency [Hz]
fs=200e6;           % Sampling frequency [Hz]
c=1540;             % Speed of sound [m/s]
lambda=c/f0;        % Wavelength [m]

N_elements = 1;
focus = [0 0 40]/1000;
R_elem = 3.5/1000;
el_size = 0.4/1000;

Th = xdc_concave(R_elem, focus(3), el_size);

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
figure; plot(impulse)

excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (Th, impulse);
figure; plot(excitation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('phantom_wave.mat')
load('phantom_amplitudes.mat')

N_elements = 1;         % Number of Elements []
fs = 200e6;             % Sampling frequency [Hz]

for i=1:1:size(phantom_wave,3)
    i
    [v,t]= calc_scat_multi(Th,Th,phantom_wave(:,:,i),phantom_amplitudes);
    signal(:,i)=v;
    Temps(i,1)=t;
end

save(['signal'],'signal');
save(['Temps'],'Temps');

test=max(signal')-min(signal');
ind1 = find(test==max(test));
ind2 = find(test==min(test));

k=1
close all
for j=10*k-9:1:10*k%size(phantom_wave,3)
    v(:) = signal(:,j);
    [N,M]=size(v);
    %v=v/max(max(v));
    t = Temps(j,1);
    figure(j)
    %figure
    [N,M]=size(v);
    v=v/max(max(v));
    for i=1:N_elements
        plot(2*(0:N-1)/fs+t,v(:,i)+i), hold on
    end
    hold off
    title('Individual traces')
    xlabel('Time [s]')
    ylabel('Normalized response')
    %axis([t t+2*N/fs 0 M+1])
    %axis([t t+2*N/fs -1e-4 1e-4])
    %xlim([t t+2*N/fs])
%     v(:) = signal(:,j);
%     [N,M]=size(v);
%     v=v/max(max(v));
%     t = Temps(j,1);
%     figure(j)
%     if size(v,2)>1
%         plot(2*(0:N-1)/fs+t,sum(v'))
%     else
%         plot(2*(0:N-1)/fs+t,v(:,i)+i)
%     end
%     title('Summed response')
%     xlabel('Time [s]')
%     ylabel('Normalized response')
%     axis([t t+2*N/fs 0 M+1])
    w = v';
    %save(['Registered_Data'],'w');
    %ScatterInfo.fs = fs;
    %ScatterInfo.t  = t;
    disp('Data saved')
end


%close all; 
for i=1:1:size(phantom_wave,3)-1
c=xcorr(signal(:,i),signal(:,i+1)); cmax(i)=max(c);
%figure(i); plot(c); axis([0 1563 -1e-38 1e-38])
end
figure; plot([0.1*(2:1:(size(phantom_wave,3)))-0.1],cmax); xlim([1 10])
