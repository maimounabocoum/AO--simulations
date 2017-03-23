function FILTER = FilterRadon(f, N , filter , Fc)
% Function written by Maïmouna Bocoum - 28/02/2017 as fully inspired by
% iradon function from matlab packages (type help iradon for more details)
%   
% N : size of fourier transform (N = 2^k) used for radon fourier transform
% order = 2^nextpow2(N) ;
% order = 2*max(64,order)
N = N/2 ;
tau = 0.1*Fc;

n = 0:N; 
filtImpResp = zeros(1,N+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/(4*tau^2); % Set the DC term 
filtImpResp(2:2:end) = -1./((tau*pi*n(2:2:end)).^2); % Set the values for odd n

% symmetry of filter around 0 axis :
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
            
%FILTER = 2*real(fft(filtImpResp)); 
FILTER = abs(f).^0.9 ;
% figure
% plot(FILTER,'r')
% 
% FILTER = abs(f) ;
% hold on
% plot(FILTER,'o')

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        FILTER = FILTER.* (sin(2*pi*f/(2*Fc))./((2*pi*f)/(4*Fc)));
        
    case 'cosine'
        %FILTER(2:end) = FILTER(2:end) .* cos(2*pi*f(2:end)/(0.5*2*Fc));
        FILTER = FILTER .* cos(2*pi*f/(4*Fc));
    case 'hamming'
        FILTER = FILTER.* (.54 + .46 * cos(2*pi*f/Fc));
    case 'hann'
        FILTER = FILTER.*(1+cos(2*pi*f./(4*Fc))) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end

FILTER(isnan(FILTER)) = 0 ;



% f = (0:size(filt,2)-1)/order;   % frequency axis up to Nyquist


end

