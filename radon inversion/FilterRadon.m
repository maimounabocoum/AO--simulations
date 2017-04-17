function FILTER = FilterRadon(f, N , filter , Fc)
% Function written by Maïmouna Bocoum - 28/02/2017 as fully inspired by
% iradon function from matlab packages (type help iradon for more details)

FILTER = abs(f) ;

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

 FILTER(abs(f) > Fc ) = 0;                         % Crop the frequency response
 FILTER = FILTER.*exp(-2*abs(f/Fc).^10);
end

