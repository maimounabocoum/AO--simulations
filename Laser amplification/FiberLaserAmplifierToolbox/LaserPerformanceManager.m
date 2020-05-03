% Run analysis on various fiber laser configurations
function [data hFig] = LaserPerformanceManager(s)

%% INITIALIZATION %%

% Resolve Inputs
if nargin<1, s = struct; end % create dummy structure if not present
if ~isfield(s,'mode'), s.mode = 'reflectance'; end % plot vs reflectance, length, losses, or 3D vs refl/length
if ~isfield(s,'plotFlag'), s.plotFlag = 1; end
if ~isfield(s,'Pp'), s.Pp = 0.200; end

% Define variables against which to run calculations
switch lower(s.mode(1:2))
    case 're' % plot versus reflectance
        if ~isfield(s,'R'), 
            Rcalc = (1:1:99)/100; 
        else
            Rcalc = s.R; % allow input of a few discrete reflectances if desired
        end
    case 'le' % plot versus length
        if ~isfield(s,'L'), s.L = 10; end
        Lcalc = s.L*(5:5:150)/100; % vary fiber length from 5% to 150% of specified
    case 'lo' % plot versus losses
        if ~isfield(s,'eps'), s.eps = 10; end
        Epscalc = s.eps*(5:5:150)/100; % vary losses from 5% to 150% of specified
    case {'bo','3d','po'} % 3D plot of both versus power
        if ~isfield(s,'R'), 
            Rcalc = (5:5:95)/100;  % vary from 5% to 95%
        elseif length(s.R)>1
            Rcalc = s.R; 
        else
            Rcalc = s.R*(10:10:100)/100; % vary from 50% to 150% of given R value
            Rcalc = Rcalc(Rcalc<=1); % of course can't go over 1
        end
        if ~isfield(s,'L'), s.L = 10; end
        Lcalc = s.L*(10:10:100)/100; % vary from 25% to 150%
end
% Pump power
Pop = s.Pp;
Pcalc = Pop*(5:5:125)/100; % vary pump power from 5% to 125% of operating point

%% CALCULATION %%

switch lower(s.mode(1:2))
    case 're' % plot versus reflectance
        Psop = zeros(1,length(Rcalc));
        wlDisc = zeros(1,length(Rcalc));
        for rr = 1:length(Rcalc)
            s.R2 = Rcalc(rr);
            s.R1 = 1;
            s.R = sqrt(s.R2 * s.R1);
            [eta(rr), Ppth(rr) Grt(rr) s] = LaserPerformance(s);
            if eta(rr)
                Psop(rr) = eta(rr)*(Pop-Ppth(rr)); 
                wlDisc(rr) = CheckWavelengthDiscrimination(Ppth(rr),s);
            end
        end
    case 'le' % plot versus length
        Psop = zeros(1,length(Lcalc));
        wlDisc = zeros(1,length(Lcalc));
        for ll = 1:length(Lcalc)
            s.L = Lcalc(ll);
            [eta(ll), Ppth(ll) Grt(ll) s] = LaserPerformance(s);
            if eta(ll)
                Psop(ll) = eta(ll)*(Pop-Ppth(ll));
                wlDisc(ll) = CheckWavelengthDiscrimination(Ppth(ll),s);
            end
        end  
    case 'lo' % plot versus cavity losses
        Psop = zeros(1,length(Epscalc));
        wlDisc = zeros(1,length(Epscalc));
        for ll = 1:length(Epscalc)   
            s.eps = Epscalc(ll);
            s.eps2db = s.eps/3;
            s.eps1db = s.eps - s.eps2db;          
            [eta(ll), Ppth(ll) Grt(ll) s] = LaserPerformance(s);
            if eta(ll)
                Psop(ll) = eta(ll)*(Pop-Ppth(ll));
                wlDisc(ll) = CheckWavelengthDiscrimination(Ppth(ll),s);
            end
        end           
    case {'bo','3d','po'} % 3D plot of both versus power
        for rr = 1:length(Rcalc)
            s.R2 = Rcalc(rr);
            for ll = 1:length(Lcalc)
                s.L = Lcalc(ll);
                [eta(rr,ll), Ppth(rr,ll) Grt(rr,ll) s] = LaserPerformance(s);
                Psop(rr,ll) = eta(rr,ll)*(Pop-Ppth(rr,ll));
                wlDisc(rr,ll) = CheckWavelengthDiscrimination(Ppth(rr,ll),s);
            end
        end
end

%% PLOTS %%

hFig = figure('windowstyle','docked');

switch lower(s.mode(1:2))
    case 're' % plot versus reflectance        
        %   Plot efficiency and threshold versus reflectance
        [ax h1 h2] = plotyy(100*Rcalc,1000*Ppth,100*Rcalc,100*eta);
        axes(ax(1)); set(gca,'fontsize',16,'linewidth',2,'ytickmode','auto'); hold on; grid on;
        xlabel('Reflectance of Output Coupler (%)');
        ylabel('Power (mW)');
        ylim([0 1000*Pop]);
        axes(ax(2)); hold on; set(gca,'fontsize',16,'linewidth',2); grid on;
        ylabel(sprintf('Slope Efficiency (%%) and\n Wavelength Discrimination (dB)'));
        yR = get(ax(2),'YLim');
        nTick1 = numel(get(ax(1),'YTick'));
        set(ax(2),'YTick',linspace(0,yR(2),nTick1));
        h4 = plot(ax(2),100*Rcalc,wlDisc); 
        h3 = plot(ax(1),100*Rcalc,1000*Psop);        
        set(h1,'linewidth',2,'color','r','marker','o');
        set(h2,'linewidth',2,'color',[0 0.5 0],'marker','o');
        set(h3,'linewidth',2,'color','b','marker','o');
        set(h4,'linewidth',2,'color','k','marker','o');
        hL = legend([h2,h1,h3,h4],{'Efficiency','Threshold P_P','Operating P_{out}','Wavelength Discrim.'},'Location','Best'); set(hL,'Color','w');
        hT = title(sprintf('Laser Performance vs Reflectance: L=%dcm, Losses=%.1fdB',round(100*s.L),s.eps1db+s.eps2db));
        data = {eta,Ppth,Grt,Rcalc,s};
    case 'le' % plot versus length
        %   Plot efficiency and threshold versus active fiber length
        [ax h1 h2] = plotyy(100*Lcalc,1000*Ppth,100*Lcalc,100*eta);
        axes(ax(1)); set(gca,'fontsize',16,'linewidth',2,'ytickmode','auto'); hold on; grid on;
        xlabel('Length of Active Fiber (cm)');
        ylabel('Power (mW)');
        ylim([0 1000*Pop]);
        axes(ax(2)); hold on; set(gca,'fontsize',16,'linewidth',2); grid on;
        ylabel(sprintf('Slope Efficiency (%%) and\n Wavelength Discrimination (dB)'));
        yR = get(ax(2),'YLim');
        nTick1 = numel(get(ax(1),'YTick'));
        set(ax(2),'YTick',linspace(0,yR(2),nTick1));
        h4 = plot(ax(2),100*Lcalc,wlDisc); 
        h3 = plot(ax(1),100*Lcalc,1000*Psop);        
        set(h1,'linewidth',2,'color','r','marker','o');
        set(h2,'linewidth',2,'color',[0 0.5 0],'marker','o');
        set(h3,'linewidth',2,'color','b','marker','o');
        set(h4,'linewidth',2,'color','k','marker','o');
        hL = legend([h2,h1,h3,h4],{'Efficiency','Threshold P_P','Operating P_{out}','Wavelength Discrim.'},'Location','Best'); set(hL,'Color','w');
        hT = title(sprintf('Laser Performance vs Fiber Length: R_{out}=%d%%, Losses=%.1fdB',round(100*s.R2),s.eps1db+s.eps2db));
        data = {eta,Ppth,Grt,Lcalc,s};  
    case 'lo' % plot versus cavity losses
        %   Plot efficiency and threshold versus active fiber length
        [ax h1 h2] = plotyy(Epscalc,1000*Ppth,Epscalc,100*eta);
        axes(ax(1)); set(gca,'fontsize',16,'linewidth',2,'ytickmode','auto'); hold on; grid on;
        xlabel('Cavity Losses (dB)');
        ylabel('Power (mW)');
        ylim([0 1000*Pop]);
        axes(ax(2)); hold on; set(gca,'fontsize',16,'linewidth',2); grid on;
        ylabel(sprintf('Slope Efficiency (%%) and\n Wavelength Discrimination (dB)'));
        yR = get(ax(2),'YLim');
        nTick1 = numel(get(ax(1),'YTick'));
        set(ax(2),'YTick',linspace(0,yR(2),nTick1));
        h4 = plot(ax(2),Epscalc,wlDisc); 
        h3 = plot(ax(1),Epscalc,1000*Psop);        
        set(h1,'linewidth',2,'color','r','marker','o');
        set(h2,'linewidth',2,'color',[0 0.5 0],'marker','o');
        set(h3,'linewidth',2,'color','b','marker','o');
        set(h4,'linewidth',2,'color','k','marker','o');
        hL = legend([h2,h1,h3,h4],{'Efficiency','Threshold P_P','Operating P_{out}','Wavelength Discrim.'},'Location','Best'); set(hL,'Color','w');
        hT = title(sprintf('Laser Performance vs Cavity Losses: R_{out}=%d%%, Length=%dcm',round(100*s.R2),round(100*s.L)));
        data = {eta,Ppth,Grt,Epscalc,s};        
    case {'bo','3d','po'} % 3D plot of both versus power
        surf(100*Lcalc',100*Rcalc,1000*Psop);
        set(gca,'fontsize',16,'linewidth',2,'ytickmode','auto'); 
        ylabel('Reflectance (%)');
        xlabel('Fiber Length (cm)');
        zlabel('Output Power (mW)');
        hT = title(sprintf('Laser Output Power vs Reflectance\n and Fiber Length: Losses=%.1fdB',s.eps1db+s.eps2db));
        data = {Psop,Lcalc,Rcalc,s};
end

end

%%% END OF MAIN FUNCTION %%%

function wlDisc = CheckWavelengthDiscrimination(Ppth,s)

% Call small-signal gain function to find wavelength-dependent gain at
% threshold
s.plotFlag = 0;
s.Pp = Ppth;
s.mode = 'wavelength';
data = SmallSignalGain(s);
G = data{1};
wl = data{2};

% Find gain at operating wavelength; compare to gain at peak
[Gmax iimax] = max(G);
wlmax = wl(iimax);
iiwlS = find(wl>=s.lamS,1,'first');
Gs = G(iiwlS);
wlDisc = Gmax-Gs;

end