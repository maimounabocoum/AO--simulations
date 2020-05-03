% Run analysis on various fiber laser configurations
function [data hFig] = AmplifierPerformanceManager(s)

% allow comparison of amplifier performances, varying by
%   1. pump power
%   2. signal power
%   3. fiber length

%% INITIALIZATION %%

% Resolve Inputs
if nargin<1, s = struct; end % create dummy structure if not present
if ~isfield(s,'mode'), s.mode = 'length'; end % plot vs length, pump, or signal -- or 3D plot of power vs length, signal
managerMode = lower(s.mode(1));
if ~isfield(s,'plotFlag'), s.plotFlag = 0; end
s.mode = 'power'; 

% Define variables against which to run calculations
switch managerMode
    case 'l' % plot versus length
        if ~isfield(s,'L'), s.L = 1; end
        if ~isfield(s,'autoLflag'), s.autoLflag = true; end
        if s.autoLflag
            Lcalc = s.L*(10:10:100)/100; % usually generate range of lengths automatically
        else
            Lcalc = s.L; % allow input of discrete lengths if desired
        end
    case 'p' % plot versus pump power
        if ~isfield(s,'Pp'), s.Pp = 0.200; end
        if ~isfield(s,'autoPflag'), s.autoPflag = true; end
        if s.autoPflag
            Pcalc = s.Pp*(10:10:100)/100;
        else
            Pcalc = s.Pp;
        end
    case 's' % plot versus signal power
        if ~isfield(s,'Ps'), s.Ps = 0.010; end
        if ~isfield(s,'autoSflag'), s.autoSflag = true; end
        if s.autoSflag
            Scalc = s.Ps*(10:10:100)/100;
        else
            Scalc = s.Ps;
        end
    case {'3'} % 3D plot of output power versus input signal, fiber length
        % set up signals
        if ~isfield(s,'Ps'), s.Ps = 0.010; end
        if ~isfield(s,'autoSflag'), s.autoSflag = true; end
        if s.autoSflag
            Scalc = s.Ps*(25:25:100)/100;
        else
            Scalc = s.Ps;
        end
        % set up lengths
        if ~isfield(s,'L'), s.L = 1; end
        if ~isfield(s,'autoLflag'), s.autoLflag = true; end
        if s.autoLflag
            Lcalc = s.L*(25:25:100)/100; % usually generate range of lengths automatically
        else
            Lcalc = s.L; % allow input of discrete lengths if desired
        end        
end

%% CALCULATION %%

switch managerMode
    case 'l' % plot versus length
        Pout = zeros(1,length(Lcalc));
        NF = zeros(1,length(Lcalc));
        for ll = 1:length(Lcalc)
            s.L = Lcalc(ll);
            [Pout(ll),~,sUsed] = AmplifierPerformance(s);
        end
    case 'p' % plot versus pump power 
        Pout = zeros(1,length(Pcalc));
        NF = zeros(1,length(Pcalc));
        for pp = 1:length(Pcalc)
            s.Pp = Pcalc(pp);
            [Pout(pp),~,sUsed] = AmplifierPerformance(s);
        end
    case 's' % plot versus signal power
        Pout = zeros(1,length(Scalc));
        NF = zeros(1,length(Scalc));
        for ss = 1:length(Scalc)
            s.Ps = Scalc(ss);
            [Pout(ss),~,sUsed] = AmplifierPerformance(s);
        end       
    case '3' % 3D plot of output power versus length, signal
        Pout = zeros(length(Lcalc),length(Scalc));
        NF = zeros(length(Lcalc),length(Scalc));
        for ll = 1:length(Lcalc)
            s.L = Lcalc(ll);
            for ss = 1:length(Scalc)
                s.Ps = Scalc(ss);
                [Pout(ll,ss),~,sUsed] = AmplifierPerformance(s);
            end
        end
end

%% PLOTS %%

hFig = figure('windowstyle','docked');

switch managerMode
    case 'l' % plot versus length        
        %   Plot output power and noise figure versus length
        [ax h1 h2] = plotyy(100*Lcalc,1000*Pout,100*Lcalc,NF);
        axes(ax(1)); set(gca,'fontsize',16,'linewidth',2,'ytickmode','auto'); hold on; grid on;
        xlabel('Active Fiber Length (cm)');
        ylabel('Output Power (mW)');
        ylim([0 1100*max(Pout)]);
        axes(ax(2)); hold on; set(gca,'fontsize',16,'linewidth',2); grid on;
        ylabel('Noise Figure (dB)');
        yR = get(ax(2),'YLim');
        nTick1 = numel(get(ax(1),'YTick'));
        set(ax(2),'YTick',linspace(0,yR(2),nTick1));    
        set(h1,'linewidth',2,'color','b','marker','o');
        set(h2,'linewidth',2,'color','g','marker','o');
        hL = legend([h1,h2],{'Output Power','Noise Figure'},'Location','Best'); set(hL,'Color','w');
        hT = title(sprintf('Amplifier Performance vs Fiber Length: P_P=%dmW, P_S=%dmW',round(1000*sUsed.Pp),round(1000*sUsed.Ps)));
        data = {Pout,Lcalc,sUsed};
    case 'p' % plot versus pump power
        %   Plot output power and noise figure versus pump power
        [ax h1 h2] = plotyy(1000*Pcalc,1000*Pout,1000*Pcalc,NF);
        axes(ax(1)); set(gca,'fontsize',16,'linewidth',2,'ytickmode','auto'); hold on; grid on;
        xlabel('Pump Power (mW)');
        ylabel('Output Power (mW)');
        ylim([0 1100*max(Pout)]);
        axes(ax(2)); hold on; set(gca,'fontsize',16,'linewidth',2); grid on;
        ylabel('Noise Figure (dB)');
        yR = get(ax(2),'YLim');
        nTick1 = numel(get(ax(1),'YTick'));
        set(ax(2),'YTick',linspace(0,yR(2),nTick1));    
        set(h1,'linewidth',2,'color','b','marker','o');
        set(h2,'linewidth',2,'color','g','marker','o');
        hL = legend([h1,h2],{'Output Power','Noise Figure'},'Location','Best'); set(hL,'Color','w');
        hT = title(sprintf('Amplifier Performance vs Pump Power: P_S=%dmW, L=%dcm',round(1000*sUsed.Ps),round(100*sUsed.L)));
        data = {Pout,Pcalc,sUsed};  
    case 's' % plot versus signal power
        %   Plot output power and noise figure versus signal power
        [ax h1 h2] = plotyy(1000*Scalc,1000*Pout,1000*Scalc,NF);
        axes(ax(1)); set(gca,'fontsize',16,'linewidth',2,'ytickmode','auto'); hold on; grid on;
        xlabel('Signal Power (mW)');
        ylabel('Output Power (mW)');
        ylim([0 1100*max(Pout)]);
        axes(ax(2)); hold on; set(gca,'fontsize',16,'linewidth',2); grid on;
        ylabel('Noise Figure (dB)');
        yR = get(ax(2),'YLim');
        nTick1 = numel(get(ax(1),'YTick'));
        set(ax(2),'YTick',linspace(0,yR(2),nTick1));    
        set(h1,'linewidth',2,'color','b','marker','o');
        set(h2,'linewidth',2,'color','g','marker','o');
        hL = legend([h1,h2],{'Output Power','Noise Figure'},'Location','Best'); set(hL,'Color','w');
        hT = title(sprintf('Amplifier Performance vs Signal Power: P_P=%dmW, L=%dcm',round(1000*sUsed.Pp),round(100*sUsed.L)));
        data = {Pout,Scalc,sUsed};        
    case '3' % 3D plot of output power versus signal, length
        surf(1000*Scalc',100*Lcalc,1000*Pout);
        set(gca,'fontsize',16,'linewidth',2,'ytickmode','auto'); 
        xlabel('Signal Power (mw)');
        ylabel('Fiber Length (cm)');
        zlabel('Output Power (mW)');
        hT = title(sprintf('Amplifier Output Power vs Fiber Length\n and Input Signal Power: Pump Power=%dmW',round(1000*sUsed.Pp)));
        data = {Pout,Lcalc,Scalc,sUsed};
end

end

%%% END OF MAIN FUNCTION %%%