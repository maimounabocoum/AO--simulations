% Semi-GUI for fiber laser and amplifier analysis
function data = FiberAnalysis

% Select which device to analyze
prompt = {sprintf(['Which fiber device is being studied?\n' ...
                      '(EDFA, EDFL, YDFA, YDFL)']); 
          'Is the device cladding-pumped?'};     
def = {'EDFA','No'};      
device = inputdlg(prompt,'Fiber Lasers and Amplifiers',1,def);
switch upper(device{1})
    case 'EDFA'
        rareEarth = 'erbium';
        islaser = 0;
    case 'EDFL'
        rareEarth = 'erbium';
        islaser = 1;
    case 'YDFA'        
        rareEarth = 'ytterbium';
        islaser = 0;
    case 'YDFL'
        rareEarth = 'ytterbium';
        islaser = 1;
end
switch lower(device{2}(1))
    case 'y'
        pump = 'clad';
    otherwise
        pump = 'core';
end

% Load default device parameters, ask user to review
s = GetDefaultSettings(rareEarth,pump);
%               % prompt    % default value
if strcmpi(s.re(1),'e')
    absStr = 'at 1530 nm';
else
    absStr = 'at 975 nm';
end
settingsBox = { ...
    'Pump wavelength (nm)',     num2str(s.lamP); ...
    'Signal wavelength (nm)',   num2str(s.lamS); ...
    'Pump power (mW)',          num2str(s.Pp); ...
    'Signal power (mW)',        num2str(s.Ps*(~islaser)); ...
    'Doped fiber length (m)',   num2str(s.L); ...
    'Doped fiber core diameter (um)', num2str(s.dCore); ...
    'Core E-field overlap with dopant (%)', num2str(100*s.Gamma); ...
    'Doping concentration (x 10^24 1/m3)', num2str(1e-24*s.Nt); ...
    ['(OR) Fiber absorption (dB/m) ' absStr], num2str(s.alph); ...
    'Pump direction',           num2str(s.direction); ...
    'Pumping approach',         num2str(s.pump); 
    };
settings = inputdlg(settingsBox(:,1), ...
    'System Parameters', 1, settingsBox(:,2));
s.lamP = 1e-9 * str2double(settings{1});
s.lamS = 1e-9 * str2double(settings{2});
s.Pp = 1e-3 * str2double(settings{3});
s.Ps = 1e-3 * str2double(settings{4});
s.L = str2double(settings{5}); 
s.dCore = 1e-6 * str2double(settings{6}); 
s.Gamma = 1e-2 * str2double(settings{7});
s.dField = s.dCore / s.Gamma;
if ~isempty(settings{8}), 
    s.Nt = 1e24 * str2double(settings{8});
else
    if str2cmpi(s.re(1),'e')
        s.Nt = ConvAbsDB2N(str2double(settings{9}),GetErSpectrum(1530));
    else
        s.Nt = ConvAbsDB2N(str2double(settings{9}),GetYbSpectrum(975));
    end
end
s.direction = lower(settings{10}(1));
s.pump = lower(settings{11});

% Ask user which functions they want to run
[fcnList iiList] = GetDefaultFunctions(islaser);
[selections,ok] = listdlg('PromptString','Select functions to run:', ...
    'SelectionMode','multiple',...
    'ListString',fcnList, ...
    'ListSize',[250 300]);  
iiSelections = iiList(selections);

if ok 
    % Run functions
    for ii = 1:length(selections)
        [data{ii} hFig(ii)] = CallDesiredFunctions(iiSelections(ii),s);
    end
end

end
%%% *** END OF MAIN FUNCTION *** %%%

%%% *** BEGIN SUPPORTING FUNCTIONS *** %%%
function s = GetDefaultSettings(re,pump)

% common, system options
s.direction = 'forward';
s.pump = pump;
s.re = lower(re(1));

switch [s.re pump]
    case 'ecore'
        s.L = 10;               % [m] length of active fiber
        s.Ps = 0.03;            % [mW] signal power
        s.Pp = 100;             % [mW] pump power
        s.lamS = 1550;          % [nm] signal wavelength
        s.lamP = 980;           % [nm] signal wavelength
        s.dCore = 5.5;            % [um] doped core diameter
        s.Gamma = 0.75;         % [1] core/field overlap ratio
        s.dField = s.dCore/s.Gamma;   % [um] field diameter
        s.alph = 6.5;          % [dB/m] signal absorption rate
        s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(1530)); % [1/m3] doping concentration
    case 'eclad'
        s.L = 50;               % [m] length of active fiber
        s.Ps = 10;            % [mW] signal power
        s.Pp = 10000;             % [mW] pump power
        s.lamS = 1550;          % [nm] signal wavelength
        s.lamP = 980;           % [nm] signal wavelength
        s.dCore = 9;            % [um] doped core diameter
        s.Gamma = 0.75;         % [1] core/field overlap ratio
        s.dField = s.dCore/s.Gamma;   % [um] field diameter
        s.alph = 100;          % [dB/m] signal absorption rate
        s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(980)); % [1/m3] doping concentration        
    case 'ycore'
        s.L = 1;       % [m] length of active fiber
        s.Ps = 5;      % [mW] signal power
        s.Pp = 225;       % [mW] pump power
        s.lamS = 1064;    % [nm] signal wavelength
        s.lamP = 975;    % [nm] signal wavelength
        s.dCore = 6;    % [um] doped core diameter
        s.Gamma = 0.75; % [1] core/field overlap ratio
        s.dField = s.dCore/s.Gamma;   % [um] field diameter
        s.alph = 350;  % [dB/m] signal absorption rate
        s.Nt = ConvAbsDB2N(s.alph,GetYbSpectrum(975)); % [1/m3] doping concentration   
    case 'yclad'
        s.L = 5;       % [m] length of active fiber
        s.Ps = 200;      % [mW] signal power
        s.Pp = 15000;       % [mW] pump power
        s.lamS = 1064;    % [nm] signal wavelength
        s.lamP = 915;    % [nm] signal wavelength
        s.dCore = 10;    % [um] doped core diameter
        s.Gamma = 0.75; % [1] core/field overlap ratio
        s.dField = s.dCore/s.Gamma;   % [um] field diameter
        s.alph = 110;  % [dB/m] signal absorption rate
        s.Nt = ConvAbsDB2N(s.alph,GetYbSpectrum(915)); % [1/m3] doping concentration         
end

end

function [fcnList iiList] = GetDefaultFunctions(islaser)

fcnList = GetFunctionList();

iiList = [];
if islaser
    for ii = 1:size(fcnList,1)
        % laser functions
        if fcnList{ii,2}==1, iiList = [iiList ii]; end
    end
else
    for ii = 1:size(fcnList,1)
        % amplifier functions
        if fcnList{ii,2}==0, iiList = [iiList ii]; end
    end
end
fcnList = fcnList(iiList);

end

function [data hFig] = CallDesiredFunctions(sel,set)

switch sel
    
    %%% BEGIN LASER FUNCTIONS %%% 
    
    case 1 % Gain vs. wavelength (below threshold)
        % Call analytical solution
        set.mode = 'wavelength';
        set.plotFlag = 1;
        set.multPflag = 1;
        set.multLflag = 1;
        [data hFig] = SmallSignalGain(set);
    case 2 % Laser performance vs. reflectance
        % Request laser architecture, cavity losses from user
        details = inputdlg({sprintf(['Plotting Performance vs. Reflectance:\n' ...
            'Please provide laser configuration details:\n' ... 
            'Resonator type (linear, ring, feedback ring):']); ...
            'Cavity losses between fiber and coupler (dB)'; ...
            'Cavity losses between coupler and fiber (dB)'; ...
            'Active fiber length (cm)'}, ...
            'Laser Configuration', 1, ...
            {'Ring';'3';'3';num2str(100*set.L)});
        set.arch = lower(details{1});
        set.eps2db = str2double(details{2});
        set.eps1db = str2double(details{3});
        set.L = str2double(details{4})/100;
        % Call analytical solution
        set.mode = 'reflectance';
        set.plotFlag = 1;
        [data hFig] = LaserPerformanceManager(set);
    case 3 % Laser performance vs. fiber length
        % Request laser architecture, cavity losses from user
        details = inputdlg({sprintf(['Plotting Performance vs. Fiber Length:\n' ...
            'Please provide laser configuration details:\n' ... 
            'Resonator type (linear, ring, feedback ring):']); ...
            'Cavity losses between fiber and coupler (dB)'; ...
            'Cavity losses between coupler and fiber (dB)'; ...
            'Output Coupler Reflectance (%)'}, ...
            'Laser Configuration', 1, ...
            {'Ring';'3';'3';'50'});
        set.arch = lower(details{1});
        set.eps2db = str2double(details{2});
        set.eps1db = str2double(details{3});
        set.R = str2double(details{4})/100;     
        % Call analytical solution
        set.mode = 'length';
        set.plotFlag = 1;
        [data hFig] = LaserPerformanceManager(set);       
    case 4 % Laser performance vs. cavity losses
        % Request laser architecture, reflectance, fiber length from user
        details = inputdlg({sprintf(['Plotting Performance vs. Cavity Losses:\n' ...
            'Please provide laser configuration details:\n' ... 
            'Resonator type (linear, ring, feedback ring):']); ...
            'Active fiber length (cm)'; ...
            'Output coupler reflectance (%)'}, ...
            'Laser Configuration', 1, ...
            {'Ring';num2str(set.L*100);'50'});
        set.arch = lower(details{1});
        set.L = str2double(details{2})/100; 
        set.R = str2double(details{3})/100;     
        % Call analytical solution
        set.mode = 'losses';
        set.plotFlag = 1;
        [data hFig] = LaserPerformanceManager(set);              
    case 5 % Laser output power vs. reflectance, length
        % Request laser architecture, reflectance, fiber length from user
        details = inputdlg({sprintf(['Plotting Power vs. Reflectance, Length:\n' ...
            'Please provide laser configuration details:\n' ... 
            'Resonator type (linear, ring, feedback ring):']); ...
            'Active fiber length (cm)'; ...
            'Output coupler reflectance (%)'; ...
            'Cavity losses between fiber and coupler (dB)'; ...
            'Cavity losses between coupler and fiber (dB)'}, ...
            'Laser Configuration', 1, ...
            {'Ring';num2str(set.L*100);'50';'3';'3'});
        set.arch = lower(details{1});
        set.L = str2double(details{2})/100; 
        set.R = str2double(details{3})/100;     
        set.eps2db = str2double(details{4});
        set.eps1db = str2double(details{5});        
        % Call analytical solution
        set.mode = 'power';
        set.plotFlag = 1;
        [data hFig] = LaserPerformanceManager(set);          
    case 6 % Pump and signal powers in fiber
        % Get more data from user        
        % Request laser architecture, cavity losses from user
        details = inputdlg({sprintf(['Plotting Performance vs. Reflectance:\n' ...
            'Please provide laser configuration details:\n' ... 
            'Resonator type (linear, ring, feedback ring):']); ...
            'Cavity losses between fiber and coupler (dB)'; ...
            'Cavity losses between coupler and fiber (dB)'; ...
            'Active fiber length (cm)'; ...
            'Output Coupler Reflectance (%%)'}, ...            
            'Laser Configuration', 1, ...
            {'Ring';'3';'3';num2str(100*set.L);'50'});
        set.arch = lower(details{1});
        set.eps2db = str2double(details{2});
        set.eps1db = str2double(details{3});
        set.L = str2double(details{4})/100;  
        set.R = str2double(details{5})/100;
        % Call numerical solution
        set.mode = 'power';
        set.plotFlag = 1;
        [data hFig] = LaserNumerical(set);          

    %%% BEGIN AMPLIFIER FUNCTIONS %%% 
        
    case 7 % Small signal gain vs wavelength
        set.mode = 'wavelength';
        set.plotFlag = 1;
        set.multPflag = 1;
        set.multLflag = 1;
        [data hFig] = SmallSignalGain(set);
    case 8 % Small signal gain vs pump power
        set.mode = 'power';
        set.plotFlag = 1;
        set.multPflag = 1;
        set.multLflag = 1;
        [data hFig] = SmallSignalGain(set);        
    case 9 % Pump and signal powers in fiber
        % Call numerical solution
        set.mode = 'power';
        set.plotFlag = 1;
        [data hFig] = AmplifierPerformance(set);          
    case 10 % Pump and signal gain in fiber
        % Call numerical solution
        set.mode = 'gain';
        set.plotFlag = 1;
        [data hFig] = AmplifierPerformance(set);           
    case 11 % Output power vs wavelength, length
        % Ask user if he wants to see even gain
        details = inputdlg({sprintf(['Plotting Power vs. Wavelength, Length:\n' ...
            'Please select a mode:\n' ... 
            '1 -- ASE in presence of strong signal or\n 2 -- ASE gain from flat input spectrum:']); ...
            'ASE input signal level (uW)'; ...
            'Active fiber length (cm)'}, ...
            'ASE Preferences', 1, ...
            {'1 -- ASE with signal';'0';num2str(set.L*100)});
        ASEmode = details{1}(1);
        set.L = str2double(details{3})/100; 
        switch lower(ASEmode)
            case {'1','a'}
                set.Pase = 0;
            case {'2','s'}
                set.Pase = str2double(details{2})*1e-6;
                set.Ps = set.Pase;
        end
        % Call numerical solution
        set.mode = 'ase';
        set.plotFlag = 1;
        [data hFig] = AmplifierPerformance(set);    
    case 12 % Output power vs fiber length (don't see profile of fiber)
        % Call numerical solution
        set.mode = 'length';
        [data hFig] = AmplifierPerformanceManager(set);
    case 13 % Output power vs pump power (don't see profile of fiber)
        % Call numerical solution
        set.mode = 'pump';
        [data hFig] = AmplifierPerformanceManager(set);        
    case 14 % Output power vs signal power (don't see profile of fiber)
        % Call numerical solution
        set.mode = 'signal';
        [data hFig] = AmplifierPerformanceManager(set);        
    case 15 % Output power vs length, signal power (don't see profile of fiber)
        % Call numerical solution
        set.mode = '3D';
        [data hFig] = AmplifierPerformanceManager(set);        
end


end

function functionList = GetFunctionList()

% list is: [displayStr   laserFlag   functionName]
%   displayStr will be seen by user
%   laserFlag indicates if it should be called when analyzing
%   a laser (1) or an amplifier (0) 
displayStr = {'Gain vs. Wavelength (Below Threshold)'       1; ...
                'Laser Performance vs. Reflectance'         1; ...
                'Laser Performance vs. Fiber Length'        1; ...
                'Laser Performance vs. Cavity Loss'         1; ...
                'Laser Output Power vs. Reflectance, Length' 1; ...
                'Pump and Signal Powers in Fiber'           1; ...
                'Small Signal Gain vs. Wavelength'          0; ...
                'Small Signal Gain vs. Pump Power'          0; ...
                'Pump and Signal Powers in Fiber'           0; ...
                'Pump and Signal Gain in Fiber'             0; ...
                'Output Power vs. Wavelength, Length'       0; ...
                'Output Power vs. Fiber Length'             0; ...                
                'Output Power vs. Pump Power'               0; ...
                'Output Power vs. Signal Power'             0; ...
                'Output Power vs. Signal, Length'           0 ...
                };
%   functionName is the function to call. each function is expected to
%   handle the call using a standardized settings structure, returning 
%   data and a figure handle
functionNames = {'SmallSignalGain'; ...
                'LaserPerformanceManager'; ...
                'LaserPerformanceManager'; ...
                'LaserPerformanceManager'; ...
                'LaserPerformanceManager'; ...
                'AmplifierPerformance'; ...
                'SmallSignalGain'; ... 
                'SmallSignalGain'; ... 
                'AmplifierPerformance'; ...
                'AmplifierPerformance'; ...
                'AmplifierPerformance'; ...
                'AmplifierPerformanceManager'; ...
                'AmplifierPerformanceManager'; ...
                'AmplifierPerformanceManager'; ...
                'AmplifierPerformanceManager' ...
                };
functionList = [displayStr  functionNames];
            
end