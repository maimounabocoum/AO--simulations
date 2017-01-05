%  Make the scatteres for a simulation and store
%  it in a file for later simulation use

%   Joergen Arendt Jensen, April 2, 1998

[phantom_positions, phantom_amplitudes] = tissue_pht(100000);
save pht_data.mat phantom_positions phantom_amplitudes
