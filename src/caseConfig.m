function config = caseConfig(case_name)
% returns the metadata for one case, or the full table if no name is given
% sampling rates and frame counts are from the experimental setup,
% mode selections are the ones settled on for the final report
%
% usage:
%   all_cases = caseConfig();
%   case_info = caseConfig('L8_G6');

% gas velocity (m/s), sampling rate (hz), expected frames, regime,
% liquid modes, gas modes. NaN modes means no spatial reconstruction
% was done for that case (too chaotic to pick a cutoff objectively)
case_table = {
  'L8_G1',   1,  500, 2500, '2D',          30,  100;
  'L8_G2',   2,  500, 2500, '2D',          30,  100;
  'L8_G3',   3,  500, 2500, '2D',          30,  400;   % transition case, air needs a lot more modes
  'L8_G4',   4, 1000, 2500, '3D',          50,  250;
  'L8_G5',   5, 1000, 2500, '3D',          50,  250;
  'L8_G6',   6, 1000, 2500, '3D',          50,  250;
  'L8_G7',   7, 1000, 2500, '3D',          50,  250;
  'L8_G8',   8, 1000, 2500, '3D',          75,  250;   % upper 3D, needs more modes
  'L8_G9',   9, 1000, 2500, 'disturbance', 75,  250;
  'L8_G10', 10, 1000, 2500, 'disturbance', NaN, NaN;
  'L8_G12', 12, 1500, 3750, 'disturbance', NaN, NaN;
  };

% build the struct array
number_of_cases = size(case_table, 1);
config = struct();

for row = 1:number_of_cases
  config(row).name            = case_table{row, 1};
  config(row).gas_velocity    = case_table{row, 2};
  config(row).fs              = case_table{row, 3};
  config(row).expected_frames = case_table{row, 4};
  config(row).regime          = case_table{row, 5};
  config(row).n_liquid_modes  = case_table{row, 6};
  config(row).n_gas_modes     = case_table{row, 7};
end

% no argument means give back the whole table
if nargin == 0
  return
end

% otherwise find the requested case
match = strcmp({config.name}, case_name);

if ~any(match)
  error('caseConfig:unknownCase', ...
    'unknown case "%s", expected something like L8_G6', case_name);
end

config = config(match);

end
