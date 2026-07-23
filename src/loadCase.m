function [S, meta] = loadCase(case_name)
% loads one case from data/Cases by its name, e.g. 'L8_G6'
% gives back the data struct plus the case metadata from caseConfig
% loadData is still there as the low level loader if you have a full path

meta = caseConfig(case_name);

% data lives in data/Cases relative to the project root, not the caller
src_dir   = fileparts(mfilename('fullpath'));
proj_root = fileparts(src_dir);
mat_path  = fullfile(proj_root, 'data', 'Cases', [case_name '.mat']);

if ~exist(mat_path, 'file')
  error('loadCase:missingFile', ...
    'no data file for %s, expected it at %s', case_name, mat_path);
end

[S, ~] = loadData(mat_path);

% keep track of what the file actually contains, the expected frame
% count is not always what was captured
meta.actual_frames = size(S.all_u_matrix_liquid, 3);

if meta.actual_frames ~= meta.expected_frames
  fprintf('note: %s has %d frames, expected %d\n', ...
    case_name, meta.actual_frames, meta.expected_frames);
end

end
