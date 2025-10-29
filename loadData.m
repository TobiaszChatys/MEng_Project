function S = loadData(filename)
    if nargin < 1
        filename = 'L8_G3.mat';
    end
    tmp = load(filename);
    S = tmp.S2P_PIV_full_mat_vars;
end
