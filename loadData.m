function S = loadData(filename)    
    filename = 'L8_G3.mat';
    tmp = load(filename);
    S = tmp.S2P_PIV_full_mat_vars;
end
