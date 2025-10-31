function [S, filename] = loadData(filename)    
    tmp = load(filename);
    S = tmp.S2P_PIV_full_mat_vars;
end
