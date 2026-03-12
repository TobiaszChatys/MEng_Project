function [Z1_masked, Z2_masked] = createMasks(X1, Y1, X2, Y2, Z1, Z2, X3, Y3)
    Y3_air    = interp1(X3, Y3, X1(1,:), 'linear', 'extrap');
    Y3_liquid = interp1(X3, Y3, X2(1,:), 'linear', 'extrap');

    Y3_matrix_air    = repmat(Y3_air,    size(Y1,1), 1);
    Y3_matrix_liquid = repmat(Y3_liquid, size(Y2,1), 1);

    Z1_masked = Z1;
    Z1_masked(Y1 < Y3_matrix_air) = NaN;   % keep only above interface

    Z2_masked = Z2;
    Z2_masked(Y2 >= Y3_matrix_liquid) = NaN; % keep only below interface
end