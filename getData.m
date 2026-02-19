function [X1, Y1, U1, V1, Z1, X2, Y2, U2, V2, Z2, X3, Y3] = getData(S, frame)
    
    [X1, Y1, U1, V1, Z1] = getAirData(S, frame);
    [X2, Y2, U2, V2, Z2] = getLiquidData(S, frame);
    [X3, Y3] = getFilmData(S, frame);
end 

function [X1, Y1, U1, V1, Z1] = getAirData(S, frame)
    X1 = S.all_transposed_x_position_matrix_air(:, :, frame) * 1e3;
    Y1 = S.all_transposed_y_position_matrix_air(:, :, frame) * 1e3;
    U1 = S.all_u_matrix_air(:, :, frame);
    V1 = S.all_v_matrix_air(:, :, frame);

    % Set zero velocities to NaN
    zero_mask = (U1 == 0) & (V1 == 0);
    U1(zero_mask) = NaN;
    V1(zero_mask) = NaN;

    Z1 = hypot(U1, V1);

    x_air = X1(1, :); % 1 x 50 (x positions)
    film_height_at_liquid_x = interp1(X3, Y3, x_air, 'linear', 'extrap'); % 1 x 50 

    liquid_above_interface = Y1 > film_height_at_liquid_x; % 99 x 50 logical matrix

    U1(liquid_above_interface) = NaN;
    V1(liquid_above_interface) = NaN;
end

function [X2, Y2, U2, V2, Z2] = getLiquidData(S, frame)
    X2 = S.all_transposed_x_position_matrix_liquid(:, :, frame) * 1e3;
    Y2 = S.all_transposed_y_position_matrix_liquid(:, :, frame) * 1e3;
    U2 = S.all_u_matrix_liquid(:, :, frame);
    V2 = S.all_v_matrix_liquid(:, :, frame);

    % Set zero velocities to NaN
    zero_mask = (U2 == 0) & (V2 == 0);
    U2(zero_mask) = NaN;
    V2(zero_mask) = NaN;

    Z2 = hypot(U2, V2);

    % For any liquid vertical positon (Y2) that is less than the film height at the same x positon,
    % set the liquid velocity to NaN

    x_liquid = X2(1, :); % 1 x 50 (x positions)
    film_height_at_liquid_x = interp1(X3, Y3, x_liquid, 'linear', 'extrap'); % 1 x 50 

    liquid_above_interface = Y2 > film_height_at_liquid_x; % 99 x 50 logical matrix

    U2(liquid_above_interface) = NaN;
    V2(liquid_above_interface) = NaN;

end

function [X3, Y3] = getFilmData(S, frame)
    X3 = S.film_x_position(1, :) * 1e3;
    Y3 = S.smoothed_film_height_matrix_out(:, frame) * 1e3;

    remove_negative_indices = Y3 >= 0;
    X3 = X3(remove_negative_indices);
    Y3 = Y3(remove_negative_indices);
end