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
    Z1 = hypot(U1, V1);
end

function [X2, Y2, U2, V2, Z2] = getLiquidData(S, frame)
    X2 = S.all_transposed_x_position_matrix_liquid(:, :, frame) * 1e3;
    Y2 = S.all_transposed_y_position_matrix_liquid(:, :, frame) * 1e3;
    U2 = S.all_u_matrix_liquid(:, :, frame);
    V2 = S.all_v_matrix_liquid(:, :, frame);
    Z2 = hypot(U2, V2);
end

function [X3, Y3] = getFilmData(S, frame)
    X3 = S.film_x_position(1, 15:1090) * 1e3;
    Y3 = S.smoothed_film_height_matrix_out(15:1090, frame) * 1e3;
end