function xlim_value = findXlimFromXticks(max_x, xticks_array)
    % round up the maximum x value
    max_x_rounded = ceil(max_x);

    % find the smallest tick value that is >= max_x_rounded
    valid_xticks = xticks_array(xticks_array >= max_x_rounded);

    if ~isempty(valid_xticks)
        xlim_value = min(valid_xticks);
    else
        xlim_value = max(xticks_array);
    end
end
