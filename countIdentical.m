function [labels, counts] = countIdentical(data)
    % Return sorted unique elements of the input vector (labels) and the
    % number of occurences of each element (counts).
    dataNoNaN = sort(data(~isnan(data)));
    NaNCount = length(data) - length(dataNoNaN);
    
    labels = unique(dataNoNaN);
    counts = diff(find([1, diff(dataNoNaN), 1]));
    
    if NaNCount
        labels = [labels, NaN];
        counts = [counts, NaNCount];
    end
end