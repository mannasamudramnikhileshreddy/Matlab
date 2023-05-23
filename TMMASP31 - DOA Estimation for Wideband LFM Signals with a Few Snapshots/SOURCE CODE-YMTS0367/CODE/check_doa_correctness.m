function correct = check_doa_correctness(actual, est, tolerance)
if nargin == 2
    tolerance = inf;
end
n = length(actual);
if n ~= length(est)
    correct = false;
    return;
end
if n == 1
    correct = abs(actual(1) - est(1)) < tolerance;
elseif n > 1
    spaces = diff(actual);
    max_deviation = min(spaces/2, tolerance);
    correct = true;
    % first
    correct = correct && ...
        (est(1) >= max(-pi/2, actual(1)-tolerance) && ...
        est(1) < actual(1) + max_deviation(1)); 
    % last
    correct = correct && ...
        (est(end) <= min(pi/2, actual(end)+tolerance) && ...
        est(end) > actual(end) - max_deviation(end));
    % middle
    for ii = 2:(n-1)
        correct = correct && ...
        (est(ii) > actual(ii) - max_deviation(ii-1) && ... 
        est(ii) < actual(ii) + max_deviation(ii));
    end
end

end

