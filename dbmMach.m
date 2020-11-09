%dbmMach.m
function M = dbmMach(delta, Beta, gamma)
    if nargin < 3
        gamma = 1.4;
    end
    d2r = pi / 180;
    M =  sqrt(2 * (cotd(Beta) + tand(delta)) / (sind(2*Beta) - tand(delta) * (gamma + cosd(2*Beta))));
end
