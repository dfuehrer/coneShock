function delta = dbmDelta(Beta, M, gamma)
    if nargin == 2
        gamma = 1.4;
    end
    %delta = atan2d(2*cotd(Beta) .* (M.^2 .* sind(Beta).^2 - 1), M.^2 .* (gamma + 1 - 2*sind(Beta).^2) + 2);
    delta = atan2d(M.^2 .* sind(2*Beta) - 2*cotd(Beta), M.^2 .* (gamma + cosd(2*Beta)) + 2);
end
