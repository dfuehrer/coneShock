function [beta, betaStrong] = dbmBeta(delta, M)
    gamma = 1.4;

    betaEq = [1, -((M.^2 + 2)./M.^2 + gamma.*sind(delta).^2), ((2*M.^2+1)./M.^4 + ((gamma+1).^2/4 + (gamma-1)./M.^2).*sind(delta).^2), -cosd(delta).^2./M.^4];
    bs = roots(betaEq);

    B = asind(sqrt(bs(bs > min(bs))));
    beta = min(B); betaStrong = max(B);
end
