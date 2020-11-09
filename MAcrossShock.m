function m2 = MAcrossShock(M1)
    gamma = 1.4;
    m2 = sqrt(((gamma - 1) .* M1.^2 + 2) ./ (2*gamma .* M1.^2 - (gamma - 1)));
end
