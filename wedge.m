
M = [1.5, 2, 5];    % list of mach numbers to run
ds = 0:.2:55;          % list of cone angles to calculate
Betas = zeros(length(M), length(ds));
for m = 1:length(M)
    M1 = M(m);
    for d = 1:length(ds)
        delta = ds(d);
        Betas(m, d) = dbmBeta(delta, M1);
    end
end

Betas(Betas ~= real(Betas)) = 0;

maxThetaS = zeros(size(M));
figure, hold on
for m = 1:length(M)
    [~, jj] = find(Betas(m, :));
    plot(ds(jj), Betas(m, jj), '.-', 'MarkerSize', 5, 'DisplayName', ['M = ' char(string(M(m)))]);
    maxThetaS(m) = max(ds(jj));
    fprintf('Max wedge angle for Mach = %.2f is:\t%.2f\n', M(m), maxThetaS(m));
end
xlabel('Wedge Angle')
ylabel('Shock Angle')
ylim([0, 75])
xlim([0, 60])
legend
        
        
        

