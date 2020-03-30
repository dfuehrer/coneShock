clear
M = [1.5, 2, 5];    % list of mach numbers to run
tc = 0:55;          % list of cone angles to calculate
maxIter = 5;        % number of times to run through fixing points
iter = 1;
tic     % remove eventually
% calculate the shock angles for the mach and cone angles listed above
[thetaS, offby, fit] = calcCSAngs(M, tc, zeros(length(M), length(tc)));
elapTime(1) = toc
% plot out what we get
figure
hold on
for m = 1:length(M)
    plot(tc, thetaS(m, :), '.', 'MarkerSize', 5, 'DisplayName', ['M = ' char(string(M(m)))]);
end
legend
% loop through readjusting points that are outside the norm
while any(offby, 'all') && (iter < maxIter)     % loop while still making changes and have tried to fix under 5 times
    tic
    % calculate things again
    [thetaS, offby, fit] = calcCSAngs(M, tc, thetaS);
    elapTime(iter+1) = toc
    % plot everything out again with the fit
    figure
    hold on
    for m = 1:length(M)
        [~, j] = find(thetaS(m, :));
        plot(tc(j), fit(m, j), '-.', 'MarkerSize', 5, 'DisplayName', ['fit ' char(string(M(m)))]);
        plot(tc(j), thetaS(m, j), '.', 'MarkerSize', 6, 'DisplayName', ['M = ' char(string(M(m)))]);
    end
    legend
    iter = iter + 1;        % iterate 
end

