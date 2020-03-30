function err = coneError(thetaShock, thetaC, M1, gamma)
    % calculate the error between a cone angle that would correspond to the given shock angle and the given shock angle
    % thetaShock is an array of shock angles to calculate the error for
    % thetaC is the desired cone angle
    % M1 is the Mach number
    if nargin < 4   % this is just so that we could use another gamma if wanted for some reason
        gamma = 1.4;
    end
    % define some scalars for angle conversion and stuffs
    r2d = 180/pi;
    d2r = pi/180;
    % scale my next lower bound by this to make it a little lower so that its less likely to still be too high
    lowScaler = 2/3;    % totally arbitrary, just chosen so i dont have to redo too many times

    % loop through the given angles in thetaShock
    % this is so that if you give an array of shock angles that you get an array out rather than having to call multiple times
    for thetaInd = 1:length(thetaShock)
        thetaS = thetaShock(thetaInd);  % number that holds the current iterations shock angle
        d2 = dbmDelta(thetaS, M1);      % delta angle corresponding to this shock angle
        M1n = M1 .* sind(thetaS);       % find mach of normal shock component
        M2n = MAcrossShock(M1n);        % find normal mach behind shock
        M2 = M2n ./ sind(thetaS - d2);  % find the mach behind shock

        % the mach directly behind the shock can be found with normal and wedge relations done above
        % the normalized velocity can be found from this eq
        vPrimes = (2 / ((gamma - 1) * M2^2) + 1) ^ (-1/2);
        % find the radial and angular velocity components from the difference in the shock angle and delta
        vrPrimeS = vPrimes * cosd(thetaS - d2);
        vthetaPrimeS = vPrimes * -sind(thetaS - d2);

        % this is the function that describes the ODE that describes this flow
        % standard: derivative of first element is the second
        % solved rest of eq for the second part
        % because of another relation, we know that the derivative of vr is vtheta
        % NOTE if i was smart i wouldnt define this function in the middle of the loop for no reason
        f = @(theta, vrp) [vrp(2); ((gamma-1)/2 * (vrp(2).^2 + vrp(1).^2 -1) .* (2*vrp(1) + vrp(2) .* cot(theta)) + vrp(2).^2 .* vrp(1)) ./ ((gamma-1)/2 * (1 - vrp(1).^2 - vrp(2).^2) - vrp(2).^2)];
        lowerBound = thetaC * lowScaler * d2r;    % arbitrary, not too low hopefully where things get too weird, but hopefully low enough to make sure we dont need to recalculate
        theta0 = -1;                        % arbitrary, less than lowerBound to make sure itll loop
        % loop thought while the cone angle we find from the ODE is less than the lower bound the ODE was calculated to
        % this is so that itll keep looping though until it finds this angle in the given data (if didnt calculate far enough then need to recalc)
        % if there is no solution it will exit later
        while theta0 < lowerBound

            % use ode45 to calculate the solution to the ODE f starting at the shock angle with the calculated velocities
            % going to lowerBound which is probably a little less than the cone angle
            % The overall goal is to find where the theta component of velocity goes to 0 becasue that would be where the cone angle should be 
            % then we can try again based on how far its off
            % NOTE this function is in radians but most other stuff is in degrees cause why not
            [thetas, vPrimes] = ode45(f, [thetaS*d2r, lowerBound], [vrPrimeS; vthetaPrimeS]);
            % TODO try on the next iterations just starting where we left off instead of calculating the whole thing again.  i dont know if that works but it could cut down on time

            % NOTE there are many ways to do this and this is a bad one but its accurate so maybe its more stable while its way slower 
            % jk i actually dont know that its more accurate
            % but what i do is instead of going halfway between the points that cross 0 or soemthing 
            % this does a poly interpolationg for the 4 points around (2 on each side of) 0
            vt = vPrimes(:, 2);                 % the theta velocities
            notBadInds = find(isfinite(vt));    % find the real number indeces
            vt = vt(notBadInds);                % then redefine vt and t (theta angles) based on only existing ones
            t = thetas(notBadInds);
            % check if we found an exact solution by some fluke and solve if thats a thing
            ind = find(vt == 0, 1);
            if ~isempty(ind)
                theta0 = thetas(ind);
            else    % when it wont be
                % find indeces that are 2 above and below the 0 point of vTheta
                inds = [find(vt < 0, 2, 'last'); find(vt > 0, 2)];
                [~, ~, mu] = polyfit(t(inds), vt(inds), length(inds)-1);    % get the std and mean of the t points to readjust for less error in the polyfit
                newT = (t - mu(1)) / mu(2);     % recalc t
                % get the polyfit for these points 
                % NOTE this is indented to be 4 points with the 0 in the middle but it can be less
                % it can be 3 if the lower bound is just barely low enough that it crosses 0 (thats ok) or an error from the bad behavior of the ODE
                % can be 2 points if doesnt cross 0 at all
                pf = polyfit(newT(inds), vt(inds), length(inds)-1); 
                if ~all(isfinite(pf))   % if all the points arent finite then just use the middleish point
                    rs = t(inds(length(inds) / 2 + 1));
                else
                    % otherwise calculate the root of that poly
                    rs = roots(pf);
                    rs = rs * mu(2) + mu(1);    % rescale cause it was scaled my mu before
                    rs = rs(rs == real(rs));    % the imaginary stuff isnt helpful
                end
                % if there were more than 2 points for the fit 
                if length(inds) > 2
                    %set the calculated cone angle as the root between the found angles (there should be 3 roots but only 1 is helpfule)
                    theta0 = rs((rs >= t(inds(3))) & (rs <= t(inds(2))));
                    % if by fluke we get more than 1 root between the points surrounding 0 then just pic the second one
                    if length(theta0) > 1
                        theta0 = theta0(2);
                    % if we dont find an angle at all (computation errors can add up for these small numbers and get roots not in between at all)
                    elseif length(theta0) < 1
                        % then set the angle as the mean of the points around the 0 (i think this will do that)
                        theta0 = mean(t(inds(floor(length(inds) / 2 + 1):floor(length(inds) * 3 / 4))));
                    end
                % if there were only 2 points to fit (fit is linear) 
                % and the root found is bigger than an arbitrary small number (so that we dont keep looping forever as it approches 0 finding nothing)
                % and as the lower bound is greater than that small num (if its less then weve checked far enough its time to go probably)
                elseif rs > eps && lowerBound > eps
                    % now calculate a quad fit over all the data to get a feel for if there will be a 0 or if this was just a bad guess
                    [~, ~, mu] = polyfit(t, vt, 2);     % get stats to even up thetas
                    newT = (t - mu(1)) / mu(2);         % even up thetas
                    pf = polyfit(newT, vt, 2);          % get quadratic polyfit over all thetas calculated
                    allrs = roots(pf);                  % find the roots of that polyfit
                    allrs = allrs * mu(2) + mu(1);      % rescale so they mean somethign
                    % weed out the roots that we dont want hopefully
                    allr = max(allrs(allrs == real(allrs)));    % pick the bigger one (fairly arbitrary but if we just havent found the 0 because we didnt calculate far enough then it might find a root close and a negative one far away so i picked the positive one)
                    % set a temp lower bound as the bigger of eps as a small num used above and the 2 options we have for roots now
                    % want to set a small lower bound but if its too small then we waste time computing a lot of poorly behaved ODE (very expensive) 
                    lb = max(eps, min(rs, allr));
                    % if the quad fit shows there may be a 0 and we got it less than the current lower bound
                    if ~isempty(allr) && lb < lowerBound
                        % then set the lower bound from lb
                        lowerBound = lb * lowScaler - eps;      % added -eps so that the points where thetaC is 0 do get less than 0 cause maybe thats important
                        theta0 = -1;    % just make sure its less so we loop again
                    else    % otherwise we dont know if there will be anything but well try again vainly based off eps (last try)
                        lowerBound = eps * lowScaler;
                        theta0 = -1;
                    end
                else
                    % if we didnt find anything ever then set the guess for the cone angle to be the average of theta component of velocities
                    % this is in hopes that the velocities will be further from 0 the worse the guess (this is to make the function better behaved for fsolve to use)
%                     theta0 = -abs(rs);
                    theta0 = mean(vt);
                    lowerBound = theta0 -1;     % gotta make sure lowerBound is less so we dont loop again
                end
            end
        end
        err(thetaInd) = theta0*r2d - thetaC;  % give error back in degrees
    end
end




% i thought this should be more efficient but it doesnt give me conclusive results about that at the moment
% but its supposed to give the derivative but precalculate things for efficiency
function dv = dvdt(t, v)
    vt2 = v(2).^2;
    c = 0.2 * (1 - v(1).^2 - vt2);
    dv = [v(2); (-c .* (2*v(1) + v(2) .* cot(t)) + vt2 .* v(1)) ./ (c - vt2)];
end
