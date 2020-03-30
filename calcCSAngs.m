% calcCSAngs.m
function [thetaS, offBy, fit] = calcCSAngs(M, tc, prevThetaS)
%     tol = .000001;
    guessIncBig   =   5;                        % how much to increment when the guess is not around the 0 point
    guessIncSmall = 0.5;                        % a small increment to add to a hopefully more precise guess (to make so we dont have to calculate again)
    fitOrder = 4;                               % the initial fit order to use for fitting
    maxFO = 8;                                  % stop fitting before using this order cause thats gonna just be weird
    bestFO = [fitOrder, length(prevThetaS)];    % store the previous best fit order as well as the number of points that are outside the std bounds
    z = 3;                                      % the scalar of the std of the data to compare points to
    proportions = [1/3, 1];   % should be 2 percentages between mu and Beta for the initial guess to surround the correct shock angle
    thetaS = prevThetaS;    % this is a waste but it makes sense to me to have the input var be prev but i dont want to use prev anywheres else
    offBy  = zeros(length(M), length(tc));      % initialize teh offBy variable (stores the value of the function where it should be 0)
    % loop through the mach numbers
    for m = 1:length(M)
        M1 = M(m);  % set M1 for this iteration
        done = 0;   % we arent done finding the angles for this iter
        [~, j, ts] = find(thetaS(m, :));    % find the values and indeces were there is a defined thetaS
        if isempty(ts)  % if its empty then we havent calculated anything yet so we need to calculate everything
            toDo = ones(size(thetaS, 2));   % set toDo to say calc everything
            noFit = 1;                      % we arent gonna do any fitting the first found cause theres no data
            fit = zeros(size(thetaS));
        else    % if we do have previous data
            % The goal here is that we get a general fit of the data and then compare the points that are too far from the fit and only recalculate those ones
            % theres some logic used that i only want to recalculate an arbitrary 5 points at once, otherwise i dont trust the fit that much to say more points are wrong
            numPointsToRecalc = 5;
            % calculate a polynomial fit of the inverse of the shock angles found (onyl for the non-zero points) of order fitOrder
            fit(m, :) = 1./polyval(polyfit(tc(j), 1./ts, fitOrder), tc);    % im using the reciprocal of the data because that seemed to straighten it out more than raw so the ends behaved better
            syx = stdOfFit(ts, fit(m, j));      % find std of fit to compare
            diff = thetaS(m, j) - fit(m, j);    % get the differences between the showck angles and the fit
            toDo = abs(diff) > syx * z;         % set to do the ones that are further than z times the std
            % loop through tightening the bounds around the fit were checking untill we have more than 5 points that need to be recalculated
            while sum(toDo) < numPointsToRecalc && z > 1
                z = z - 1;  % decrement z to tighten
                % recalculate the fit, std, and toDo for tighter bounds
                fit(m, :) = 1./polyval(polyfit(tc(j), 1./ts, fitOrder), tc);
                syx = stdOfFit(ts, fit(m, j));
                diff = thetaS(m, j) - fit(m, j);
                toDo = abs(diff) > syx * z;
            end
            % loop through using a tighter fit untill we have 5+ points that need to be recalculated
            while sum(toDo) >= numPointsToRecalc && fitOrder < maxFO
                fitOrder = fitOrder + 1;    % increment fit order for tighter fit
                % recalculate the fit, std, and toDo for tighter bounds
                fit(m, :) = 1./polyval(polyfit(tc(j), 1./ts, fitOrder), tc);
                syx = stdOfFit(ts, fit(m, j));
                diff = thetaS(m, j) - fit(m, j);
                toDo = abs(diff) > syx * z;
                % if this fit has more points to redo than previous best then reset best to this fit
                if sum(toDo) > bestFO(2)
                    bestFO = [fitOrder, sum(toDo)];
                end
            end
            % finally calculate the fit etc again after all that fit moving stuffs
            fit(m, :) = 1./polyval(polyfit(tc(j), 1./ts, bestFO(1)), tc);
%             figure
%             hold on
%             plot(tc(j), fit(m, j), '-.', 'MarkerSize', 6);
%             plot(tc, thetaS(m, :), '.', 'MarkerSize', 6);
            syx = stdOfFit(ts, fit(m, j));
            diff = thetaS(m, j) - fit(m, j);
            toDo = abs(diff) > syx * z;     % do the ones where diff is too big
            toDo(thetaS(m, :) == 0) = 0;    % dont redo points that were 0 becuase they dont exist
                            
            noFit = 0;      % we do want to base things off the fit
        end
        % loop through the cone angles to calculate the shock angles
        for c = 1:length(tc)
            if ~toDo(c)     % go to next iteration if we dont need to do this point again
                continue
            end
            thetaC = tc(c); % set the cone angle variable for this iteration

            if noFit    % if not based on fit base init guess on mu and Beta
                mu = asind(1/M1);
                Beta = dbmBeta(thetaC, M1);
                % Beta for wedges stops existign before the shock angle for cones so itll give complex answeres
                if real(Beta) ~= Beta
                    % past existance of beta then just use the last result + the increment as upper bound for init guess
                    Beta = thetaS(m, c-1) + guessIncBig;
                end

                guess = Beta * proportions + mu * (1 - proportions);    % set guess array with lower bound from mu and upper from Beta or last sol
            else
                % if basing off the fit then guess based off the point at the fit and the point on the other side of the fit from the data point
                % this is to sorta force the point thats too far from the fit to go up to the fit 
                guess = sort([fit(m, c), fit(m, c) - sign(diff(c)) * guessIncBig]);     % sort it cause i want it in order and dont want to use logic
            end
            % calculate the cone angle error from these guesses and then add a point in the middle that will just be the upper bound for now
            cE = coneError(guess, thetaC, M1);
            guess(3) = guess(2);
            cE(3) = cE(2);
            firstTime = 1;  % neccessary for loop
            % loop through while either the lower bound isnt negative or the upper isnt positive since we know that if the guessed angle is too high it should give too high of an output
            while (cE(1) > 0 || cE(3) < 0)
                % if both points are wrong 
                if cE(1) > 0 && cE(3) < 0
                    % then make a new middle point in the center or the bounds so that it doesnt go back and forth on old points (hopefully this data point gives us something better)
                    guess(2) = mean(guess([1, 3]));
                    cE(2) = coneError(guess(2), thetaC, M1);    % get error at that point
                    guessIncSign = -sign(cE(2));                % set the direction as the one opposite the middle point (if the middle is neg then move the upper) to avoid errors swapping points
                % if the bottom is pos then reset the bottom, otherwise move the upper
                elseif cE(1) > 0    
                    guessIncSign = -1;
                else
                    guessIncSign = 1;
                end
                % set which bound to move
                bound = 2 + guessIncSign;
                notBound = 4 - bound;
                % if its the first iter (with garbage middle point) 
                % or the line between the bound and middle point crosses 0 on the other side of the bound 
                % (eg if its the lower bound and the middle point is less than the lower that indicates that the 0 would be higher not lower which isnt helpful)
                if firstTime || guessIncSign * (cE(bound) - cE(2)) < 0
                    % then just move the bound by the big increment

                    % first move further point closer (make the bounds tighter)
                    cE(notBound) = cE(2);
                    guess(notBound) = guess(2);
                    % then move the middle point to the bound
                    cE(2) = cE(bound);
                    guess(2) = guess(bound);
                    % reset the bound from its current state + the big increment
                    guess(bound) = guess(bound) + guessIncBig * guessIncSign;
                    firstTime = 0;  % now its not the first time
                else
                    % otherwise try to guess where the 0 will be based on the line between the middle point and this bound 
                    guessInc = - cE(bound) * (guess(2) - guess(bound)) / (cE(2) - cE(bound)) + guessIncSmall * guessIncSign;    % use the small inc to go further so that its hopefully on the other side
                    % if the move was bigger than big inc then just use big inc cause its probably a bad guess
                    if ~(abs(guessInc) < guessIncBig)
                        guessInc = guessIncBig * guessIncSign;
                    end                     
                    % move the other bound in and middle point to bound
                    cE(notBound) = cE(2);
                    guess(notBound) = guess(2);
                    cE(2) = cE(bound);
                    guess(2) = guess(bound);   
                    guess(bound) = guess(bound) + guessInc; % set bound's guess
                end
                % if had to guess more than 90 degrees then this point doesnt work so were done
                if guess(3) > 90
                    done = 1;
                    break
                end
                % calculate the error and this bounds new guess
                cE(bound) = coneError(guess(bound), thetaC, M1);
            end
            % if done then move on to the next mach numb
            if done
                break
            end

            % use 0 solver to find where the cone angle error is 0 (where the guess of shock angle is correct for this cone angle)
%             [thetaS(m, c), offBy(m, c)] = badZero(@(thetaS) coneError(thetaS, thetaC, M1), guess, tol);
            [thetaS(m, c), offBy(m, c)] = fzero(@(thetaS) coneError(thetaS, thetaC, M1), guess([1, 3]));
        end
    end
    % calculate the fit again (was not calculated for the very first time above)
    [~, j, ts] = find(thetaS(m, :));
    fit(m, :) = 1./polyval(polyfit(tc(j), 1./ts, bestFO(1)), tc);
end

% the standard deviation of the fit
function [syx, dif] = stdOfFit(y, yc)
    % y is the point, yc is the fit at that point
    dif = y - yc;
    syx = sqrt(sum(dif.^2) / (length(y) - 2));
end

% i should comment this but basically this is a bad 0 solver i wrote for this application
% it uses the 3 point guess array from above
% it uses a quadratic fit for the 3 guess points for where the 0 might be (similar to secant method) to move the bounds in
% because of error it has to fall back on bisection sometimes
% this is probably bad cause i dont know what im doing really and i calculate a poly fit every iteration which isnt great
function [x, fx] = badZero(func, guess, tol)
    fErr = func(guess);
    while (abs(range(guess)) > tol)
        rs = roots(polyfit(guess, fErr, 2));
        rs = rs(rs == real(rs));
        r = rs((rs >= guess(1)) & (rs <= guess(3)));
        if length(r) ~= 1
            r = mean(guess([1, 3]));
            cond = fErr(2) * fErr(1) <= 0;
            if xor(cond, r > guess(2))
                temp = r; r = guess(2); guess(2) = temp;
            end
        else
            cond = r > guess(2);
        end
        if cond
            guess(1) = guess(2);
            fErr(1) = fErr(2);
        else
            guess(3) = guess(2);
            fErr(3) = fErr(2);
        end
        guess(2) = r;
        fErr(2) = func(guess(2));
    end
    x = guess(2);
    fx = fErr(2);
end
