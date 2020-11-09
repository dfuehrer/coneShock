%MfromCSAngs.m

function [M, offby] = MfromCSAngs(thetaC, thetaS)

    % vars
    proportions = [.4, .6];
    guessIncBig   =  .5;                        % how much to increment when the guess is not around the 0 point
    guessIncSmall = 0.1;                        % a small increment to add to a hopefully more precise guess (to make so we dont have to calculate again)

    % lower bound is Mach for delta=0 -> mu
    Md0 = 1 / sind(thetaS);
    % the oblique shock Mach is an upper bound
    Mo = dbmMach(thetaC, thetaS);
    % Beta for wedges stops existign before the shock angle for cones so itll give complex answeres
    if ~isreal(Mo)
        % oblique shock detaches around 46 deg delta so wont find Mach
        % TODO figure out a good case here
        Mo = 1 / sind(thetaS) - thetaC -  guessIncBig;
    end

    guess = Mo * proportions + Md0 * (1 - proportions);    % set guess array with lower bound from mu and upper from Beta or last sol
    cE(1) = coneError(thetaS, thetaC, guess(1));
    cE(2) = coneError(thetaS, thetaC, guess(2));
    guess(3) = guess(2);
    cE(3) = cE(2);
    firstTime = 1;  % neccessary for loop
    % loop through while either the lower bound isnt negative or the upper isnt positive since we know that if the guessed angle is too high it should give too high of an output
    while (cE(1) > 0 || cE(3) < 0)
        % if both points are wrong 
        if cE(1) > 0 && cE(3) < 0
            % then make a new middle point in the center or the bounds so that it doesnt go back and forth on old points (hopefully this data point gives us something better)
            guess(2) = mean(guess([1, 3]));
            cE(2) = coneError(thetaS, thetaC, guess(2));    % get error at that point
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
        % calculate the error and this bounds new guess
        cE(bound) = coneError(thetaS, thetaC, guess(bound));
    end
    [M, offBy] = fzero(@(M) coneError(thetaS, thetaC, M), guess([1, 3]));
end
