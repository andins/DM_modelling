function [v, time_steps] = DDM_race(e, sig, b, rho, nu, duration)


dt = 0.1;  % ms
if isnumeric(duration)
    % init
    time_steps = 0:dt:duration;
    v = zeros(2,length(time_steps));

    for t=2:length(time_steps)
        W = randn(1,3);  % noise instances
        % first integrator
        v(1,t) = ;
        % second integrator
        v(2,t) = ;
    end
elseif strcmp(duration, 'free')
    % initial values
    t=;
    v=;
    time_steps(t) = ;
    while v(1,t)<=b && v(1,t)>=-b && v(2,t)>=-b && v(2,t)<=b  % decision criterion
        t = ;
        time_steps(t) = ; 
        W = randn(1,3);
        v(1,t) = ;
        v(2,t) = ;
    end
else
    error('duration must be either a number for fixed duration in ms or the string ''free'' to run simulatation until decision criterion is met')
end
