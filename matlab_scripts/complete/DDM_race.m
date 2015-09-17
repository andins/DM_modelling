function [v, time_steps] = DDM_race(e, sig, b, rho, nu, duration)


dt = 0.1;  % ms
if isnumeric(duration)
    % init
    time_steps = 0:dt:duration;
    v = zeros(2,length(time_steps));

    for t=2:length(time_steps)
        W = randn(1,3);  % noise instances
        % first integrator
        v(1,t) = v(1,t-1) + dt*e(1) + sig* sqrt(dt)* (sqrt(1-rho)*W(1) + sqrt(rho)*W(3));
        % second integrator
        v(2,t) = v(2,t-1) + dt*e(2) + sig* sqrt(dt)* (sqrt(1-rho)*W(2) + sqrt(rho)*nu*W(3));
    end
elseif strcmp(duration, 'free')
    % initial values
    t=1;
    v=ones(2,1).*0;
    time_steps(t) = 0;
    while v(1,t)<=b && v(1,t)>=-b && v(2,t)>=-b && v(2,t)<=b  % decision criterion
        t = t+1;
        time_steps(t) = time_steps(t-1) + dt; 
        W = randn(1,3);
        v(1,t) = v(1,t-1) + dt*e(1) + sig* sqrt(dt)* (sqrt(1-rho)*W(1) + sqrt(rho)*W(3));
        v(2,t) = v(2,t-1) + dt*e(2) + sig* sqrt(dt)* (sqrt(1-rho)*W(2) + sqrt(rho)*nu*W(3));
    end
else
    error('duration must be either a number for fixed duration in ms or the string ''free'' to run simulatation until decision criterion is met')
end
