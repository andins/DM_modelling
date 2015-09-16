function [v, time_steps] = DDM(drift, sig, b, duration)

% integration step
dt = 0.1;  % ms

% if duration is a numeric value integrate until the time is over
if isnumeric(duration)
    time_steps = 0:dt:duration;  % the time stamp in each step
    v = zeros(1,length(time_steps));  % init v

    for t=2:length(time_steps)  % integration (Euler)
        v(t) = v(t-1) + dt*drift + sqrt(dt)*randn()*sig;
    end
elseif strcmp(duration, 'free')  % if duration is 'free' integrate until the decision criterion is met
    t=1;  % first time step
    v(t)=0;  % first v
    time_steps(t) = 0;  % first time stamp
    while v(t)<=b && v(t)>=-b  % decision criterion
        t = t+1;  % increment t
        time_steps(t) = time_steps(t-1) + dt; % increment time_step
        v(t) = v(t-1) + dt*drift + sqrt(dt)*randn()*sig; % integration (Euler)     
    end
else  % raise an error if duration is set wrong
    error('duration must be either a number for fixed duration in ms or the string ''free'' to run simulatation until decision criterion is met')
end