function [v, time_steps] = DDM(drift, sig, b, duration)

% integration step
dt = 0.1;  % ms

% if duration is a numeric value integrate until the time is over
if isnumeric(duration)
    time_steps = 0:dt:duration;  % the time stamp in each step
    v = zeros(1,length(time_steps));  % init v

    for t=2:length(time_steps)  % integration (Euler)
        v(t) = ;
    end
elseif strcmp(duration, 'free')  % if duration is 'free' integrate until the decision criterion is met
    t=;  % first time step
    v(t)=;  % first v
    time_steps(t) = ;  % first time stamp
    while   % decision criterion
        t = ;  % increment t
        time_steps(t) = ; % increment time_step
        v(t) = ; % integration (Euler)     
    end
else  % raise an error if duration is set wrong
    error('duration must be either a number for fixed duration in ms or the string ''free'' to run simulatation until decision criterion is met')
end