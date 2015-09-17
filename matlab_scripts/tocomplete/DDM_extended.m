function [v, time_steps] = DDM_extended(mu, sig_mu, sig, b, range_z, duration)

z = ;  % sample initial point for v from uniform
drift = ;  % sample drift from normal
dt = 0.1;  % ms

if isnumeric(duration)
    time_steps = 0:dt:duration;
    v = zeros(1,length(time_steps));
    v(1) = z;
    for t=2:length(time_steps)
        v(t) = ;
    end
elseif strcmp(duration, 'free')
    t=1;
    v(t)=z;
    time_steps(t) = 0;
    while 
        t = ;
        time_steps(t) = ; 
        v(t) = ;
    end
else
    error('duration must be either a number for fixed duration in ms or the string ''free'' to run simulatation until decision criterion is met')
end