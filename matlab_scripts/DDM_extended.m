function [v, time_steps] = DDM_extended(mu, sig_mu, sig, b, range_z, duration)

z = rand() * (range_z(2)-range_z(1)) + range_z(1);  % sample initial point for v from uniform
drift = randn()*sig_mu + mu;  % sample drift from normal
dt = 0.1;  % ms

if isnumeric(duration)
    time_steps = 0:dt:duration;
    v = zeros(1,length(time_steps));
    v(1) = z;
    for t=2:length(time_steps)
        v(t) = v(t-1) + dt*drift + sqrt(dt)*randn()*sig;
    end
elseif strcmp(duration, 'free')
    t=1;
    v(t)=z;
    time_steps(t) = 0;
    while v(t)<=b && v(t)>=-b
        t = t+1;
        time_steps(t) = time_steps(t-1) + dt; 
        v(t) = v(t-1) + dt*drift + sqrt(dt)*randn()*sig;
    end
else
    error('duration must be either a number for fixed duration in ms or the string ''free'' to run simulatation until decision criterion is met')
end