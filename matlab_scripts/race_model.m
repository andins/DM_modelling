function [v, time_steps] = race_model(e, sig, b, duration)

if all(size(e)==[2,1])
elseif all(size(e)==[1,2])
    e = e';
else
    error('e must be a two-elements vector!')
end

dt = 0.1;  % ms
time_steps = 0:dt:duration;
% init
v = zeros(2,length(time_steps));

for t=2:length(time_steps)  % integration
    v(:,t) = v(:,t-1) + dt*e + sqrt(dt)*randn(2,1)*sig;
end