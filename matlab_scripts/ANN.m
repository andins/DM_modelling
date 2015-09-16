function [popA, popB, traceTimes] = ANN(params, stimulus)

% runSingleStim accepts a function that takes time as input and returns a
% current. This is used for all the non-stimulus currents (useful to simulate targets).
% Of course, this could be changed to incorporate the stimulus current too,
% but then the stimulus noise would need to be handled separately.

% In the simple no-target condition, this can be implemented as an
% anonymous function returning a constant. More generally, you can spin it
% off into a separate file and pass the function handle
Ibg=.3255; % Background input current 
nonStimCurrents=@(currTime) Ibg;

singleTrialTrace=runSingleSim(params,stimulus,nonStimCurrents);
popA=binTraces(singleTrialTrace(1,:),params.binWidth); 
popB=binTraces(singleTrialTrace(2,:),params.binWidth);


if ~isnumeric(params.trialLen)
    traceTimes=binTraces(params.dt:params.dt:params.maxTrialLen,params.binWidth); % These are the times at which the traces were stored (i.e. the centers of the time bins)
else
    traceTimes=binTraces(params.dt:params.dt:params.trialLen,params.binWidth); % These are the times at which the traces were stored (i.e. the centers of the time bins)
end
