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


end


% additional functions


function result=binTraces(inpTraces,binWidth)
% Takes a trace sampled at some higher resolution and reduces it by
% averaging over non-overlapping windows of size 'binWidth'
% Trace should be one-dimensional and have length divisible by binWidth
result=mean(reshape(inpTraces,binWidth,[]));
end


function result=stepNoise(Inoise,s,d,sigma)
    % August 23rd 2010
% Steps the filtered white noise forward
result=Inoise*s+d*sigma*randn(size(Inoise));
end    


function result=phi(I)
% August 23rd 2010
% Function for f-I curve
% Parameters for f-I curve, values in Hz
a=270.0;
b=108.0;
c=0.154; 

x=a*I-b;

if(abs(c*x)<1e-6)
    result=.5*x+1.0/c; % This is just for numerical stability
else
    result=x/(1-exp(-c*x));
end
end

function result=runSingleSim(params,stimulus,nonStimCurrents)
% August 23rd 2010
% Rishidev Chaudhuri
%
%
% Run a single trial of the two population reduced mean-field model
% See reducedMF.m for how to set up parameters
% Note that there is no recurrent AMPA in this model

% Assign some variables for readability
if isnumeric(params.trialLen)
    trialLen=params.trialLen;
elseif strcmp(params.trialLen, 'free')
    trialLen=params.maxTrialLen;
end
dt=params.dt;
tau=params.tau;
nmdaGamma=params.nmdaGamma;
thr = params.thr;

storedRates=zeros(2,round(trialLen/dt)); % Results go here

% Pre-calculate and store stimulus current
% If you want a time varying stimulus you could pass more parameters in the
% structure 'stimulus' and have the program calculate the currents at each
% step, or you could pass a stimulus function (similar to nonStimCurrents)
% or you could incorporate the stimulus into nonStimCurrents (and change
% the name)
stimulus.currents=zeros(2,1);
stimulus.currents(1)=stimulus.JAext*stimulus.mu0*(1+stimulus.inputGain*stimulus.mu);
stimulus.currents(2)=stimulus.JAext*stimulus.mu0*(1-stimulus.inputGain*stimulus.mu);

% Noise parameters
% noiseS and noiseD are used to exactly calculate the Ornstein-Uhlenbeck
% process that gives us Inoise

noiseS=exp(-dt/params.tauNoise); noiseD=sqrt(0.5*(1-noiseS^2)); 
sigma=params.sigmaBG;

% Variables used during the trial
s=zeros(2,1); currRates=zeros(2,1);
Inoise=zeros(2,1); Iext=zeros(2,1); Itotal=zeros(2,1);

currStep=1;

for currTime=dt:dt:trialLen
    
    % Find the external current
    IExt=nonStimCurrents(currTime);
    
    % Add in the stimulus, if needed
    if(currTime>=stimulus.tOn)
        IExt=IExt+stimulus.currents;
    end;
    
    % We also allow the stimulus to have some noise of its own, which can
    % be accounted for by increasing sigma
    if((currTime>=stimulus.tOn)&&(currTime<(stimulus.tOn+dt)))
        sigma=sqrt(params.sigmaBG^2+stimulus.sigma^2);
    end;
    
    % This function just updates a noise current as a zero-mean
    % Ornstein-Uhlenbeck process
    Inoise=stepNoise(Inoise,noiseS,noiseD,sigma);
    
    Itotal=params.J*s+Inoise+IExt;
    
    % Phi is a separate function containing the f-I relationship
    currRates(1)=phi(Itotal(1));
    currRates(2)=phi(Itotal(2));
    
    % Update the NMDA gating variable
    % This equation is the heart of the model 
    s=s+dt*(-s/tau+nmdaGamma*currRates.*(1-s));
    
    storedRates(:,currStep)=currRates;
    currStep=currStep+1;    
    if strcmp(params.trialLen, 'free')  % only check dec criterion if required by user
        if any(currRates >= thr)  % if dec criterion is met return
            break
        end
    end
end

result=storedRates;
end