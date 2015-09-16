% August 23rd 2010
% Rishidev Chaudhuri
%
%
% Run a single trial of the two population reduced mean-field model
% See reducedMF.m for how to set up parameters
% Note that there is no recurrent AMPA in this model

function result=runSingleSim(params,stimulus,nonStimCurrents)

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
   