params.dt=0.2e-3; % Time step (s)
params.trialLen='free'; % Trial length (s) or 'free' to run until decision criterion is met
params.maxTrialLen=10.0; % maximum trial length if trialLen is 'free' (s)
params.thr = 35;  % decision criterion [Hz]

% Recurrent connections
w_r=1.;
w_i=1.;
params.J=[0.2609*w_r,-0.0497*w_i; -0.0497*w_i, 0.2609*w_r]; % nA

params.tau=.1; % seconds
params.nmdaGamma=0.641;
params.tauNoise=2e-3; % seconds
params.sigmaBG=0.02; % nA

% Stimulus equation as in Wong 2007. 
stimulus.mu=0;  % differential input to the decision pools
stimulus.inputGain=.1; % this modulates the gain of the stimulus strength 
stimulus.JAext= 5.2e-4; % nA/Hz
stimulus.mu0=40; % Hz
stimulus.s_n = params.sigmaBG^2; % fluctuations of the stimulus: if you want them iqual to the prestimulus period set it to params.sigmaBG^2
stimulus.sigma=sqrt(stimulus.s_n-params.sigmaBG^2); % nA. This allows to separate stimulus and background noise
stimulus.tOn=.5; % Stimulus onset time (s)

% We don't necessarily want to store the firing rates at the same
% resolution at which they were calculated. binWidth is the number of data
% points that should be averaged together for each recorded data point.
% Note that this is not a sliding window (windows are non overlapping!).
params.binWidth=50; 
    
