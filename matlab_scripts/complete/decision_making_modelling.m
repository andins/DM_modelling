% Andrea Insabato
% 17/09/15
% email: andrea.insabato@upf.edu

%% DDM integration

mu = .01;  % drift rate
sig = 0.2;  % noise
b = 5;  % boundary
duration = 'free'; % can be 'free' or a numeric value in ms

[v, time_steps] = DDM(mu, sig, b, duration);

figure(10), hold on
plot(time_steps, v, 'k', 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')
plot([time_steps(1) time_steps(end)], [-b -b],'--r')

%% DDM prob distr
clear all

mu = .0;
sig = 1;
b = 5;
num_trials = 1000;
for t=1:num_trials
    [v(:,t), time_steps] = DDM(mu, sig, b, 100);
end
% create a grid on dimensions v and time
ll = min(v(:)); 
ul = max(v(:))+0.1;
vpoints = linspace(ll,ul,100);
hpoints = linspace(time_steps(1), time_steps(end)+0.1, 200);
pv = zeros(length(vpoints)-1, length(hpoints)-1);
exact_pv = zeros(size(pv));

% calc P(v,t) for each point on the grid
for i=1:length(hpoints)-1
    % flatten v array
    v_flat = reshape(v(time_steps>=hpoints(i) & time_steps< hpoints(i+1),:),1,[]);
    % calc exact analytical probability
    exact_pv(:,i) = normpdf(vpoints(1:end-1), mu*hpoints(i), sig*sqrt(hpoints(i)));
    % calculate prob numerically
    for j=1:length(vpoints)-1
        pv(j,i) = sum(v_flat>= vpoints(j) & v_flat<vpoints(j+1))/num_trials;        
    end
end
figure(11)
hold on
% plot numerical prob with color 
imagesc(hpoints, vpoints, pv)
colormap jet
% and analytical prob with contour
[X, Y] = meshgrid(hpoints(1:end-1), vpoints(1:end-1));
contour(X,Y,exact_pv, 10,'color', [.7 1 .7])
xlim([hpoints(1), hpoints(end)])
ylim([vpoints(1), vpoints(end)])

%% create a synthetic dataset
clear all
clc

% stimuli statistics

% mean evidence supplied by stimulus for decision 1
e1_mu = [0, 1, 2, 4];
% mean evidence supplied by stimulus for decision 2
e2_mu = 0;
% std of evidence supplied by stimulus for decision 1
e1_sigma = ones(size(e1_mu))*0.3;
% std of evidence supplied by stimulus for decision 2
e2_sigma = 0.3;

figure(1), hold on
cols = get(gca,'colororder');
xx = linspace(-1,5,100);
% plot distributions of evidence
plot(xx,normpdf(xx, e2_mu, e2_sigma),'color', cols(1,:),'linewidth',2)
for i=1:length(e1_mu)
    plot(xx,normpdf(xx, e1_mu(i), e1_sigma(i)),'--','color',cols(i+1,:),'linewidth',2)
end

% choices
num_trials = 500; % for each condition
% theorical probability of correct choices
prob_correct = [0.5, 0.75, 0.9, 0.98];
% create choices as Bernoulli trials
for i=1:length(e1_mu)
    choices(i,:) = binornd(ones(1,num_trials),prob_correct(i));
end
% prob of correct in this sample
p_corr = sum(choices,2)/num_trials;
% std of prob of correct (according to binomial distr)
std_p_corr = sqrt(p_corr .* (1-p_corr) ./ num_trials);
figure(2)
subplot(1,2,1), hold on
errorbar(e1_mu, p_corr, std_p_corr,'ok','markersize',10,'linestyle','none')

% create RTs as Gaussian Mixtures
% one component will correspond to correct trials and the other to errors
% define some arbitrary (but plausible) values for mean and std
meanRT = [900, 700, 600, 500;
          900, 800, 700, 600];  % ms
stdRT = [700, 680, 600, 570];  % ms
mixing = prob_correct;  % mixing coefficients
figure(3)
for i=1:length(e1_mu)
    obj = gmdistribution(meanRT(:,i),stdRT(i), [mixing(i), 1-mixing(i)]);
    RTs(i,:) = obj.random(num_trials);
    subplot(2,2,i)
    hist(RTs(i,:),20)
end
figure(2)
subplot(1,2,2), hold on
errorbar(e1_mu, mean(RTs,2), std(RTs,1,2),'ok','markersize',10,'linestyle','none')
%% analytical error rate and RTs of DDM
% the input gain
a = .015;
% drift rate
mu = a * (e1_mu + .001);
% noise
sig = .7;
b = 20;  % boundary
t_nd = 100;  % non decision time
err_rate = 1./ (1+exp(2*mu.*b./sig.^2));  % error rate
fpt = (b./mu .* tanh(mu.*b./sig.^2) + t_nd);  % first passage time (RT)
figure(2)
subplot(1,2,1)
plot(e1_mu, 1-err_rate, 'r')
subplot(1,2,2)
plot(e1_mu, fpt, '--r')

%% numerical error rate and RTs of DDM
num_trials = 500;
% init vars
Rts = zeros(num_trials, length(mu));
choice = zeros(num_trials, length(mu));

for m=1:length(mu)
    for t=1:num_trials
        [v, time_steps] = DDM(mu(m), sig, b, 'free');
        Rts(t,m) = time_steps(end) + t_nd;  % find RT
        choice(t,m) = sign(v(end));  % find choice
    end
end
p_corr = sum(choice==1)/num_trials;
mRt = mean(Rts);
figure(2)
subplot(1,2,1)
plot(e1_mu, p_corr, '*r', 'linestyle','none')
subplot(1,2,2)
plot(e1_mu, mRt, '*r', 'linestyle','none')

%% supermodel
[model opt_params] = supermodel(e1_mu, RTs);  % create model and optimize params 
plot(e1_mu, model(opt_params, e1_mu),'--g')
%% supermodel validation
validation_set = sort(rand(1, 100)*(e1_mu(end)-e1_mu(1)) + e1_mu(1));
plot(validation_set, model(opt_params, validation_set), 'm')
mu = a * (validation_set + .001);  % drift rate associated to validation set
plot(validation_set,  (b./mu .* tanh(mu.*b./sig.^2) + t_nd), 'r')  % plot FPT for validation set

%% predicted RT distributions for correct and error trials
clear Rts choice  % clear some vars
num_trials = 1000;
drift = mu(2);  % take an arbitrary drift

for t=1:num_trials
    [v, time_steps] = DDM(drift, sig, b, 'free');
    Rts(t) = time_steps(end) + t_nd;
    choice(t) = sign(v(end));
end

Rt_corr = Rts(choice==1);
Rt_err = Rts(choice==-1);
mRt_corr = mean(Rt_corr);
mRt_err = mean(Rt_err);
figure, hold on
[f,x] = hist(Rt_corr,20);
stairs(x, f/trapz(x,f),'b','linewidth',2)
[f,x] = hist(Rt_err,20);
stairs(x, f/trapz(x,f),'r','linewidth',2)
plot([mRt_corr mRt_corr], get(gca,'ylim'),'b--','linewidth',2)
plot([mRt_err mRt_err], get(gca,'ylim'),'r--','linewidth',2)

%% variability in z and mu

mu = .02;  % mean drift rate
sig = 0.2;  % noise
b = 5;  % boundary
duration = 1000;  % ms
range_z = [-0, 0];  % range for initial point distribution (uniform distr)
sig_mu = .01;  % std for drift (normal distr)

[v, time_steps] = DDM_extended(mu, sig_mu, sig, b, range_z, duration);

figure(20), hold on
plot(time_steps, v, 'k', 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')
plot([time_steps(1) time_steps(end)], [-b -b],'--r')

%% predicted RT distributions for correct and error trials (extended model)
clear Rts choice
num_trials = 1000;
drift = 0.02;
range_z = [-0, 0];
sig_mu = .01;

for t=1:num_trials
    [v, time_steps] = DDM_extended(drift, sig_mu, sig, b, range_z, 'free');
    Rts(t) = time_steps(end) + t_nd;
    choice(t) = sign(v(end));
end

Rt_corr = Rts(choice==1);
Rt_err = Rts(choice==-1);
mRt_corr = mean(Rt_corr);
mRt_err = mean(Rt_err);
figure, hold on
[f,x] = hist(Rt_corr,20);
stairs(x, f/trapz(x,f),'b','linewidth',2)
[f,x] = hist(Rt_err,20);
stairs(x, f/trapz(x,f),'r','linewidth',2)
plot([mRt_corr mRt_corr], get(gca,'ylim'),'b--','linewidth',2)
plot([mRt_err mRt_err], get(gca,'ylim'),'r--','linewidth',2)

%% race model
e = [.0, .0];  % evidence to each integrator
sig = 0.2;  % noise
b = 5;  % boundary
duration = 1000;  % ms

[v, time_steps] = race_model(e, sig, b, duration);

figure, hold on
plot(time_steps, v, 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')

%% continuum race/DDM

e = [.01, -.01];  % drift for each integrator
sig = 0.2;  % noise
b = 5;  % boundary
duration = 'free';  % ms
rho = 0;  % correlation
nu = -1;  % sign of correlation

[v, time_steps] = DDM_race(e, sig, b, rho, nu, duration);

figure
subplot(1,2,1), hold on
plot(time_steps, v, 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')
plot([time_steps(1) time_steps(end)], [-b -b],'--r')
subplot(1,2,2), hold on
plot(time_steps, (v(1,:)-v(2,:))/2, 'linewidth', 2)
plot(time_steps, v(1,:)+v(2,:), 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')
plot([time_steps(1) time_steps(end)], [-b -b],'--r')

%% ANN: rate model
clear all
% define some parameters
run('wong_params.m')
% run one trial
[popA, popB, time_steps] = ANN(params, stimulus);

figure, hold on
plot(time_steps, popA(:,:),'b','linewidth',2)
plot(time_steps, popB(:,:),'r','linewidth',2)
% limit axis to avoid showing 0 values
if strcmp(params.trialLen, 'free')
    xlim([0 time_steps(find(popA==0,1,'first')-1)])
end

%% P(corr) and RTs with ANN
clear all
num_trials = 100;
t_nd = 0;  % non decision time
e1_mu = [0, 1, 2, 4];  % mean evidence for stimulus 1
% NOTE: for this one set params.trialLen to 'free'
run('wong_params.m')  

Rts = zeros(num_trials, length(e1_mu));
choice = zeros(num_trials, length(e1_mu));

for c=1:length(e1_mu)
    stimulus.mu = e1_mu(c);  % use the corresponding evidence value
    for t=1:num_trials
        [popA, popB, time_steps] = ANN(params, stimulus);
        % find RT
        last_step = find(popA==0,1,'first')-1;
        Rts(t,c) = time_steps(last_step) + t_nd;
        % find choice
        if ~isempty(last_step)
            if popA(last_step) > popB(last_step)
                choice(t,c) = 1;
            elseif popB(last_step) > popA(last_step)
                choice(t,c) = -1;
            elseif popA(last_step) == popB(last_step)  % undecided trials
                choice(t,c) = 0;
            end
        end
    end
end
p_corr = sum(choice==1)/num_trials;
mRt = mean(Rts);
figure(2)
subplot(1,2,1)
plot(e1_mu, p_corr, 'sb', 'linestyle','none')
subplot(1,2,2)
plot(e1_mu, mRt*1000, 'sb', 'linestyle','none')
Rt_corr = Rts(choice==1);
Rt_err = Rts(choice==-1);
mRt_corr = mean(Rt_corr);
mRt_err = mean(Rt_err);
figure, hold on
[f,x] = hist(Rt_corr,20);
stairs(x, f/trapz(x,f),'b','linewidth',2)
[f,x] = hist(Rt_err,20);
stairs(x, f/trapz(x,f),'r','linewidth',2)
plot([mRt_corr mRt_corr], get(gca,'ylim'),'b--','linewidth',2)
plot([mRt_err mRt_err], get(gca,'ylim'),'r--','linewidth',2)

%% compare dynamics of RM/DDM/ANN
% DDM
num_trials = 1;
e = [.0, -.0];
sig = 0.2;
b = 5;
rho = 1.;
nu = -1;

figure
subplot(1,3,1), hold on
for t=1:num_trials
    [v, time_steps] = DDM_race(e, sig, b, rho, nu, 'free');
    % plot the dynamics in the v1/v2 plane
    plot_2Dtimeseries(v(1,:), v(2,:),time_steps)
end

% RM
subplot(1,3,2), hold on
rho = 0.;
for t=1:num_trials
    [v, time_steps] = DDM_race(e, sig, b, rho, nu, 'free');
    plot the dynamics in the v1/v2 plane
    plot_2Dtimeseries(v(1,:), v(2,:),time_steps)
end

% ANN
subplot(1,3,3), hold on
run('wong_params.m')
params.trialLen = 3.;

for i=1:num_trials
    [popA, popB, time_steps] = ANN(params, stimulus);
    plot the dynamics in the v1/v2 plane
    plot_2Dtimeseries(popA, popB, time_steps)
end


colormap('jet')
