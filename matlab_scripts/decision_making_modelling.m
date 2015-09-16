% Andrea Insabato
% email: andrea.insabato@upf.edu

%% DDM integration

mu = .01;
sig = 0.2;
b = 5;
duration = 'free'; %1000;  % ms

[v, time_steps] = DDM_num(mu, sig, b, duration);

figure(10), hold on
plot(time_steps, v, 'k', 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')
plot([time_steps(1) time_steps(end)], [-b -b],'--r')

%% DDM prob distr
clear all

mu = .1;
sig = 1;
b = 5;
num_trials = 1000;
for t=1:num_trials
    [v(:,t), time_steps] = DDM_num(mu, sig, b, 100);
end
ll = min(v(:));
ul = max(v(:))+0.1;
vpoints = linspace(ll,ul,100);
hpoints = linspace(time_steps(1), time_steps(end)+0.1, 200);
pv = zeros(length(vpoints)-1, length(hpoints)-1);
exact_pv = zeros(size(pv));

for i=1:length(hpoints)-1
    v_flat = reshape(v(time_steps>=hpoints(i) & time_steps< hpoints(i+1),:),1,[]);
    exact_pv(:,i) = normpdf(vpoints(1:end-1), mu*hpoints(i), sig*sqrt(hpoints(i)));
    for j=1:length(vpoints)-1
        pv(j,i) = sum(v_flat>= vpoints(j) & v_flat<vpoints(j+1))/num_trials;        
    end
end
figure(11)
hold on
imagesc(hpoints, vpoints, pv)
colormap jet
[X, Y] = meshgrid(hpoints(1:end-1), vpoints(1:end-1));
contour(X,Y,exact_pv, 10,'color', [.7 1 .7])
xlim([hpoints(1), hpoints(end)])
ylim([vpoints(1), vpoints(end)])

%% create a synthetic dataset
clear all
clc

% stimuli statistics

e1_mu = [0, 1, 2, 3];
e2_mu = 0;
e1_sigma = ones(size(e1_mu))*0.3;
e2_sigma = 0.3;

figure(1), hold on
cols = get(gca,'colororder');
xx = linspace(-1,5,100);
plot(xx,normpdf(xx, e2_mu, e2_sigma),'color', cols(1,:),'linewidth',2)
for i=1:length(e1_mu)
    plot(xx,normpdf(xx, e1_mu(i), e1_sigma(i)),'--','color',cols(i+1,:),'linewidth',2)
end

% choices
num_trials = 500; % for each condition
prob_correct = [0.5, 0.75, 0.85, 0.95];
for i=1:length(e1_mu)
    choices(i,:) = binornd(ones(1,num_trials),prob_correct(i));
end
figure(2)
subplot(1,2,1), hold on
plot(e1_mu, sum(choices,2)/num_trials,'ok','markersize',10,'linestyle','none')
% RTs
meanRT = [900, 700, 600, 500;
          900, 800, 700, 600];  % ms
stdRT = [700, 680, 600, 570];  % ms
mixing = prob_correct;
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
a = .015;
mu = a * (e1_mu + .001);
sig = .7;
b = 20;
t_nd = 100;
err_rate = 1./ (1+exp(2*mu.*b./sig.^2));
fpt = (b./mu .* tanh(mu.*b./sig.^2) + t_nd);
figure(2)
subplot(1,2,1)
plot(e1_mu, 1-err_rate, 'r')
subplot(1,2,2)
plot(e1_mu, fpt, '--r')

%% numerical error rate and RTs of DDM
num_trials = 500;
Rts = zeros(num_trials, length(mu));
choice = zeros(num_trials, length(mu));

for m=1:length(mu)
    for t=1:num_trials
        [v, time_steps] = DDM_num(mu(m), sig, b, 'free');
        Rts(t,m) = time_steps(end) + t_nd;
        choice(t,m) = sign(v(end));
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
c = randn(1,7);
poly4 = @(c, x) c(7) + c(1)*x.^1 + c(2)*x.^2 + c(3)*x.^3 + c(4)*x.^4 + c(5)*x.^5 + c(6)*x.^6;
% rt_poly = a_1*mu + a_2*mu.^2 + a_3*mu.^3 + a_4*mu.^4 + a_5*mu.^5;
fun = @(c) sum((poly4(c,e1_mu)' - mean(RTs,2)).^2);
c_opt = fminunc(fun, c);
plot(e1_mu, poly4(c_opt, e1_mu),'--g')
validation_set = sort(rand(1, 100)*(e1_mu(end)-e1_mu(1)) + e1_mu(1));
plot(validation_set, poly4(c_opt, validation_set), 'm')
mu = a * (validation_set + .001);
plot(validation_set,  (b./mu .* tanh(mu.*b./sig.^2) + t_nd), 'r')

%% predicted RT distributions for correct and error trials
clear Rts choice
num_trials = 1000;
drift = mu(2);

for t=1:num_trials
    [v, time_steps] = DDM_num(drift, sig, b, 'free');
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

%% "real" RT distribution for correct and error trials



%% variability in z and mu

mu = .0;
sig = 0.2;
b = 5;
duration = 1000;  % ms
range_z = [-3, 3];
sig_mu = .01;

[v, time_steps] = DDM_extended(mu, sig_mu, sig, b, range_z, duration);

figure(20), hold on
plot(time_steps, v, 'k', 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')
plot([time_steps(1) time_steps(end)], [-b -b],'--r')

%% predicted RT distributions for correct and error trials (extended model)
clear Rts choice
num_trials = 1000;
drift = mu(20);
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
e = [.0, .0];
sig = 0.2;
b = 5;
duration = 1000;  % ms

[v, time_steps] = race_model(e, sig, b, duration);

figure, hold on
plot(time_steps, v, 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')

%% continuum race/DDM

e = [.0, -.0];
sig = 0.2;
b = 5;
duration = 1000;  % ms
rho = 1;
nu = -1;

[v, time_steps] = DDM_race(e, sig, b, rho, nu, duration);

figure, hold on
plot(time_steps, v, 'linewidth', 2)
plot([time_steps(1) time_steps(end)], [b b],'--r')
plot([time_steps(1) time_steps(end)], [-b -b],'--r')

%% ANN: rate model