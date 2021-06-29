function [rtot,xcorr1,xcorr2,partCorr,xpartCorr] = Kuramoto_ode(n_osc,K,c)

%% Set the initial Kuramoto model parameters
options = odeset('MaxStep',.02); % ODE options
fs = 50; 
L = 60; % In seconds
t_span = (1/fs):(1/fs):L;
T = fs*L; % Total number of timesteps
tau = 1/fs; %Time step size (fs = 1/tau Hz)

%% SOlve ODE
% Simulate the Kuramoto model
w = sqrt(.2)*randn([n_osc 1]);
p_cdf = rand([n_osc 1]);
phi0 = (pi/4 + 0.2*tan(pi*(p_cdf-0.5))); %solve cdf eqn for x
% alpha = 10;
% w = trnd(1,n_osc,1) + alpha;
% p_cdf = rand([n_osc 1]);
% phi0 = (alpha + tan(pi*(p_cdf-0.5))); %solve cdf eqn for x
[t,theta] = ode45(@(t,y)Kuramoto_model(t,y,n_osc,w,K),t_span,phi0);% repalce w with r

%% Calculate Phase Coherence
r = abs((1/n_osc)*sum(exp(1i*theta),2));
rtot = mean(r);

%% Calculate PC
osc = sin(theta((5*fs+1):end,:));
EEG = zeros(size(osc));
for i = 1:n_osc
    for j = (i-c):(i+c)
        if j<1
            index = j+n_osc;
            EEG(:,i) = EEG(:,i) + osc(:,index);
        elseif j<=n_osc
            index = j;
            EEG(:,i) = EEG(:,i) + osc(:,index);
        elseif j>n_osc
            index = j-n_osc;
            EEG(:,i) = EEG(:,i) + osc(:,index);
        end
    end
end
EEG = EEG./(2*c +1);

[~,~,~,~,pc,xpc,cc1,cc2] = fn_measure(fs,EEG');


 m1 = triu(ones(size(pc)),1);
 v1 = pc(logical(m1));

 m2 = triu(ones(size(xpc)),1);
 v2 = xpc(logical(m2));


 m3 = triu(ones(size(cc1)),1);
 v3 = cc1(logical(m3));

 m4 = triu(ones(size(cc2)),1);
 v4 = cc2(logical(m4));

 partCorr = mean(v1);
 xpartCorr = mean(v2);
 xcorr1 = mean(v3);
 xcorr2 = mean(v4);