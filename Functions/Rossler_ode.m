function [xcorr1,xcorr2,partCorr,xpartCorr] = Rossler_ode(n_osc,conf)
%% Set the initial Rossler model parameters
fs = 60;
N=6000; % number of samples used
tspan=linspace(0,fs,N); % 50 second and sample rate is 60 Hz
Initial=[1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0]; % initial condition for solving Rossler oscillator equation

c=random('norm',3,1,6,1); % additional noise which is selected from random normal distibution
a=find(c<0);
if ~isempty(a)
c=random('norm',3,1,6,1); 
end

%% Solve ODE
[t,y] = ode45(@(t,y) Rossler_model(t,y,conf,c),tspan,Initial); % solving rossler model
x=y(:,[1,4,7,10,13,16]);  % extracting X segments

%% computing IP
tmp=zeros(N,n_osc);
for i=1:6
sig=x(:,i);
z = hilbert(sig); % Analytic signal
tmp(:,i) = (unwrap(angle(z)));
end
tmp = tmp(3001:end,:);% discard the first 50s

%%  Calculate CC,PC and PCC
[~,~,~,~,pc,xpc,cc1,cc2] = fn_measure(fs,tmp');

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

end

