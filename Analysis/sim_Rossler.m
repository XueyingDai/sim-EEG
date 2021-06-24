% This script evaluates multivariate phase synchrony (MPS) measures on the Rossler oscillators. 
% The script generates Table 2 of the paper.
% Before running this script, Function folder needs to be added to MATLAB paths.
%
% Reference: 
%
% Payam Shahsavari Baboukani, Ghasem Azemi, Boualem Boashash, Paul Colditz, Amir Omidvarnia,
% A novel multivariate phase synchrony measure: Application to multichannel newborn EEG analysis,
% Digital Signal Processing, Volume 84, 2019, Pages 59-68, ISSN 1051-2004,
% https://doi.org/10.1016/j.dsp.2018.08.019.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edited by Xueying, Jun 2021

clc;
clear variables;
% solving equation with ODE45 and initial condition of [1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0]
CH=8; % number of configuration
n_osc = 6;
n_run=500; % number of runs for monte-carlo method

partCorr=zeros(CH,n_run);
xpartCorr=zeros(CH,n_run);

xcorr1=zeros(CH,n_run);
xcorr2=zeros(CH,n_run);

for ii=1:n_run
    for j=1:CH
        conf = j;
         [xcorr1(j,ii),xcorr2(j,ii),partCorr(j,ii),xpartCorr(j,ii)] = Rossler_ode(n_osc,conf);
    end
    disp(['Run: ' num2str(ii) ' was completed!'])
end

% save('rossler_8_model_16value_for_1st_order_xpc.mat')
%%
% calculating mean and std values over 500 runs
s_pc1=zeros(CH,1);
m_pc1=mean(partCorr,2);
for i=1:CH
    s_pc1(i)=std(partCorr(i,:));
end

s_pcN=zeros(CH,1);
m_pcN=mean(xpartCorr,2);
for i=1:CH
    s_pcN(i)=std(xpartCorr(i,:));
end

s_cc=zeros(CH,1);
m_cc=mean(xcorr1,2);
for i=1:CH
    s_cc(i)=std(xcorr1(i,:));
end

subplot(2,2,1)
boxplot(xcorr1')
ylim([0 1])
title('Cross Correlation')

subplot(2,2,2)
boxplot(xcorr2')
ylim([0 1])
title('Cross Correlation (exclude values at zero lag)')

subplot(2,2,4)
boxplot(partCorr')
ylim([0 1])
title('1st-order Partial Correlation')

subplot(2,2,3)
boxplot(xpartCorr')
ylim([0 1])
title('1st-order Partial Cross Correlation')

% s_pc1=zeros(CH,1);
% m_pc1=mean(z_partCorr,2);
% for i=1:CH
%     s_pc1(i)=std(z_partCorr(i,:));
% end
% 
% s_pcN=zeros(CH,1);
% m_pcN=mean(z_xpartCorr,2);
% for i=1:CH
%     s_pcN(i)=std(z_xpartCorr(i,:));
% end
% 
% s_cc=zeros(CH,1);
% m_cc=mean(z_xcorr1,2);
% for i=1:CH
%     s_cc(i)=std(z_xcorr1(i,:));
% end
% 
% subplot(2,2,1)
% boxplot(z_xcorr1')
% % ylim([0 1])
% title('Cross Correlation')
% 
% subplot(2,2,2)
% boxplot(z_xcorr2')
% % ylim([0 1])
% title('Cross Correlation (exclude values at zero lag)')
% 
% subplot(2,2,3)
% boxplot(z_partCorr')
% % ylim([0 1])
% title('1st-order Partial Correlation')
% 
% subplot(2,2,4)
% boxplot(z_xpartCorr')
% % ylim([0 1])
% title('1st-order Partial Cross Correlation')