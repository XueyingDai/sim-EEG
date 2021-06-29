% This script (1) Simulates Kuramoto oscillators 50 times for each different K and c (i0) values.
%             (2) For each oscillator system, calculates the average level of snychrony with different metrices, 
%                 including Cross Correlation, Cross Correlation (set to 0 if maximum at 0 lag),
%                 1st Order Partial Correlation, and 1st order Partial Cross Correlation.
%             (3) Plot the average snychrony across 50 times as a function of K and c(i0)

% This script is adapted from 
% https://www.mathworks.com/matlabcentral/fileexchange/69807-multivariate-phase-synchrony-analysis-of-eeg?s_tid=srchtitle

% Reference: 
% Payam Shahsavari Baboukani, Ghasem Azemi, Boualem Boashash, Paul Colditz, Amir Omidvarnia,
% A novel multivariate phase synchrony measure: Application to multichannel newborn EEG analysis,
% Digital Signal Processing, Volume 84, 2019, Pages 59-68, ISSN 1051-2004,
% https://doi.org/10.1016/j.dsp.2018.08.019.

% Created by Xueying Dai
% June 2021
% Lopouratory
% University of California, Irvine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Monte Carlo test
n_run = 50;

n_osc = 19; % Total number of oscillators of Kuramoto Model

k = 0:0.5:6; % coupling strength
k_tot = length(k);

c = 0:2:6; % i0 in paper, models volume conduction
c_tot = length(c);

r = zeros(n_run,k_tot,c_tot);% phase coherence
pc = zeros(n_run,k_tot,c_tot);% 1st order partial correlation
xpc = zeros(n_run,k_tot,c_tot);% 1st order partial cross correlation
corr1 = zeros(n_run,k_tot,c_tot);% cross correlation
corr2 = zeros(n_run,k_tot,c_tot);% cross correlation (set to 0 at zero lag)

for cindex = 1:c_tot    
    vc = c(cindex);    
    for kindex = 1:k_tot
        kt = k(kindex);        
        for i = 1:n_run            
            [r(i,kindex,cindex),corr1(i,kindex,cindex),corr2(i,kindex,cindex),pc(i,kindex,cindex),xpc(i,kindex,cindex)] = Kuramoto_ode(n_osc,kt,vc);           
        end        
        disp(['k = ', num2str(kt),' is completed.'])
    end
    disp(['c = ',num2str(vc),' is completed.'])
end
%%
subplot (2,3,1)
rmean = mean(r(:,:,1),1);
rmean2 = (rmean-min(rmean))/(max(rmean)-min(rmean));

plot(k,rmean,'LineWidth',1)
xlim([0 7])
ylim([0 1])
xlabel('Coupling Strength K')
title('Phase Coherence')

subplot(2,3,2)
rmean = mean(corr1(:,:,1),1);
rmean3 = mean(corr1(:,:,2),1);
rmean4 = mean(corr1(:,:,3),1);
rmean5 = mean(corr1(:,:,4),1);

plot(k,rmean,k,rmean3,k,rmean4,k,rmean5,'LineWidth',1)
legend('i = 0','i = 2','i = 4','i = 6')
% plot(k,rmean)
% legend('i = 0')
xlim([0 7])
ylim([0 1])
xlabel('Coupling Strength K')
title('Cross Correlation')

subplot(2,3,3)
rmean = mean(corr2(:,:,1),1);
rmean3 = mean(corr2(:,:,2),1);
rmean4 = mean(corr2(:,:,3),1);
rmean5 = mean(corr2(:,:,4),1);

plot(k,rmean,k,rmean3,k,rmean4,k,rmean5,'LineWidth',1)
legend('i = 0','i = 2','i = 4','i = 6')
% plot(k,rmean)
% legend('i = 0')
xlim([0 7])
ylim([0 1])
xlabel('Coupling Strength K')
title('Cross Correlation (Set to 0 at zero lag)')

subplot (2,3,6)
rmean = mean(pc(:,:,1),1);
rmean3 = mean(pc(:,:,2),1);
rmean4 = mean(pc(:,:,3),1);
rmean5 = mean(pc(:,:,4),1);

plot(k,rmean,k,rmean3,k,rmean4,k,rmean5,'LineWidth',1)
legend('i = 0','i = 2','i = 4','i = 6')
% plot(k,rmean)
% legend('i = 0')
xlim([0 7])
ylim([0 1])
xlabel('Coupling Strength K')
title('1st order Partial Correlation')

subplot (2,3,5)
rmean = mean(xpc(:,:,1),1);
rmean3 = mean(xpc(:,:,2),1);
rmean4 = mean(xpc(:,:,3),1);
rmean5 = mean(xpc(:,:,4),1);

plot(k,rmean,k,rmean3,k,rmean4,k,rmean5,'LineWidth',1)
legend('i = 0','i = 2','i = 4','i = 6')
% plot(k,rmean)
% legend('i = 0')
xlim([0 7])
ylim([0 1])
xlabel('Coupling Strength K')
title('1st order Partial Cross Correlation')