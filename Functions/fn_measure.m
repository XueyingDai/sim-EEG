% Function:
% Calculate (1)1st order Partial Correlation, 
%           (2)1st order cross Partial Correlation,
%           (3)cross Correlation including zero lags,
%           (4)cross Correlation excluding zero lags.
% Yields results of z and original value for each measurement

%INPUTS:
%     maxLag = maximum number of lags for partial cross correlation between two vectors(in
%     time points)
%     data = EEG signal (Rows = channels, Column = Time points)

%     
% OUTPUTS:
%     Z = nElec x nElec, partial cross correlation value after z-transfer
%     lagTime = lag time (in time points) at max cross-correlation

% NOTE: For a typical cross-correlation calculation, data1 and data1 will
% be exactly the same.  This function calculates cross correlation between
% all channel pairs. The matrices data1 and data1 will be different for the
% permutation resampling calculation, where data1 represents
% time-shifted data.


% Adapt from crossPartCorrFn5
% Originally from Derek
% Edited by Xueying on Apr/22/2021
% Lopouratory


function [z_pc,z_xpc,z_cc1,z_cc2,pc,xpc,cc1,cc2] = fn_measure(fs,data)

%identiiify useful values
nElec = size(data,1); 
nTime = size(data,2);

data = data - repmat(mean(data,2),1,nTime); % subtract mean
data = data./(repmat(std(data,0,2),1,nTime)); % divide by standard dev

lag = .2*fs; %lag time is 200 ms;
lagVec = -lag:lag;
lagInd = length(lagVec);

pc = zeros(nElec,nElec);% 1st order partial correlation
xpc = zeros(nElec,nElec);% 1st order partial cross correlation
cc1 = zeros(nElec,nElec);% include 0 lag
cc2 = zeros(nElec,nElec);% exclude 0 lag

z_pc = zeros(nElec,nElec);% 1st order partial correlation
z_xpc = zeros(nElec,nElec);
z_cc1 = zeros(nElec,nElec);% include 0 lag
z_cc2 = zeros(nElec,nElec);


% Calculate partial correlations: time shifted and without time shifted
for i=1:nElec
    for j=(i+1):nElec
        elec1 = data(i,:);
        elec2 = data(j,:);
        
        data_rest = data;
        data_rest([i j],:) = [];   
        
        %% Partial correlation
        % collect all unpartial and 1st order partial correlation, and find the minimum
        temp = zeros(nElec-1,1);
        temp(1) = corr(elec1',elec2','Type','Pearson');
        
        for k = 1:nElec-2     
            ts = [elec1;elec2;data_rest(k,:)];
            pp = fn_pc(ts);
            temp(k+1) = pp(1,2);% partial correlation between elec1 and elec2
        end        
        temp = abs(temp);
        pc(i,j) = min(temp);
        
        %% cross Partial correlation
        temp = zeros(nElec-2,1);
        xpp = zeros(2*lag+1,1);        
%         A = xcorr(elec1',elec2','normalized',lag);
%         temp(1) = max(abs(A));
        
        for k = 1:nElec-2
            for tlag = 1:lagInd
                ll = lagVec(tlag);
                data_rest_shifted = circshift(data_rest,ll,2);    
                ts = [elec1;elec2;data_rest_shifted(k,:)];
                pp = fn_pc(ts);
                xpp(tlag) = pp(1,2);% partial correlation between elec1 and elec2
            end
%             temp(k+1) = max(abs(xpp));
            temp(k) = max(abs(xpp));
        end  
        
        temp = abs(temp);
        xpc(i,j) = min(temp);
        
        %% Fisher transformation for pc
        z_pc(i,j) = 0.5*log((1+pc(i,j))/(1-pc(i,j)));
        z_pc(j,i) = z_pc(i,j);

        z_xpc(i,j) = 0.5*log((1+xpc(i,j))/(1-xpc(i,j)));
        z_xpc(j,i) = z_xpc(i,j);       
    end
end

%% cross Correlation (include 0 lag)
B= xcorr(data',lag,'normalized');
[C,ind] = max(abs(B));
cc1 = reshape(C,[nElec,nElec]);
cc1 = triu(cc1,1);
ind2 = reshape(ind,[nElec, nElec]);

%% Cross Correlation (exclude 0 lag)
lagVec = -lag:lag;
maxLag = lagVec(ind2);
realLag = abs(sign(maxLag));
cc2 = cc1.*realLag;
cc2 = triu(cc2,1);

%% Fisher transformation for cc
z_cc1 = 0.5.*log((1+cc1)./(1-cc1));
z_cc2 = 0.5*log((1+cc2)./(1-cc2));

end

