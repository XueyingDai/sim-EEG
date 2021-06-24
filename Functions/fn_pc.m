function [partCorr] = fn_pc(data)
% Calculate partial correlation for the data.

% INPUTS:
%     fs = sampling frequency
%     data = EEG signal of all channels  (Rows = channels, Column = Time points)
% OUTPUTS:
%     partCorr = partial correlation calculated using elements in the
%     precision matrix

nElec = size(data,1);
nTime = size(data,2);

% Ensure that the data has zero mean and unit variance
% data = data - repmat(mean(data,2),1,nTime); % subtract mean
% data = data./(repmat(std(data,0,2),1,nTime)); % divide by standard dev

%Calculate covariance matrix by definition
[coVar,pValue] = corr(data','Type','Pearson');

% If the coVar matrix is positive definite and invertible, calculate
% the percision matrix which is the inversiton of coVar matrix
precision = eye(nElec,'double');
condition = ~(rcond(coVar) > 1e-14);
if condition == 0 % 0 = positive definite and invertible
    precision = inv(coVar);
end

% Calculate partial correlation
partCorr = zeros(nElec);
for i = 1:nElec
    for j = i+1:nElec
        partCorr(i,j) = -precision(j,i)/(sqrt(precision(i,i)*precision(j,j)));
%         partCorr(j,i) = partCorr(i,j);
    end
end

end

