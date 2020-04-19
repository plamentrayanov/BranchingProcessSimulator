function [Z_mean, Z_lower, Z_upper, Z_median]=confInterval(Z, alpha)
% [Z_MEAN, Z_LOWER, Z_UPPER, Z_MEDIAN]=confInterval(Z, ALPHA) calculates the confidence intervals of the branching 
% process Z with confidence level ALPHA.
%
% Z_MEAN - the expectation of the branching process
% 
% Z_LOWER and Z_UPPER - the lower and upper confidence bounds, resprectively
%
% Z_MEDIAN - the median of the branching process

Z_lower=zeros(1,size(Z,2));
Z_upper=zeros(1,size(Z,2));
for t=1:size(Z,2)
    [F,X]=ecdf(Z(:,t));
    Z_lower(t)=X(find(F>=alpha/2, 1, 'first'))';
    Z_upper(t)=X(find(F>=1-alpha./2, 1, 'first'))';
end
Z_mean=mean(Z,1)';
Z_median=median(Z,1)';
end