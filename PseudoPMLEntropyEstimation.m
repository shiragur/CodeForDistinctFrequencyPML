function [EstimatedEntropy] = PseudoPMLEntropyEstimation(sample,threshold,N)
% This function implements the algorithm of [CSS19] (PseudoPML paper) using the algorithm 
%to compute approximate PML distribution as a sub routine.

succ=1;
[~,n]=size(sample);
distinct=unique(sample);
histogram=histc(sample,distinct); % Set of frequency of elements
vfreq=unique(histogram);
hoh=histc(histogram,vfreq);    % Discretizing the observed frequencies

vfreqbt = vfreq(vfreq(:)<=threshold);
hohbt = hoh(vfreq(:)<=threshold);
k = sum(size(hohbt))-1;
n2 = sum(histogram(histogram(:)<=threshold));
rmin = 1/(10^4*n);
rmax = (2*threshold)/n;

vfreqat = vfreq(vfreq(:)>threshold);
hohat = hoh(vfreq(:)>threshold);
n1 = n-n2;

fracofemp=n1/n

emp_est = 0;
pml_est = 0;

% fracvalue=0;
% fracvalue=n1/n;
if n1>0
emp_prob = vfreqat./n;
emp_prob_entr = (-1)*emp_prob.*log(emp_prob);
m = sum(size(hohat))-1;
bias = m/(2*n);
emp_est = emp_prob_entr*hohat' + bias;
end

if n2>0
pml_entr = ApproximatePMLACSS21(vfreqbt,hohbt,k,rmin,rmax,threshold);
pml_est = (n2/n)*pml_entr - (n2/n)*log(n2/n) ;
end
EstimatedEntropy = pml_est+emp_est;
end