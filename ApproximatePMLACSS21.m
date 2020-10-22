function [pml_entr] = ApproximatePMLACSS21(vfreq,histofhist,k,rmin,rmax,threshold)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
rmax = min(rmax,1);
alpha=1/(threshold);
templ=floor(log(rmax/rmin)/log(1+alpha));
vr=[rmin * (1+alpha).^(0:templ),rmax];
l=templ+2; %Actuall l


%Generating the histogram and frequency vector
% distinct=unique(sample);
% histogram=histc(sample,distinct); % Set of frequency of elements
% vfreq=unique(histogram);
% histofhist=histc(histogram,vfreq);    % Discretizing the observed frequencies
% k=sum(size(histofhist))-1; 


C=-log(vr)'*vfreq; %C_ij=-n_j log p_i

%constrainterror=exp(-20)
cvx_begin quiet
cvx_solver mosek
%cvx_precision best
variable X(l,k) nonnegative;
variable NX(l,1) nonnegative;
minimize(qsdpmlobj(X,NX,C,k))
subject to
sum(X,1) == histofhist;
vr*(sum(X,2)+NX) <= 1;
cvx_end

%Rounding Algorithm
FinalX = max([X,NX],zeros(l,k+1));

Rowsum = sum(FinalX,2); %X_i values
Rowsum(1) = Rowsum(1)+ceil(sum(Rowsum))-sum(Rowsum);
for i = 1:l-1
    rem = Rowsum(i)-floor(Rowsum(i));
    Rowsum(i) = floor(Rowsum(i));
    Rowsum(i+1)=Rowsum(i+1)+rem;
end
Rowsum(l) = floor(Rowsum(l));
c = vr*Rowsum;
vr = vr./c;
vr_entr = (-1)*vr.*log(vr);
pml_entr = vr_entr*Rowsum;
end
