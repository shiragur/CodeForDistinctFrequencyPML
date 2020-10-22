%% This a test file, which primarily invokes function PseudoPMLEntropyEstimation


%Domain size
N=10^5;

%Entropy of Zipf(-0.6) distribution
s1=0.5;
Hns1 = sum(1./((1:N).^s1));
TrueEntropy1=log(Hns1)+(s1/Hns1)*sum(log((1:N)).*((1:N).^(-s1)));


%Entropy of Zipf(-1) distribution
s2=1;
Hns2 = sum(1./((1:N).^s2));
TrueEntropy2=log(Hns2)+(s2/Hns2)*sum(log((1:N)).*((1:N).^(-s2)));


%Entropy of uniform distribution
TrueEntropy3=log(N);

%Mix of two uniforms
% N1=N/10;
% N2=N-N1;
% TrueEntropy3=0.5*log(2*N1)+0.5*log(2*N2);
%TrueEntropy4=-[(1/2)*log(1/(2*N1))+(1/2)*log(1/(2*N2))];


%Sample sizes
ds=9;
K=[10^2,5*10^2,10^3,5*10^3,10^4,5*10^4,10^5,5*10^5,10^6];
zipf06entropy=zeros(ds,1);
zipf1entropy=zeros(ds,1);
uniformentropy=zeros(ds,1);

% timetorun=zeros(ds,2);
% 
% zipf06entropypml=zeros(ds,1);
% zipf1entropypml=zeros(ds,1);
% uniformentropypml=zeros(ds,1);
% 
% zipf06entropymle=zeros(ds,1);
% zipf1entropymle=zeros(ds,1);
% uniformentropymle=zeros(ds,1);

zipf06frac=zeros(ds,1);
zipf1frac=zeros(ds,1);
uniformfrac=zeros(ds,1);

trial=50;

for j=1:ds
    zipf06avg=0;
    zipf1avg=0;
    uniformavg=0;
%     succtraill1=0;
%     succtraill2=0;
%     succtraill3=0;
%     succtraill1pml=0;
%     succtraill2pml=0;
%     succtraill3pml=0;

    
%     zipf06avgpml=0;
%     zipf1avgpml=0;
%     uniformavgpml=0;
%     
%     zipf06avgmle=0;
%     zipf1avgmle=0;
%     uniformavgmle=0;
    
    zipf06fracavg=0;
    zipf1fracavg=0;
    uniformfracavg=0;

for i=1:trial
    fprintf('n=%i, trial=%i',K(j),i)
% Generate a sample of size 10,000 from the uniform distribution of support 100,000
clearvars -except N1 N2 uniformfracavg zipf06avgmle zipf06entropymle zipf1avgmle zipf1entropymle uniformavgmle uniformentropymle timetorun estimateEntropy1pml estimateEntropy2pml estimateEntropy3pml zipf06avgpml zipf1avgpml uniformavgpml zipf06entropypml zipf1entropypml uniformentropypml zipf06frac zipf1frac uniformfrac threshold aa para trial ds i j N K s1 s2 zipf06itersupp zipf1itersupp uniformitersupp zipf06avg zipf1avg uniformavg TrueEntropy1 TrueEntropy2 TrueEntropy3 zipf06entropy zipf1entropy uniformentropy zipf06support zipf1support uniformsupport stdzipf06support stdzipf1support stduniformsupport zipf06fracavg zipf1fracavg uniformfracavg

%succtraill1 succtraill2 succtraill3 succtraill1pml succtraill2pml succtraill3pml
aa=0;

threshold=18;

tic;
samp1=zipf_rand(N,s1,K(j));
estimateEntropy1=PseudoPMLEntropyEstimation(samp1,threshold,N);
%succtraill1=succtraill1+succ;
zipf06avg=zipf06avg+(estimateEntropy1-TrueEntropy1)^2;
zipf06fracavg=zipf06fracavg+aa;
%toc;
% estimateEntropy1pml=PseudoPMLEntropyEstimation(samp1,K(j)+100,N);
% %succtraill1pml=succtraill1pml+succ;
% zipf06avgpml=zipf06avgpml+(estimateEntropy1pml-TrueEntropy1)^2;
% %aa=0;
% 
% hist_vec=int_hist(samp1);
% estimateEntropy1mle=entropyOfDistribution(hist_vec);
% zipf06avgmle=zipf06avgmle+(estimateEntropy1mle-TrueEntropy1)^2;


% tic;
samp2=zipf_rand(N,s2,K(j));

% tic;
estimateEntropy2=PseudoPMLEntropyEstimation(samp2,threshold,N);
% timetorun(j,1)=timetorun(j,1)+toc;
%succtraill2=succtraill2+succ;
zipf1avg=zipf1avg+(estimateEntropy2-TrueEntropy2)^2;
zipf1fracavg=zipf1fracavg+aa;
% toc;
% tic;
% [estimateEntropy2pml,~,succ]=PseudoPMLEntropyEstimation(samp2,K(j)+100,N);
% %[estimateEntropy2pml,~]=estEntroPMLapproximate(samp2,N);
% timetorun(j,2)=timetorun(j,2)+toc;
% zipf1avgpml=zipf1avgpml+(estimateEntropy2pml-TrueEntropy2)^2;
% %succtraill2pml=succtraill2pml+succ;
% 
% hist_vec=int_hist(samp2);
% estimateEntropy2mle=entropyOfDistribution(hist_vec);
% zipf1avgmle=zipf1avgmle+(estimateEntropy2mle-TrueEntropy2)^2;



%aa=0;
% tic;
samp3=randi(N,1,K(j));
%--------
% index=randi([0 1],1,K(j));
% samp3=randi([1 N1],1,K(j)).*index;
% samp3=samp3+randi([N1+1 N],1,K(j)).*(1-index);
%--------
% fprintf('hi1')
% pause(1)
estimateEntropy3=PseudoPMLEntropyEstimation(samp3,threshold,N);
%succtraill3=succtraill3+succ;
uniformavg=uniformavg+(estimateEntropy3-TrueEntropy3)^2;
uniformfracavg=uniformfracavg+aa; 
toc;
% fprintf('hi2')
% pause(1)
% [estimateEntropy3pml,~,succ]=PseudoPMLEntropyEstimation(samp3,K(j)+100,N);
% %[estimateEntropy3pml,~]=estEntroPMLapproximate(samp3,N);
% %succtraill3pml=succtraill3pml+succ;
% uniformavgpml=uniformavgpml+(estimateEntropy3pml-TrueEntropy3)^2;
% 
% hist_vec=int_hist(samp3);
% [m,~]=size(hist_vec(hist_vec(:)>=1));
% estimateEntropy3mle=entropyOfDistribution(hist_vec)+m/(2*K(j));
% uniformavgmle=uniformavgmle+(estimateEntropy3mle-TrueEntropy3)^2;


end


zipf06entropy(j,1)=sqrt(zipf06avg/trial);
% zipf06entropypml(j,1)=sqrt(zipf06avgpml/trial);
% zipf06entropymle(j,1)=sqrt(zipf06avgmle/trial);

zipf1entropy(j,1)=sqrt(zipf1avg/trial);
% zipf1entropypml(j,1)=sqrt(zipf1avgpml/trial);
% zipf1entropymle(j,1)=sqrt(zipf1avgmle/trial);

uniformentropy(j,1)=sqrt(uniformavg/trial);
% uniformentropypml(j,1)=sqrt(uniformavgpml/trial);
% uniformentropymle(j,1)=sqrt(uniformavgmle/trial);

fprintf('n=%i, UniformRMSE=%.4f, Zipf06RMSE=%.4f, Zipf1RMSE=%.4f',K(j),uniformentropy(j,1),zipf06entropy(j,1),zipf1entropy(j,1))

zipf06frac(j,1)=zipf06fracavg/trial;
zipf1frac(j,1)=zipf1fracavg/trial;
uniformfrac(j,1)=uniformfracavg/trial;

% timetorun(j,1)=timetorun(j,1)/trial;
% timetorun(j,2)=timetorun(j,2)/trial;


% zipf06entropy(j,1)=sqrt(zipf06avg/succtraill1);
% zipf06entropypml(j,1)=sqrt(zipf06avgpml/succtraill1pml);
% 
% zipf1entropy(j,1)=sqrt(zipf1avg/succtraill2);
% zipf1entropypml(j,1)=sqrt(zipf1avgpml/succtraill2pml);
% 
% uniformentropy(j,1)=sqrt(uniformavg/succtraill3);
% uniformentropypml(j,1)=sqrt(uniformavgpml/succtraill3pml);
% 
% zipf06frac(j,1)=zipf06fracavg/succtraill1;
% zipf1frac(j,1)=zipf1fracavg/succtraill2;
% uniformfrac(j,1)=uniformfracavg/succtraill3;
% 
% timetorun(j,1)=timetorun(j,1)/succtraill2;
% timetorun(j,2)=timetorun(j,2)/succtraill2pml;
end

save('/Users/shiragur/Documents/MATLAB/JuneFolder/ACSS_PositiveUniform.mat')
