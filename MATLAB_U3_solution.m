%solution to TSAna_coursework3
%--------------------------------------------------------------------------
%(c) 2017, Georgia Koppe, Dept. Theoretical Neuroscience, CIMH, Heidelberg
%University
%for comments, questions, errors, please contact
%Georgia.Koppe@zi-mannhiem.de

clear all; close all; clc

%--------------------------------------------------------------------------
%Task 1 - univariate ARMA
%--------------------------------------------------------------------------
load('Tut3_file1.mat');


%For this task, use the first time series (i.e. row 1) of DLPFC.
%--------------------------------------------------------------------------
%1.1 Estimate an AR(3) model. Determine which of the preceding 3 time steps
%(i.e., coefficients {a1, a2 a3}) are significantly influencing the current
%time step using the t statistic (compute t-statistics explicitly,
%look up value in tdistribution (Matlab function tcdf)).
%What do the results tell you about the process?

y=DLPFC(1,:);               %@E and C: what you want to predict, y_i in this case one of your cluster time series...

%linear regression on past p=3 time steps
T=length(y); p=3;
Xp=zeros(T-p,p+1);          %@E and C: predictor matrix (here AR model, but you can fill in whatever)
Xp(:,1)=ones(T-p,1);
for i=0:p-1
    Xp(:,p-i+1)=y(i+1:T-p+i)';
end
yT=y(p+1:end)';             %@E and C: here we cut off one value b/c it's an AR model, you don't have to do this!!!
XpXp=(Xp'*Xp)^-1;
a=XpXp*Xp'*yT;              %@E and C: solution to the regression weights 
xpred1=Xp*a;
eps=yT-xpred1;

%t-statistic with df= T-2p-1
sdeps=std(eps,1);
v11=XpXp(2,2); v22=XpXp(3,3); v33=XpXp(4,4);        
a1=a(2); a2=a(3); a3=a(4);
t1=a1/(sdeps*sqrt(v11));        %@E and C: T-value for parameters
t2=a2/(sdeps*sqrt(v22));
t3=a3/(sdeps*sqrt(v33));    
df=T-2*p-1;                 %@E and C: degrees of freedom, should be I think T- number of predictors, but check that in the book
tcdf(t1,df)
tcdf(t2,df)
tcdf(t3,df)
