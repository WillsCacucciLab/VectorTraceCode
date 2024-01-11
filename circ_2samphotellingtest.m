function [pval, F] = circ_2samphotellingtest(alpha1,alpha2,r1,r2)
%   Parametric Hotelling 2-independent sample test for equal angular means of 2nd order data.
%
%   Details of this test can be found in Zar Section 27.11, where is introduced as follows:
%   'Batschelet (1978, 1981) explained how the Hotelling (1931) procedure can to extended 
%   to [..] the equality of means ..'.
%
%   The test has 2 by N-3 degrees of freedom where N is the sum of phase samples 1 and 2. 
%
%
%   H0: the two populations have equal mean phase angles
%   HA: the two populations have unequal mean phase angles
%   
%   Usage:
%     [pval, F] = circ_htest(alpha1, alpha2, r1 ,r2)
%
%   Input:
%     alpha1:  angles in radians 
%     alpha2:  angles in radians 
%     r1: resultant vector lengths of alpha1, in case the values
%     in alpha1 consist of mean values of distributions themselves
%     r2: resultant vector lengths of alpha2, in case the values
%     in alpha2 consist of mean values of distributions themselves
%
%   Notes:
%     Inputs can be either a row or column vector. r1 and r2 need to be the same 
%     length as respective alpha, and each value in r should have the same index as 
%     their alpha counterpart.
%
%   Output:
%     pval    p-value of the Hotelling paired-sample test. Discard H0 if
%             pval is small.
%     F       the test statistic of the Hottelling test
%
% References:
%   Biostatistical Analysis, J. H. Zar (1999)
%
% TW adaption of 'circ_htest' (for paired samples), by:
% RL van den Brink, 2014
% r.l.van.den.brink@fsw.leidenuniv.nl

%%% check the input
if nargin == 0
    error('not enough input arguments')
elseif nargin == 4
    if length(r1) ~= length(alpha1) || length(r2) ~= length(alpha2)
        error('vector of resultant lengths must be the same length as phase values')
    end    
    %make sure they are row vectors
    if size(r1,2) > size(r1,1); r1 = r1'; end
    if size(r2,2) > size(r2,1); r2 = r2'; end
else
    error('too many input arguments')
end

%%% transform and normalize if needed
%make sure these are all row vectors
if size(alpha1,2) > size(alpha1,1); alpha1 = alpha1'; end
if size(alpha2,2) > size(alpha2,1); alpha2 = alpha2'; end

% TW edit: Remove Nans
nanInd = isnan(alpha1) | isnan(r1) ; 
alpha1 = alpha1( ~nanInd );
r1     = r1( ~nanInd );
nanInd = isnan(alpha2) | isnan(r2) ; 
alpha2 = alpha2( ~nanInd );
r2     = r2( ~nanInd );


% normalize to between 0 and 2pi
%rad if needed (for the calculations below the phases have to be positive
%and in radians)
alpha1(alpha1 < 0) = alpha1(alpha1 < 0) + (2*pi);
alpha2(alpha2 < 0) = alpha2(alpha2 < 0) + (2*pi);

%%% Run the test
k1 = length(alpha1); %number of phase pairs
k2 = length(alpha2); %number of phase pairs

%get normalized rectangular coordinates
X1 = r1.*cos(alpha1);
X2 = r2.*cos(alpha2);
Y1 = r1.*sin(alpha1);
Y2 = r2.*sin(alpha2);

%get (normalized) grand mean rectangular coordinates
Xbar1 = mean(X1);
Xbar2 = mean(X2);
Ybar1 = mean(Y1);
Ybar2 = mean(Y2);

% Get sigma x2, y2, xy, for each sample (i.e. eq. 27.24 - 27.26).
sigmaxsq1 = sum(X1.^2) - (sum(X1)^2 ./ k1);
sigmaysq1 = sum(Y1.^2) - (sum(Y1)^2 ./ k1);
sigmaxy1  = sum(X1.*Y1) - ((sum(X1).*sum(Y1)) / k1);
sigmaxsq2 = sum(X2.^2) - (sum(X2)^2 ./ k2);
sigmaysq2 = sum(Y2.^2) - (sum(Y2)^2 ./ k2);
sigmaxy2  = sum(X2.*Y2) - ((sum(X2).*sum(Y2)) / k2);

% Get combined sigma quantities (ie. eq. 27.31 - 27.33).
sigmaxsqc = sigmaxsq1 + sigmaxsq2;
sigmaysqc = sigmaysq1 + sigmaysq2;
sigmaxyc  = sigmaxy1  + sigmaxy2;

% Get F-statistic, see Biostatistical Analysis, Zar (1999), equation 27.34
const =  (k1+k2-3) / (2*((1/k1) + (1/k2)));
F     =  const  *  (((( Xbar1-Xbar2 )^2*sigmaysqc) - (2*( Xbar1-Xbar2 ) * ( Ybar1-Ybar2 ) * sigmaxyc) + (( Ybar1-Ybar2 )^2 * sigmaxsqc)) / ...
                            ((sigmaxsqc * sigmaysqc) - (sigmaxyc^2)) );
%get p-value                    
pval = 1 - fcdf(F,2,k1+k2-3); 


