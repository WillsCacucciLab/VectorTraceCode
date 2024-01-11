function [ma,mr] = circ_mean2ndOrder( alpha, r )
% Calculate the mean of a group of mean angles.
%
%       [ma,mr] = circ_mean2ndOrder( alpha, r );
%
% Where alpha is the list of mean angles, and r is the corresponding list of
% rayleigh vectors from the underlying samples.


% Remove NaNs
nanInd = isnan(alpha) | isnan(r) ; 
alpha = alpha( ~nanInd );
r     = r( ~nanInd );


% Make sure angles are in radians, lie in range 0-2pi.
alpha(alpha<0) = alpha(alpha<0) + (2*pi);


% Calculate Xbar and Ybar (Zar eq. 26.29 & 26.30)
XBar = mean( r.*cos(alpha) );
YBar = mean( r.*sin(alpha) );

% And then get the mean angle and RV
mr   = sqrt( XBar^2 + YBar^2 );
ma   = atan2( YBar, XBar );
if ma<0;  ma = ma + (2*pi);  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%converts degrees to radians
function alpha = ang2rad(alpha)
alpha = alpha * pi /180;
