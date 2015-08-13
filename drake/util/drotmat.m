function dMdtheta = drotmat(theta)
% derivative of the 2D rotation matrix (with respect to theta)

c=cos(theta); 
s=sin(theta);
dMdtheta = [-s,-c; c,-s];
