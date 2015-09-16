% August 23rd 2010
% Function for f-I curve

function result=phi(I)

% Parameters for f-I curve, values in Hz
a=270.0;
b=108.0;
c=0.154; 

x=a*I-b;

if(abs(c*x)<1e-6)
    result=.5*x+1.0/c; % This is just for numerical stability
else
    result=x/(1-exp(-c*x));
end;

