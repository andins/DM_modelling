% August 23rd 2010
% Steps the filtered white noise forward
function result=stepNoise(Inoise,s,d,sigma)
    
result=Inoise*s+d*sigma*randn(size(Inoise));
    
