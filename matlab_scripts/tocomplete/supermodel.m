function [model, c_opt] = supermodel(input, data)

c = randn(1,7);  % random init params
% the model
model = @(c, x) c(7) + c(1)*x.^1 + c(2)*x.^2 + c(3)*x.^3 + c(4)*x.^4 + c(5)*x.^5 + c(6)*x.^6;
fun = @(c) sum((model(c,input)' - mean(data,2)).^2);  %  function to minimize: SSE
options = optimset('Display', 'off', 'LargeScale', 'off');
% optimization procedure
c_opt = fminunc(fun, c, options);