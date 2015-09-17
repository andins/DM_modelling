function plot_2Dtimeseries(x, y, time)

% Andrea Insabato 2014
% the function is expecting a 2xN or Nx2 matrix X and a vector time with N
% elements representing the evolution of time

if size(x,2)~=1
    if size(x,1)==1
        x=x'; % make x a column vector 
    else
        error('plot_2Dtimeseries: x should be a vector of length N')
    end
end
if size(y,2)~=1
    if size(y,1)==1
        y=y'; % make y a column vector 
    else
        error('plot_2Dtimeseries: y should be a vector of length N')
    end
end
if size(time,2)~=1
    if size(time,1)==1
        time=time'; % make time a column vector 
    else
        error('plot_2Dtimeseries: time should be a vector of length N')
    end
end
if length(x)~=length(y) || length(time)~=length(y)
    error('plot_2Dtimeseries: x, y and time should be the same length')
end

z = [zeros(size(x)) zeros(size(y))];
        surface([x x],...
                [y y],...
                z,...
                [time time],...
                'facecol','no','edgecol','interp','linew',2)
