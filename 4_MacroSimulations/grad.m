function [y]=grad(x)
for i=1:length(x)
    if i>1 && i<length(x)
        y(i)=(x(i+1)-x(i-1))/2;
    elseif i==1
        y(i)=(x(i+1)-x(end))/2;
    elseif i==length(x)
        y(i)=(x(1)-x(i-1))/2;
    end
end
