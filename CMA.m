%Constant Modulus Algorithm
%w = initial weight vector
%mu = step size
%x = array input
function [w, err] = CMA(w, mu, x)

%run length
[~,R] = size(x);

for i = 1:R
    %calculate current array output
    y(i) = w(:,i)'*x(:,i);
    %calculate error 
    err(i) = y(i) - y(i)*abs(y(i))^2;
    %update weight vector
    w(:,i+1) = w(:,i) + mu*x(:,i)*conj(err(i));
end