%orthogonalized Constant Modulus Algorithm
%w = initial weight vector
%R_inv = initial inverse correlation matrix
%mu = step size
%alpha = forgetting factor
%x = array input
function [w, err, R_inv] = oCMA(w, R_inv, mu, alpha, x)

%run length
[~,R] = size(x);

for i = 1:R
    %calculate current array output
    y(i) = w(:,i)'*x(:,i);
    %calculate error
    err(i) = y(i) - y(i)*abs(y(i))^2;
    %Calculate inverse correlation matrix
    R_inv = R_inv/(1-alpha) - 1/(1-alpha)*(alpha*R_inv*(x(:,i)*x(:,i)')*R_inv)/(1-alpha + alpha*x(:,i)'*R_inv*x(:,i));
    %update weight vector
    w(:,i+1) = w(:,i) + mu*R_inv*x(:,i)*conj(err(i));
end