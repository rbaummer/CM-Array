%Constant Modulus Algorithm
%w = initial weight vector
%mu = step size
%x = array input
%C = constraint vector (direction)
function [w, w_c, B, err] = LCCMA(w, mu, x, C)

%run length
[~,R] = size(x);
%array size
N = length(C);

%generate the blocking matrix for a generalized sidelobe canceller using
%the constraint C
%projection matrix onto constraint subspace C
Pc = C*(C'*C)^-1*C';
Pc_orth = diag(ones(N,1)) - Pc;

%find the orthonormalization of Pc_orth
[Q,~] = qr(Pc_orth);
%Blocking matrix B is the first N-1 columns of Q
B = Q(:,1:N-1);

%weights for constraint vector
w_c = C*(C'*C)^-1;

for i = 1:R
    x_reduced = B'*x(:,i);
    %calculate current array output
    y = w_c'*x(:,i) - w(:,i)'*x_reduced;
    %calculate error 
    err(i) = (abs(y)^2-1)*conj(y);
    %update weight vector
    w(:,i+1) = w(:,i) + mu*x_reduced*err(i);
end