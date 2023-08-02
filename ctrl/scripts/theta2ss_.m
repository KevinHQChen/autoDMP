function [A, B] = theta2ss_(theta, sign, no, ni, np, ns, dt, mdlNum)
%#codegen
prms = zeros(np, 2*no); % parameter matrix
A = zeros(2*ns, 2*ns);
B = zeros(2*ns, ni);

for i=1:2*no % apply correct sign to each np segment of theta and set it as a column vector of prms
    prms(:,i) = sign( ((i-1)*np + 1):i*np ) .* theta( ((i-1)*np + 1):i*np );
    A(i, i) = [1 - prms(1,i)];
    B(i, :) = prms(2:end, i)';
end
