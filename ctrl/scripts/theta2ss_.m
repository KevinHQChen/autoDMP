function [A, B] = theta2ss_(theta, sign, y0, no, np)
%#codegen
prms = zeros(np, 2*no); % parameter matrix
for i=1:2*no % apply correct sign to each np segment of theta and set it as a column vector of prms
    prms(:,i) = sign( ((i-1)*np + 1):i*np ) .* theta( ((i-1)*np + 1):i*np );
end

A = diag(ones(1,2*no)+prms(1,:));
B = prms(2:end, :)';
