function [A, B] = theta2ss_(theta, sign, y0)
%#codegen
no = size(y0, 1);
np = size(theta, 1)/no;
prms = zeros(np, no);
for i=1:no
    prms(:,i) = sign( ((i-1)*np + 1):i*np ) .* theta( ((i-1)*np + 1):i*np );
end

A = diag(ones(1,no)+prms(1,:));
B = prms(2:end, :)';
