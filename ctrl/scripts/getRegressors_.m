function [z, phi] = getRegressors_(y, yPrev, u, sign_, no, ni, np, dt, mdlNum)
%#codegen
z = zeros(no, 1);
phi = zeros(no*np, no);

regs = [y'; repmat(u, 1, no)]; % normalized regressor matrix
for i=1:no % concatenate each col vector of regs into blkdiag form and apply sign
    phi( ((i-1)*np + 1):i*np, i ) = sign_( ((i-1)*np + 1):i*np ).*regs(:,i);
end

z = y - yPrev;
