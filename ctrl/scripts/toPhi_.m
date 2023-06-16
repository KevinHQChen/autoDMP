function phi_ = toPhi_(y, u, sign)
%#codegen
no = size(y,1);
ni = no;
np = ni+1;
phi_ = zeros(np*no, no);

regs = [y'; repmat(u, 1, no)]; % normalized regressor matrix
for i=1:no % concatenate each row vector of regs into blkdiag form and apply sign
    phi_( ((i-1)*np + 1):i*np, i ) = sign( ((i-1)*np + 1):i*np ).*regs(:,i);
end
