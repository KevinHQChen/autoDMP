function [ywt, ywtT, uwt, uwtT] = wtMod_(ytotal, yDest, ywtT, uwtT, dt, no, ni, k_2)
%#codegen
ywt = zeros(1,2*no);
uwt = zeros(1,ni);

% time-scaled sigmoid (around 0->1 in 0->k_2 seconds, k_2 being a time scale factor)
r = 0.2; % 10% to 90% rise time
dwt = dt/k_2;
k_1 = 2.197/(r*k_2);
x0 = 0.5*k_2; % midpoint
sigmoid = @(x) 1/(1 + exp(-k_1*(x-x0)));

for i = 1:2*no
    if yDest(i) ~= 0
        % drive ywt to 1
        if (ywtT(i) <= 1)
            ywtT(i) = ywtT(i) + dwt;
        end
    else
        % drive ywt to 0
        if (ywtT(i) > 0)
            ywtT(i) = ywtT(i) - dwt;
        end
    end
    if ywtT(i) <= 0
        ywt(i) = 0;
    else
        ywt(i) = sigmoid(ywtT(i)*k_2);
    end
end

% for i = 1:ni
%     if yDest(i) ~= 0
%         % drive uwt to 0
%         if (uwtT(i) > 0)
%             uwtT(i) = uwtT(i) - dwt;
%         end
%     else
%         % drive uwt to 1
%         if (uwtT(i) <= 1)
%             uwtT(i) = uwtT(i) + dwt;
%         end
%     end
%     uwt(i) = sigmoid(uwtT(i)*k_2);
% end
