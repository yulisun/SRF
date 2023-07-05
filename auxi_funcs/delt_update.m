function [delt] = delt_update(Q,delt1,delt2,lambda,eta,u,delt_type,fuse_type)
if strcmp(fuse_type,'dot_fixed') == 1
    alfa = - eta * sqrt(sum(delt2.^2,1)) + lambda + max(sqrt(sum(delt2.^2,1)))*eta/8; 
    alfa = alfa/ u;
elseif strcmp(fuse_type,'dot_fixed_Gradient') == 1
    step_size = 0.01;
    dy = sqrt(sum(delt2.^2,1));
    for i = 1:100        
        dx = sqrt(sum(delt1.^2,1));
        Gradient = repmat((- eta * dy + lambda)./(dx+eps),[size(delt1,1),1]).*delt1 + u*(delt1-Q);
        delt1 = delt1 - step_size * Gradient;
    end
    delt = delt1;
    delt_type = 3;
end

switch delt_type
    case 1
        delt = shrinkage(Q,alfa);
    case 21
        Q1 = sqrt(sum(Q.^2,1));
        Q1(Q1==0) = alfa(Q1==0);
        Q2 = (Q1 - alfa) ./ Q1;
        delt = Q * diag((Q1 > alfa) .* Q2);
end



