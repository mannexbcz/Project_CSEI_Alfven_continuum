function int = derive(vals, r)
% Compute the derivative of the values vals, using Finite differences
    int = zeros(length(r),1);
    int(2:end-1,1)= (vals(3:end)-vals(1:end-2))./(r(3:end)-r(1:end-2));
    int(1,1)= (-3*vals(1)+4*vals(2)-vals(3))/(r(3)-r(1));
    int(end,1)=(vals(end-2)-4*vals(end-1)+3*vals(end))/(r(end)-r(end-2));

return