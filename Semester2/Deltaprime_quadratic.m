function Deltaprime = Deltaprime_quadratic(r,R0,F,q,pprime)
    % This function computes the value of Deltaprime for the quadratic
    % approximation of q (see equation (35) of the report)
    Deltaprime = zeros(1,length(r));
    
    for i = 2:length(r)
        Deltaprime(1,i) = (r(i)/R0)*(betap(i,r,R0,F,q,pprime)+locind(i,r,R0,F,q,pprime)/2);
    end


end

function b = betap(i,r,R0,F,q,pprime)
% local poloidal beta
    mu0 =4*pi*1e-7;
    prefactor = -2*(R0^3*q(i))/(F(i)^1*r(i)^4);
    int = trapz(r(1:i),r(1:i).^3.*pprime(1:i)*mu0,1);
    b = prefactor.*int;
end

function l = locind(i,r,R0,F,q,pprime)
% local internal inductance
    prefactor = -2*(q(i)^2)/(F(i)^2*r(i)^4);
    int = trapz(r(1:i),r(1:i).^3.*F(1:i).^2./q(1:i).^2,1);
    l = prefactor.*int;
end