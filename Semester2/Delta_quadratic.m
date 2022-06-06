% This function computes the value of the Shafranov shift, for the
% quadratic approximation of the q-profile

function Delta = Delta_quadratic(Deltaprimequad,r,Delta0)
    Delta = zeros(1,length(r));
    for i = 1:length(r)
        Delta(1,i) = trapz(r(1:i),Deltaprimequad(1:i)',1)+Delta0;
    end
end