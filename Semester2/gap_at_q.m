function [gap,relativegap] = gap_at_q(valq,q,w1,w2)
% Returns the value of the gap between modes w1 and w2, at q=valq
% Inputs:
%   valq : double
%   q : vector of doubles, q-profile
%   w1,w2: vector of doubles
index = index_q_r(q,valq);
gap = abs(w2(index)-w1(index));
relativegap = gap/abs(w2(index)+w1(index));
end

