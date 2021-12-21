function w = eigenmodes_shifted(size,delta_p,a,n,R0,B0,nr,q)
% Computes the Alfvén eigenmodes by solving the matrix form of the Alfvén
% continuum equation, in the case of a circular and shifted equilibriuum. 
% Inputs : 
%       size : rank of the matrices considered, int
%       a :  minor radius, double
%       delta_p : derivative of the shift, double
%       n : toroidal mode number, int
%       R0 : major radius, double
%       B0 : magnetic field, double
%       nr : number of magnetic surfaces considered, int
%       q : safety factor, handle function
% Outputs: 
%       w : Alfvén eigenmodes, matrix of size (size,nr)

F = R0*B0;
r=linspace(0,a,nr);
w=zeros(size,nr);
w2=zeros(size,nr);

for i=1:nr 
    epsilon=r(i)/R0; % inverse aspect ratio
    
    % Choice of the considered values of the poloidal number m. 
    m_min=round((-n*q(r(i))-size/2+1/2));
    m_max=m_min+size-1;
    
    % Construction of matrices M and N
    [M,N] = build_MN_shifted(epsilon,delta_p,n,m_min, m_max, q(r(i)),F,size);
    
    % solution of the matrix form of the Alfvén contiuum equation,
    % (w^2M-N)\xi = 0.
    w2(:,i)=eig(N,M); 
    
    % Remove negatives values of w2 (non physical solutions)
    for j=1:size  
        if w2(j,i)<0
            w2(j,i)=nan;
        end
    end
    
    w2(:,i)=sort(w2(:,i)); % sort the obtained eigenmodes
    w(:,i)=sqrt(w2(:,i));
end

end