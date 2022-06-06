function w = eigenmodes_textbookcase(M,N,size,band,epsilon,epsilon2,epsilon3,a,n,F,nr,q,npoints)
% Computes the Alfvén eigenmodes by solving the matrix form of the Alfvén
% continuum equation, in the textbook case (eq.(78)).
% Inputs : 
%       M,N : equilibrium coefficients, handle functions
%       size : size of the matrices considered, int
%       band : number of band of the matrices, int
%       epsilon, epsilon2, epsilon3 : coefficients, doubles
%       a :  minor radius, double
%       n : toroidal mode number, int
%       F : toroidal flux current, double
%       nr : number of magnetic surfaces considered, int
%       q : safety factor, handle function
%       npoints : number of mesh points for the numerical integration, int
% Outputs: 
%       w : Alfvén eigenmodes, matrix of size (size,nr)

w=zeros(size,nr);
w2=zeros(size,nr);

r=linspace(0,a,nr);

for i=1:nr
    
    m_min=round((-n*q(r(i))-size/2+1/2));
    m_max=m_min+size+1;

    m_min_coeff=0;
    m_max_coeff=(band+1)/2;
    
    % Fourier coefficients
    [coeffs_M] = get_fourier_coeff_pt_milieu(M,epsilon,epsilon2,epsilon3,m_min_coeff,m_max_coeff,npoints);
    [coeffs_N] = get_fourier_coeff_pt_milieu(N,epsilon,epsilon2,epsilon3,m_min_coeff,m_max_coeff,npoints);
    
    % Construction of matrices M and N
    [Mmat,Nmat] = build_MN_bands_Fourier(coeffs_M, coeffs_N, m_min, m_max, n, q(r(i)), F, size, band);

    %solution du système
    w2(:,i)=eig(Nmat,Mmat); 
    for j=1:size %Enlever les valeurs de omega^2 qui sont negatives, donc non physiques
        if w2(j,i)<0
            w2(j,i)=nan;
        end
    end
    w2(:,i)=sort(w2(:,i));
    w(:,i)=sqrt(w2(:,i));
end
end