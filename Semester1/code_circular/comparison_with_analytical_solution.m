% -------------------------------------------------------------------------
% This code compares the numerical solution of the Alfvén continuum with 
% the analytical solution of the matrix system, in the case of a circular, 
% concentric, and large aspect ratio equilibriuum; and for matrices of rank
% 2. The analytical solution is given by equation (77) of the report.
% -------------------------------------------------------------------------
%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultLineLinewidth', 1);
format long;

%% Parameters

size=2;     % rank of the matrices A and B
a=0.1;      % minor radius
n=1;        % toroidal mode number
R0=1;       % major radius
nr = 1000;  % number of magnetic surfaces considered
q = @(r) 1+2*(r/a).^2; % safety factor

%% Construction of the matrices and solution of the linear system
r=linspace(0,a,nr);

w=zeros(size,nr);
w2=zeros(size,nr);
K1=zeros(size,nr);
K2=zeros(size,nr);
w_analytic=zeros(size,nr);

for i=1:nr  
    epsilon=r(i)/R0; % inverse aspect ratio
    
    % Choice of the considered values of the poloidal number m. 
    m_min=round((-n*q(r(i))-size/2+1/2));
    m_max=m_min+size-1;

    % Construction of matrices M and N
    M = build_M(epsilon, size);
    N = build_N(n,m_min,m_max,q(r(i)),size);
    k1 = N(1,1);
    k2 = N(2,2);
    
    % solution of the matrix form of the Alfvén contiuum equation
    w2(:,i)=eig(N,M);
    
    % Remove negatives values of w2 (non physical solutions)
    for j=1:size 
        if w2(j,i)<0
            w2(j,i)=nan;
        end
    end
    
    w2(:,i)=sort(w2(:,i));
    w(:,i)=sqrt(w2(:,i));
    
    % Analytical solution, see eq. (77)
    w_analytic(1,i)=sqrt((k1+k2-sqrt((k1-k2).^2+16.*k1.*k2.*epsilon.^2))/(2.*(1-4.*epsilon.^2)));
    w_analytic(2,i)=sqrt((k1+k2+sqrt((k1-k2).^2+16.*k1.*k2.*epsilon.^2))/(2.*(1-4.*epsilon.^2)));

end


%% Figure
% plots the absolute error of the numerical solution with respect to the
% analytical one, as well as associated eigenmodes.

save=0; %Set 1 to save the figure, 0 otherwise
namefig=['comparaison_analytique'];
path='C:\Users\manon\Desktop\projet CSE I\figures\circulaire\';

figure
for i=1:size
semilogy(r./a,abs(w(i,:)-w_analytic(i,:)))
hold on
yyaxis right
plot(r./a,w(i,:))
hold on
yyaxis left
end

grid on
xlabel('$r/a$')
ylabel('Absolute error on $\tilde{\omega}$')
legend('Error on $\tilde{\omega_1}$','Error on $\tilde{\omega_2}$','$\tilde{\omega_1}$','$\tilde{\omega_2}$')
yyaxis right
ylabel('$\tilde{\omega}$')

if save
saveas(gcf,[path,namefig],'epsc')
saveas(gcf,[path,namefig,'.fig'])
end
