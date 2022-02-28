function t = theta_star_mid(r,thmid,kprime,k,delta,deltaprime,d,dprime,q,qbar,B0,R0,npoints)

    Zeta = @(theta,d) theta + asin(d).*sin(theta);
    sd = @(r,d,dprime) (dprime.*r)./(sqrt(1-d.^2)) ;
    R = @(r,theta,delta,d) R0+delta+r.*cos(Zeta(theta,d));
    D = @(r,theta,k,kprime,deltaprime,d,dprime) k*cos(theta).*(deltaprime+cos(Zeta(theta,d))-sd(r,d,dprime).*sin(theta).*sin(Zeta(theta,d)))+(kprime.*r+k).*sin(theta).*(1+asin(d).*cos(theta)).*sin(Zeta(theta,d)) ;
    BdotGradPhi =  @(r,theta,delta,d) B0*R0./(R(r,theta,delta,d).^2);
    BdotGradTheta = @(r,theta,k,kprime,deltaprime,delta,d,dprime,qbar) B0.*k./(qbar.*R(r,theta,delta,d).*D(r,theta,k,kprime,deltaprime,d,dprime));

    integrand = @(theta) BdotGradPhi(r,theta,delta,d)./BdotGradTheta(r,theta,k,kprime,deltaprime,delta,d,dprime,qbar);
    
    t=zeros(1,npoints);
    
    t(1)=(1/q)*midpoint_composite_quadrature(integrand,0,thmid(1),1);
    sum = t(1);
    for i=2:npoints
        sum = sum + (1/q)*midpoint_composite_quadrature(integrand,thmid(i-1),thmid(i),1) ;
        t(i) = sum ; 
    end
    disp(t)
return