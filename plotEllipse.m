function plotEllipse(M,C)
    N = 50;
    
    th = linspace(-pi,pi,N);
    Basis = [1 0 0; 0 1 0]; %xy
    E = (Basis*inv(M)*Basis');
    
    ell = E^(1/2)*[cos(th); sin(th)];
    
    fill(C(1) + ell(1,:),C(2) + ell(2,:),[0.8 0.8 0.8],'FaceAlpha',0.3);
end