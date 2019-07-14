function res=ParticuleAbsEllipso(xo,yo,beta,R1,R2,X,Y)
masque_particule = R2.^2*(1-(X-xo).^2./R1.^2-(Y-yo).^2./R2.^2)>0;
res = exp(-2*beta.*sqrt(R2.^2*(1-(X-xo).^2./R1.^2-(Y-yo).^2./R2.^2)).*masque_particule);
%imagesc(res)