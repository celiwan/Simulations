function SignalIntegral = SignalAbsEllipso(omega0,R1,R2,Io,xo,yo,beta,X,Y)
Faisceau = Io*exp(-2*(X.^2+Y.^2)/(omega0.^2));
bille = ParticuleAbsEllipso(xo,yo,beta,R1,R2,X,Y);
BilleFaisceau =Faisceau.*bille;
dx=X(1,2)-X(1,1);
dy=Y(2,1)-Y(1,1);
%imagesc(BilleFaisceau);
SignalIntegral = sum(sum(BilleFaisceau*dx*dy));
