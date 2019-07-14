function res_tot=ProfileIntensiteEllipso(omega0,R1,R2,I0,xo,yo,beta,size_x,size_y,Nsampling)
x=linspace(-size_x/2,size_x/2,Nsampling); %bornes d'intégration
y=linspace(-size_y/2,size_y/2,Nsampling);
[X,Y]=meshgrid(x,y);
cmap = jet(length(beta)); % colormap pour graphe  
cmap2 = parula(length(beta));
res_mult = []; % Matrice 3D selon ((xo,Itransmis),R)
cell_res_mult = {};
cell_minAbs = {[]}; % Tableau minimums d'absorption
minAbs=[];
minMat=[];
betaTab = [];
minMatTab = [];
fwhmTab = [];

cell_HFWHM = {}; % Tableau FWHM
HFWHM=[];
figure; % on appelle la figure 1
for iiii=1:length(beta)
    disp(beta(iiii));
for iii=1:length(R1)
    R1(iii)
    tic
    for tt=1:length(yo)  
        for ii=1:length(xo) % déplacement sphéroide
            SignalIntegral(ii) = SignalAbsEllipso(omega0,R1(iii),R2(iii),I0,xo(ii),yo(tt),beta(iiii),X,Y); 
        end
        plot(xo,SignalIntegral/max(SignalIntegral),'.-','LineWidth',1,'Color',cmap2(iiii,:)) %plot du passage du sphéroide
        hold on % on retient les plots pour tracer à la fin de la boucle
        res_tot =[xo; SignalIntegral]; % résultat pour un sphéroide de rayon défini
        
        %size(res_tot)
        res_mult = cat(3,res_mult,res_tot); % stockage pour un R dans la matrice 3D selon R
        %size(res_mult)        
    end
    toc   
end
for jjj=1:size(res_mult,3) %on parcourt la matrice 3D selon R
    minAbs = [minAbs min(res_mult(2,:,jjj)./max(res_mult(2,:,jjj)))]; % on cherche le minimum et on le stocke dans un tableau
    FWHMpos = 0.5 * (max(res_mult(2,:,jjj))+min(res_mult(2,:,jjj))); % valeur de I associée à la FWHM
    indexFWHM = find(res_mult(2,:,jjj) < FWHMpos); % recherche des indices correspondants aux valeurs les plus proches
    FWHM1 = xo(indexFWHM(length(indexFWHM))+1);
    FWHM2 = xo(indexFWHM(length(indexFWHM))-1);
    HFWHM = [HFWHM abs(FWHM1+FWHM2)./2]; % calcul de la FWHM et stockage dans un tableau
end
disp(minAbs);
cell_minAbs{iiii} = minAbs; % on cherche le minimum et on le stocke dans un tableau
cell_HFWHM{iiii} = HFWHM; % calcul de la FWHM et stockage dans un tableau
res_mult = [];
minMat(iiii) = min(minAbs(:,:));
minAbs = [];
HFWHM=[];
minAbsCell={};
end
hold off
figure('DefaultAxesFontSize',14)
for iiii=1:length(beta)
plot(R1, cell_minAbs{iiii}, 'o-','LineWidth',1,'Color',cmap(iiii,:)) % on trace la donnée minAbs vs R
hold on
end
%save('minAbsR_inf_2000/minAbs0_001.dat','minAbs','-ascii');
grid on
xlabel('R (µm)')
ylabel('Transmittance (norm.)')
hold off
figure(3);
for iiii=1:length(beta)
plot(R1, 2*cell_HFWHM{iiii}, 'o-','LineWidth',1,'Color',cmap(iiii,:)) % on trace la donnée minAbs vs R
hold on
end
grid on
xlabel('R (µm)')
ylabel('FWHM (µm)')
hold off
%save('HWHMR_inf_2000/HWHM0_001.dat','HFWHM','-ascii');
figure(4);
for iiii=1:length(beta)
    betaTab = [betaTab beta(iiii)];
    fwhmTab = [fwhmTab 2*cell_HFWHM{iiii}];
end
plot(betaTab, fwhmTab, 'o-','LineWidth',1,'Color','red') % on trace la donnée minAbs vs R
hold on
grid on
xlabel('\beta (cm^{-1})')
ylabel('FWHM (�m)')
hold off
%save('HWHMR_inf_2000/HWHM0_001.dat','HFWHM','-ascii');
%disp(minMat);
figure(5);
%miniMat=cell2mat(cell_minAbs{1});
%disp(miniMat);
for iiii=1:length(beta)
%assignin('base',minAbsCell,cell_minAbs);

%if mod(iiii,2)==0
    %betaTab = [betaTab beta(iiii)];
minMatTab = [minMatTab minMat(iiii)];    
%end

end
disp(size(minMatTab));
disp(size(betaTab));
plot(betaTab, minMatTab, 'o-','LineWidth',1,'Color','red') % on trace la donnée minAbs vs R
hold on
grid on
xlabel('\beta (cm^{-1})')
ylabel('Transmittance (norm.)')
hold off
%save('HWHMR_inf_2000/HWHM0_001.dat','HFWHM','-ascii');


