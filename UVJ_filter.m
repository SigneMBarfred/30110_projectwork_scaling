
clear; close all; clc;

%%
data = load('seds_1.txt');
t_data = readtable('seds_1.txt');

% mu meter
lambda = data(:,1);

% color vector
Vcol = jet(12);

figure
hold on
fig_1 =plot(lambda,data(:,2:13));
set(fig_1, {'color'}, num2cell(jet(12),2));
set(gca, 'yScale', 'log')
set(gca, 'xScale', 'log')
set(gca, 'XTick', [0.1 0.2 0.3 0.5 0.8 1.2 1.8 2.7 4]);
set(gca, 'XLim', [0 5]);
xlabel('Wavelength [\mum]')
ylabel('Flux [\muJy]')
title('Flux as a function of wavelength')
legend('G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12',Location='best',Orientation='vertical')
hold off


%%

% Opdeling i filter - Der antages tophat
% Averaging over U, V og J band

% U - Ultraviolet
U_c = 0.3709; U_w = 0.0518/2;
U = table2array(varfun(@(x)((x>=(U_c-U_w)) & (x<=(U_c+U_w))), t_data(:,1)));      % Logical Vector
U_band = data(U,2:13);

% V - Visual
V_c = 0.5487; V_w = 0.0954/2;
V = table2array(varfun(@(x)((x>=(V_c-V_w)) & (x<=(V_c+V_w))), t_data(:,1)));      % Logical Vector
V_band = data(V,2:13);  

% j - ulta deep j
J_c = 1.2525; J_w = 0.1718/2;
J = table2array(varfun(@(x)((x>=(J_c-J_w)) & (x<=(J_c+J_w))), t_data(:,1)));      % Logical Vector
J_band = data(J,2:13);

% Ave SEDs
Uave = sum(U_band)/size(U_band,1);
Vave = sum(V_band)/size(V_band,1);
Jave = sum(J_band)/size(J_band,1);

% Ave wl
UwlAve = sum(lambda(U))/size(lambda(U),1);
VwlAve = sum(lambda(V))/size(lambda(V),1);
JwlAve = sum(lambda(J))/size(lambda(J),1);

L1 = [linspace(UwlAve,UwlAve,12),linspace(VwlAve,VwlAve,12),linspace(JwlAve,JwlAve,12)];


%%

% flux som bølgelængde med average flux for band filter
figure;
hold on
p1 = plot(lambda,data(:,2:13));
set(p1, {'color'}, num2cell(jet(12),2));
p2 = scatter(L1(:,1:12),Uave,'filled',MarkerFaceColor=[0 0 0]);
p3 = scatter(L1(:,13:24),Vave,'filled',MarkerFaceColor=[0 0 1]);
p4 = scatter(L1(:,25:36),Jave,'filled',MarkerFaceColor=[1 0 0]);
set(gca, 'yScale', 'log')
set(gca, 'xScale', 'log')
set(gca, 'XTick', [0.1 0.2 0.3 0.5 0.8 1.2 1.8 2.7 4]);
set(gca, 'XLim', [0 5]);
xlabel('Wavelength [\mum]')
ylabel('Flux [\muJy]')
title('Flux as a function of wavelength')
hold off

leg1=legend(p1,'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','Location','best');
ah1=axes('position',get(gca,'position'),'visible','off');
leg2=legend(ah1,[p2,p3,p4],'U band', 'V band','J band','Location','northeastoutside');set(leg2,'FontSize',9);
% legendTitle(leg1, 'Broadband type' );
% legendTitle(leg2, 'Galaxy name' );

% saveas(f2,'Flux_wl_aveUVJ','fig');

%%

% ave mag for each band type
% First ave wavelength then convert to mag

m_U = 23.9 - 2.5*log10(Uave);
m_V = 23.9 - 2.5*log10(Vave);
m_J = 23.9 - 2.5*log10(Jave);

% convert entire data to mag

mag = 23.9 - 2.5*   log10(data(:,2:13));

%%

% mag som bølgelængde med average mag for band filter
figure;
hold on
p5 = plot(lambda,mag);
set(p5, {'color'}, num2cell(jet(12),2));
p6 = scatter(L1(:,1:12),m_U,'filled',MarkerFaceColor=[0 0 0]);
p7 = scatter(L1(:,13:24),m_V,'filled',MarkerFaceColor=[0 0 1]);
p8 = scatter(L1(:,25:36),m_J,'filled',MarkerFaceColor=[1 0 0]);
set(gca, 'xScale', 'log')
set(gca, 'YDir','reverse')
set(gca, 'XTick', [0.1 0.2 0.3 0.5 0.8 1.2 1.8 2.7 4]);
set(gca, 'XLim', [0 5]);
xlabel('Wavelength [\mum]')
ylabel('AB magnitude')
title('AB magnitude as a function of wavelength')
hold off

leg3=legend(p5,'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','Location','best');
ah2=axes('position',get(gca,'position'),'visible','off');
leg4=legend(ah2,[p6,p7,p8],'U band', 'V band','J band','Location','northeastoutside');set(leg4,'FontSize',9);


%%

% Determine galaxy at z=1

% Opdeling i filter - Der antages tophat
% Averaging over R, H and K band

% Define lambda obs at z=1
z=1;
t_data.Var1 = (1+z).*t_data.Var1;
lambdaZ_one = (1+z).*data(:,1);

% R - Red
R_c = 0.6219; R_w = 0.1547/2;
R = table2array(varfun(@(x)((x>=(R_c-R_w)) & (x<=(R_c+R_w))), t_data(:,1)));      % Logical Vector
R_band = data(R,2:13);

% H - H ultra deep
H_c = 1.6466; H_w = 0.2905/2;
H = table2array(varfun(@(x)((x>=(H_c-H_w)) & (x<=(H_c+H_w))), t_data(:,1)));      % Logical Vector
H_band = data(H,2:13);  

% K - K ulta deep
K_c = 2.1557; K_w = 0.3074/2;
K = table2array(varfun(@(x)((x>=(K_c-K_w)) & (x<=(K_c+K_w))), t_data(:,1)));      % Logical Vector
K_band = data(K,2:13);

% Ave SEDs
Rave = sum(R_band)/size(R_band,1);
Have = sum(H_band)/size(H_band,1);
Kave = sum(K_band)/size(K_band,1);

% Ave wl
RwlAve = sum(lambdaZ_one(R))/size(lambdaZ_one(R),1);
HwlAve = sum(lambdaZ_one(H))/size(lambdaZ_one(H),1);
KwlAve = sum(lambdaZ_one(K))/size(lambdaZ_one(K),1);

L2 = [linspace(RwlAve,RwlAve,12),linspace(HwlAve,HwlAve,12),linspace(KwlAve,KwlAve,12)];

%%

% flux som bølgelængde med average flux for band filter ved z = 1

figure;
hold on
% p1 = plot(lambda,data(:,2:13));
% p2 = scatter(L1(:,1:12),Uave,'filled',MarkerFaceColor=[1 0.5 0]);
% p3 = scatter(L1(:,13:24),Vave,'filled',MarkerFaceColor=[0 1 0]);
% p4 = scatter(L1(:,25:36),Jave,'filled',MarkerFaceColor=[0 1 1]);
p9 = plot(t_data.Var1,data(:,2:13));
set(p9, {'color'}, num2cell(jet(12),2));
p10 = scatter(L2(:,1:12),Rave,'filled',MarkerFaceColor=[0 0 0]);
p11 = scatter(L2(:,13:24),Have,'filled',MarkerFaceColor=[0 0 1]);
p12 = scatter(L2(:,25:36),Kave,'filled',MarkerFaceColor=[1 0 0]);
set(gca, 'yScale', 'log')
set(gca, 'xScale', 'log')
set(gca, 'XTick', [0.1 0.2 0.3 0.5 0.8 1.2 1.8 2.7 4]);
set(gca, 'XLim', [0 5]);
xlabel('Wavelength [\mum]')
ylabel('Flux [\muJy]')
title('Flux as a function of wavelength for galaxies at z=1')
hold off

leg6=legend(p9,'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12','Location','best');
ah4=axes('position',get(gca,'position'),'visible','off');
leg7=legend(ah4,[p10,p11,p12],'R band', 'H band','K band','Location','northeastoutside');set(leg7,'FontSize',9);


%%
% VUJ diagram

UV = m_U - m_V;
VJ = m_V - m_J;

%%

%eq1 UV
nUV = UV > 1.3;
%eq2 VJ
nVJ = VJ < 1.6;
%eq3
UVVJ = 0.88 * VJ+0.49;

%%
figure
hold on
plot(sort(VJ),sort(UVVJ))
plot([0 0.58 1],[1.3 1.3 1.3])
plot([1.6,1.6],[0,2.5])
hold off


%%

% Skæringspunkt

xVJline = [1.6 1.6];
yVJline = [0 2.5];

xUVline = [0 0.58 1];
yUVline = [1.3 1.3 1.3];

xUVVJ = VJ;
yUVVJ = UVVJ;

[x1,y1] = polyxpoly(xUVline,yUVline, xUVVJ,yUVVJ);
[x2,y2] = polyxpoly(xVJline,yVJline, xUVVJ,yUVVJ);

%%

% Hældning

a=(y2(1)-y1(1))/(x2(1)-x1(1));
b=y2(1)-a*x2(1);

oll = a*1.5+b;

%%

figure
hold on
plot(VJ,UVVJ)
scatter(VJ,UVVJ)
scatter(x1,y1)
scatter(x2,y2)
grid on
hold off

%%
% Define QSGs boundary
% 
xv = transpose([0,0.9205,1.6,1.6,0,0]);
yv = transpose([1.3,1.3,1.8980,2.5,2.5,1.3]);

P = [xv,yv];
pgon = polyshape(P,'SolidBoundaryOrientation','ccw' );
plot(pgon)

TFin = isinterior(pgon,VJ,UV);

numel(VJ(TFin))
numel(VJ(~TFin))


%%

% plot VUJ diagram

Mcol =(jet(12));

figure
plot(xv,yv)
axis equal

hold on
set(gca, 'XLim', [0 2.5]);
set(gca, 'YLim', [0 2.5]);
ylabel('U-V [AB_{mag}]')
xlabel('V-J [AB_{mag}]')
title('VUJ diagram')
legend('G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12',Location='best',Orientation='vertical')
% fig_2 = plot(VJ(TFin),UV(TFin),'*'); % points inside
% set(fig_2, {'color'}, Mcol(TFin,:));
% fig_3 = plot(VJ(~TFin),UV(~TFin),'+'); % points outside
% set(fig_1, {'color'}, Mcol(~TFin,:));
scatter(VJ(1),UV(1),"MarkerFaceColor",Mcol(1,:));
scatter(VJ(2),UV(2),"MarkerFaceColor",Mcol(2,:));
scatter(VJ(3),UV(3),"MarkerFaceColor",Mcol(3,:));
scatter(VJ(4),UV(4),"MarkerFaceColor",Mcol(4,:));
scatter(VJ(5),UV(5),"MarkerFaceColor",Mcol(5,:));
scatter(VJ(6),UV(6),"MarkerFaceColor",Mcol(6,:));
scatter(VJ(7),UV(7),"MarkerFaceColor",Mcol(7,:));
scatter(VJ(8),UV(8),"MarkerFaceColor",Mcol(8,:));
scatter(VJ(9),UV(9),"MarkerFaceColor",Mcol(9,:));
scatter(VJ(10),UV(10),"MarkerFaceColor",Mcol(10,:));
scatter(VJ(11),UV(11),"MarkerFaceColor",Mcol(11,:));
scatter(VJ(12),UV(12),"MarkerFaceColor",Mcol(12,:));
hold off


