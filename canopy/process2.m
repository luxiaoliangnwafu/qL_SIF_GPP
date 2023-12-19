clear;
clc;
close all;
load  dataset_canopy.mat
% elimate zeros in SIF and GPP;
index  = find(GPP_temp.PPFD <= 5) ; 
SIF_dataset1(index,:)=[];
SIF_dataset2(index,:)=[];
SIF_dataset3(index,:)=[];
SIF_dataset4(index,:)=[];
SIF_dataset5(index,:)=[];
SIF_dataset6(index,:)=[];
GPP_temp(index,:)=[];


GPP_MLR_PAR_6set = []; 
GPP_MLR_SIF_6set = []; 
qL_PAR_6set = []; 
qL_SIF_6set = []; 
SIF_PSII_TOT_PPFD_6set = []; 
SIF_TOC_6set = [];
SIF_TOT_6set = [];


for i =1:6

    clear SIF_dataset
    eval(['SIF_dataset=','SIF_dataset',num2str(i),';']);

SIF_dataset = table2timetable(SIF_dataset);
SIF_dataset(SIF_dataset.SFM_O2A<=0,2:end)=table(nan);

PAR = GPP_temp.PPFD;
Tair = GPP_temp.Tair;
rH = GPP_temp.rH/100;
VPD = (1-rH).*(0.6107*exp(17.27.*Tair./(Tair+237.3)));

t = GPP_temp.DateTime;
A = GPP_temp.Ca_temp;
[B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(48),'SamplePoints',t);
GPP_temp.Ca_temp = B;
clear A B t

t = SIF_dataset.Datetime;
A = SIF_dataset.EVI2;
[B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(48),'SamplePoints',t);
SIF_dataset.EVI2 = B;
clear A B t

Ca = interp1(GPP_temp.DateTime, GPP_temp.Ca_temp, SIF_dataset.Datetime,'linear');
fesc_LP = 0.9;
fesc_CP =  fesc_LP.*SIF_dataset.NDVI.*SIF_dataset.NIR_v_O2A./SIF_dataset.fPAR;
 Tau_star = 36.9+1.18.*(Tair-25)+0.036.*(Tair-25).^2; % unit is ppm
 Tau_star = 10.^-3.*Tau_star;% unit is hPa
 Ca = 10.^-3.*Ca;% unit is hPa
VPD(VPD<0) = nan; % unit is hPa
VPD =  10.*VPD;
lamda = 900; % unit is mol H2O molC
Cc = (3.*Ca.*Tau_star - Tau_star.*1.6.*VPD./lamda - Ca.*(3.*Tau_star.*1.6.*VPD./lamda).^0.5)./(3.*Tau_star - 1.6.*VPD./lamda);
Cc = 10.^3.*Cc;% unit is ppm
Tau_star = 10.^3.*Tau_star;% unit is ppm
Ca = 10.^3.*Ca;% unit is ppm

load('RC_ALL_NEW.mat')
[~,~,V2] = svd(RC_PSII);% [~,S1,V1] = svd(RC_PSI);V1_N = V1(:,1)./sum(V1(:,1));
V2_N = V2(:,1)./sum(V2(:,1));
 SIF_WL = 760;
RC_PSII_WL = V2_N(:,1)'./V2_N(SIF_WL-640,1);

clear V2 V2_N SIF_WL RC_PSI RC_PSII SIF_WL

% total PSII SIF
fPSII_wheat = 0.6024; % need check

SIF_PSII_TOT_WL = 10.*SIF_dataset.SFM_O2A .*fPSII_wheat ./ fesc_CP .*RC_PSII_WL; % SFM_O2A  is uw cm-2 nm-1 sr-1 time 10 go to mw m-2 nm-1 sr-1
SIF_PSII_TOT_PPFD = PPFD_ChlF(640:850, SIF_PSII_TOT_WL./10); % PPFD_ChlF is used to mw m-2 nm-1 sr-1 so that it need to divide 10.

SIF_TOC = 10.*SIF_dataset.SFM_O2A;
SIF_TOC(SIF_PSII_TOT_PPFD<0)=nan;
SIF_TOC(SIF_PSII_TOT_PPFD==0)=nan;

SIF_TOT = 10.*SIF_dataset.SFM_O2A./ fesc_CP;
SIF_TOT(SIF_PSII_TOT_PPFD<0)=nan;
SIF_TOT(SIF_PSII_TOT_PPFD==0)=nan;


SIF_PSII_TOT_PPFD(SIF_PSII_TOT_PPFD<0)=nan;
SIF_PSII_TOT_PPFD(SIF_PSII_TOT_PPFD==0)=nan;

%qL's parameter
x = Tair+273.15;
load qL_4parameters.mat 
m = b(1).*(b(2).*exp((b(3).*(x-b(4))./(x.*8.314.*b(4)))))./(b(2)-b(3).*(1-exp((b(2).*(x-b(4))./(x.*8.314.*b(4))))));

SIF_dataset_daily = retime(SIF_dataset,'daily',@(x)mean(x,'omitnan'));% Use 'hourly' timestep with 'mean' method to get the hourly average

GPP_temp_daily = retime(GPP_temp,'daily',@(x)mean(x,'omitnan'));% Use 'hourly' timestep with 'mean' method to get the hourly average
Tair_ave = GPP_temp_daily.Tair;
NDVI_temp= interp1(SIF_dataset_daily.Datetime, SIF_dataset_daily.NDVI, SIF_dataset.Datetime,'linear');

bigleaf_index_min = 0;
bigleaf_index_max = 1;
bigleaf_index = bigleaf_index_min+(NDVI_temp-min(NDVI_temp)).*(bigleaf_index_max-bigleaf_index_min)./(max(NDVI_temp)-min(NDVI_temp));
qL_SIF = m./((SIF_PSII_TOT_PPFD./bigleaf_index).^(1./m)+m);
KDF = 9;
P0 = 0.86; % Pfündel 2021 0.862(SW)-0.883(LW)
GPP_MLR_SIF=(Cc - Tau_star) ./ (4.5.*Cc +10.5.*Tau_star).*(1+KDF).*(P0./(1-P0)).*qL_SIF.*SIF_PSII_TOT_PPFD;


% QL -BASED ON PAR from FigureSX1 qL_PAR 
aqL = 0.7132;% 
bqL = 8.4674e-04; % 
PAR(PAR<0) = 0;
qL_PAR = aqL.*exp(-bqL.*PAR);
GPP_MLR_PAR=(Cc - Tau_star) ./ (4.5.*Cc +10.5.*Tau_star).*(1+KDF).*(P0./(1-P0)).*qL_PAR.*SIF_PSII_TOT_PPFD;


GPP_MLR_PAR_6set = [GPP_MLR_PAR_6set GPP_MLR_PAR]; 
GPP_MLR_SIF_6set = [GPP_MLR_SIF_6set GPP_MLR_SIF]; 
qL_PAR_6set = [qL_PAR_6set qL_PAR]; 
qL_SIF_6set = [qL_SIF_6set qL_SIF]; 
SIF_PSII_TOT_PPFD_6set = [SIF_PSII_TOT_PPFD_6set SIF_PSII_TOT_PPFD]; 

SIF_TOC_6set = [SIF_TOC_6set SIF_TOC];
SIF_TOT_6set = [SIF_TOT_6set SIF_TOT];

end


GPP_MLR_PAR_mean= mean(GPP_MLR_PAR_6set,2,'omitnan');
GPP_MLR_SIF_mean= mean(GPP_MLR_SIF_6set,2,'omitnan');
qL_PAR_mean= mean(qL_PAR_6set,2,'omitnan');
qL_SIF_mean= mean(qL_SIF_6set,2,'omitnan');
SIF_PSII_TOT_PPFD_mean= mean(SIF_PSII_TOT_PPFD_6set,2,'omitnan');
SIF_TOC_mean= mean(SIF_TOC_6set,2,'omitnan');
SIF_TOT_mean= mean(SIF_TOT_6set,2,'omitnan');
GPP_EC = GPP_temp.GPP_DT_uStar;

LAI_min = 0.5;
LAI_max = 4.9;%LAI_obs = 4.9;
EVI2_temp= interp1(SIF_dataset_daily.Datetime, SIF_dataset_daily.EVI2, SIF_dataset.Datetime,'linear');
LAI_all = LAI_min+(EVI2_temp-min(EVI2_temp)).*(LAI_max-LAI_min)./(max(EVI2_temp)-min(EVI2_temp));
clear SIF_dataset_daily

Datetime = SIF_dataset.Datetime;
DoY = GPP_temp.DoY;
Hour = GPP_temp.Hour;

dataset_figure8_10 = table2timetable([table(Datetime)   table(qL_SIF_mean) table(GPP_MLR_SIF_mean)...
    table(qL_PAR_mean) table(GPP_MLR_PAR_mean) table(SIF_PSII_TOT_PPFD_mean) table(GPP_EC)...
    table(LAI_all) table(DoY) table(PAR) table(Hour) table(Tair) table(SIF_TOC_mean) table(SIF_TOT_mean)]);


%% Figure S6 LAI

ax1= figure('visible','on');
set(gcf,'Position',[50 50 500 450]);
set(gca,'position', [0.12 0.15 0.84 0.8]);   % [x0 y0 width height]

x_input = day(dataset_figure8_10.Datetime,'dayofyear') + hour(dataset_figure8_10.Datetime)/24;
y_input = dataset_figure8_10.LAI_all;   

g2 = plot(x_input,y_input,'-','Color','k', 'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

ylim([0 6])
set(gca,'YTick',[0: 1 : 6],'ycolor','k')
ylabel(['LAI  (m^2 m^-^2)'])

xlim([40 160])
set(gca,'XTick',[40: 20 : 160],'xcolor','k')
xlabel(['Day of year (Doy, 2022)'],'Color','k')

legend(g2,'LAI', 'Location','northwest','NumColumns',1,'FontSize',12)
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',16);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
txt = {'Winter wheat'};
text(0.75*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',16)
%print('FigureS6 LAI','-djpeg','-r300')

%% Figure 8-1 Seaonal_GPPSIF_GPPEC

 dataset_figure8_10(dataset_figure8_10.DoY == 71 | dataset_figure8_10.DoY == 156,:)=[];
 dataset_figure8_10(dataset_figure8_10.DoY == 71 | dataset_figure8_10.DoY == 156,:)=[];
dataset_figure8_10_daily = retime(dataset_figure8_10,'daily','mean');% Use 'hourly' timestep with 'mean' method to get the hourly average
dataset_figure8_10_daily_se = retime(dataset_figure8_10,'daily',@(x)std(x,'omitnan'));

figure('visible','on');
set(gcf,'Position',[50 50 950 1650]);

subplot(3,1,1);
set(gca,'position', [0.15 0.7 0.7 0.28]);   % [x0 y0 width height]

x_input = dataset_figure8_10_daily.DoY ;
y_input = dataset_figure8_10_daily.SIF_PSII_TOT_PPFD_mean;
y_input_se = dataset_figure8_10_daily_se.SIF_PSII_TOT_PPFD_mean;
g1 = plot(x_input,y_input,'-o','Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

ylim([0 20])
set(gca,'YTick',[0: 4 : 20],'ycolor','k')
ylabel(['SIF_T_O_T_\__F_U_L_L_\__P_S_I_I (\mumol m^-^2 s^-^1)'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

legend(g1,'SIF_T_O_T_\__F_U_L_L_\__P_S_I_I  Daily','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')

txt = {'(a)'};
text(0.95*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)


subplot(3,1,2);
set(gca,'position', [0.15 0.38 0.7 0.28]);   % [x0 y0 width height]
x_input = dataset_figure8_10_daily.DoY;
y_input = dataset_figure8_10_daily.qL_SIF_mean;
y_input_se = dataset_figure8_10_daily_se.qL_SIF_mean;
g1 = plot(x_input,y_input,'-s','Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);

xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'ycolor','k')
ylim([0 0.8])
set(gca,'YTick',[0: 0.2 : 0.8],'ycolor','k')
ylabel(['{\it{q}}_{L}_\__S_I_F'])

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

legend(g1,'{\it{q}}_{L}_\__S_I_F Daily','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {'(c)'};
text(0.95*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)
hold off

subplot(3,1,3);
set(gca,'position', [0.15 0.06 0.7 0.28]);   % [x0 y0 width height]
x_input = dataset_figure8_10_daily.DoY;
y_input = dataset_figure8_10_daily.GPP_EC  ;
y_input_se = dataset_figure8_10_daily_se.GPP_EC;
g3 = plot(x_input,y_input,'-d','Color','k', 'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k');

hold on
g2 = area(dataset_figure8_10_daily.DoY,y_input+y_input_se);
g2.FaceColor = 'k';
g2.EdgeColor = [1 1 1];
g2.FaceAlpha = 0.2;

g2 = area(dataset_figure8_10_daily.DoY,y_input-y_input_se);
g2.FaceColor = 'w';
g2.EdgeColor = [1 1 1];
g2.FaceAlpha = 1;

g3 = plot(x_input,y_input,'-d','Color','k', 'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k');

x_input = dataset_figure8_10_daily.DoY;
y_input = dataset_figure8_10_daily.GPP_MLR_SIF_mean;
y_input_se = dataset_figure8_10_daily_se.GPP_MLR_SIF_mean;
g4 = plot(dataset_figure8_10_daily.DoY,y_input,'-s','Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);

ci = [y_input-y_input_se,y_input+y_input_se];

xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')
ylim([0 60])
set(gca,'YTick',[0: 10 : 60],'ycolor','k')
xlabel(['Day of year (Doy, 2022)'],'Color','k')

ylabel(['GPP (\mumol m^-^2 s^-^1)'],'Color','k')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
legend([g3, g4],'GPP_E_C Daily', 'GPP_{\it{q}}_{L}_\__S_I_F Daily','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {'(e)'};
text(0.95*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)
hold off
%print('Figure8_1','-djpeg','-r300')



%% Figure 8-2 Seaonal_PAR_qLPAR_GPPSIF_GPPEC

figure('visible','on');
set(gcf,'Position',[50 50 950 1650]);

subplot(3,1,1);
set(gca,'position', [0.15 0.7 0.7 0.28]);   % [x0 y0 width height]

x_input = dataset_figure8_10_daily.DoY ;
y_input = dataset_figure8_10_daily.PAR   ;
y_input_se = dataset_figure8_10_daily_se.PAR;
g1 = plot(x_input,y_input,'-o','Color','b', 'MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);

xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

ylim([0 1200])
set(gca,'YTick',[0: 300 : 1500],'ycolor','k')
ylabel(['PAR (\mumol m^-^2 s^-^1)'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

legend(g1,'PAR  Daily','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')

txt = {'(b)'};
text(0.95*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)


subplot(3,1,2);
set(gca,'position', [0.15 0.38 0.7 0.28]);   % [x0 y0 width height]
x_input = dataset_figure8_10_daily.DoY;
y_input = dataset_figure8_10_daily.qL_PAR_mean;
y_input_se = dataset_figure8_10_daily_se.qL_PAR_mean;
g1 = plot(x_input,y_input,'-s','Color','b', 'MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);

xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'ycolor','k')
ylim([0 0.8])
set(gca,'YTick',[0: 0.2 : 0.8],'ycolor','k')
ylabel(['{\it{q}}_{L}_\__P_A_R'])

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

legend(g1,'{\it{q}}_{L}_\__P_A_R Daily','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {'(d)'};
text(0.95*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)
hold off


subplot(3,1,3);
set(gca,'position', [0.15 0.06 0.7 0.28]);   % [x0 y0 width height]
x_input = dataset_figure8_10_daily.DoY;
y_input = dataset_figure8_10_daily.GPP_EC  ;
y_input_se = dataset_figure8_10_daily_se.GPP_EC;
g3 = plot(x_input,y_input,'-d','Color','k', 'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on
g2 = area(dataset_figure8_10_daily.DoY,y_input+y_input_se);
g2.FaceColor = 'k';
g2.EdgeColor = [1 1 1];
g2.FaceAlpha = 0.2;

g2 = area(dataset_figure8_10_daily.DoY,y_input-y_input_se);
g2.FaceColor = 'w';
g2.EdgeColor = [1 1 1];
g2.FaceAlpha = 1;

g3 = plot(x_input,y_input,'-d','Color','k', 'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k');

x_input = dataset_figure8_10_daily.DoY;
y_input = dataset_figure8_10_daily.GPP_MLR_PAR_mean;
y_input_se = dataset_figure8_10_daily_se.GPP_MLR_PAR_mean;
g4 = plot(dataset_figure8_10_daily.DoY,y_input,'-v','Color','b', 'MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);

ci = [y_input-y_input_se,y_input+y_input_se];

xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')
ylim([0 60])
set(gca,'YTick',[0: 10 : 60],'ycolor','k')
xlabel(['Day of year (Doy, 2022)'],'Color','k')

ylabel(['GPP (\mumol m^-^2 s^-^1)'],'Color','k')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
legend([g3, g4],'GPP_E_C Daily', 'GPP_{\it{q}}_{L}_\__P_A_R Daily','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {'(f)'};
text(0.95*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)
hold off
%print('Figure8-2','-djpeg','-r300')

%% Figure9-1 Diurnal_SIF_QL_GPPSIF_GPPEC_halfhourly

dataset_figure8_10(dataset_figure8_10.Hour<7,:)=[];
dataset_figure8_10(dataset_figure8_10.Hour>19,:)=[];

figure('visible','on');
set(gcf,'Position',[50 50 950 1650]);

subplot(3,1,1);
set(gca,'position', [0.15 0.7 0.7 0.28]);   % [x0 y0 width height]

[C,~,~] = unique(dataset_figure8_10.Hour);
x_input = C;
clear y_input y_input_se y_input2 y_input_se2

for i  = 1: size(C,1)
    y_input(i) = mean(dataset_figure8_10.PAR(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se(i) = std(dataset_figure8_10.PAR(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input2(i) = mean(dataset_figure8_10.SIF_PSII_TOT_PPFD_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se2(i) = std(dataset_figure8_10.SIF_PSII_TOT_PPFD_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
end

g1 =  errorbar(x_input,y_input2,y_input_se2,'-o', 'Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
xlim([min(x_input)-1 max(x_input)+1])
ylim([0 30])
set(gca,'YTick',[0: 5 : 30],'ycolor','k')
ylabel(['SIF_T_O_T_\__F_U_L_L_\__P_S_I_I (\mumol m^-^2 s^-^1)'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

legend([g1],'SIF_T_O_T_\__F_U_L_L_\__P_S_I_I','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {['(a)']};
text(0.95*max(xlim),0.96*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)

subplot(3,1,2);
set(gca,'position', [0.15 0.38 0.7 0.28]);   % [x0 y0 width height]
[C,~,~] = unique(dataset_figure8_10.Hour);
x_input = C;
clear y_input y_input_se y_input2 y_input_se2

for i  = 1: size(C,1)
    y_input(i) = mean(dataset_figure8_10.qL_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se(i) = std(dataset_figure8_10.qL_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input2(i) = mean(dataset_figure8_10.qL_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se2(i) = std(dataset_figure8_10.qL_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
end

g1 = errorbar(x_input,y_input,y_input_se,'-s','Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
xlim([min(x_input)-1 max(x_input)+1])
ylim([0 1])
set(gca,'YTick',[0: 0.2 : 1],'ycolor','k')
ylabel(['{\itq}_L'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
legend([g1],'{\it{q}}_{L}_\__S_I_F','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {'(c)'};
text(0.95*max(xlim),0.96*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)

subplot(3,1,3);
set(gca,'position', [0.15 0.06 0.7 0.28]);   % [x0 y0 width height]
[C,~,~] = unique(dataset_figure8_10.Hour);
x_input = C;
clear y_input y_input_se y_input2 y_input_se2 y_input3 y_input_se3
for i  = 1: size(C,1)
    y_input(i) = mean(dataset_figure8_10.GPP_EC(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se(i) = std(dataset_figure8_10.GPP_EC(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input2(i) = mean(dataset_figure8_10.GPP_MLR_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se2(i) = std(dataset_figure8_10.GPP_MLR_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input3(i) = mean(dataset_figure8_10.GPP_MLR_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se3(i) = std(dataset_figure8_10.GPP_MLR_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
end

g3 = errorbar(x_input,y_input,y_input_se,'-d','Color','k', 'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1);
hold on
g4 = errorbar(x_input,y_input2,y_input_se2,'-^','Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
hold off

xlim([min(x_input)-1 max(x_input)+1])
ylim([0 60])
set(gca,'YTick',[0: 10 : 70],'ycolor','k')
ylabel(['GPP (\mumol m^-^2 s^-^1)'],'Color','k')
xlabel(['Hour of day, 2022'],'Color','k')

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
legend([g3,g4],'GPP_E_C', 'GPP_{\it{q}}_{L}_\__S_I_F', 'Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {['(e)']};
text(0.95*max(xlim),0.96*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)
%print('Figure9-1','-djpeg','-r300')

%% Figure9-2 bf Diurnal_SIF_PAR_GPPSIF_GPPEC_halfhourly

dataset_figure8_10(dataset_figure8_10.Hour<7,:)=[];
dataset_figure8_10(dataset_figure8_10.Hour>19,:)=[];

figure('visible','on');
set(gcf,'Position',[50 50 950 1650]);

subplot(3,1,1);
set(gca,'position', [0.15 0.7 0.7 0.28]);   % [x0 y0 width height]

[C,~,~] = unique(dataset_figure8_10.Hour);
x_input = C;
clear y_input y_input_se y_input2 y_input_se2

for i  = 1: size(C,1)
    y_input(i) = mean(dataset_figure8_10.PAR(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se(i) = std(dataset_figure8_10.PAR(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input2(i) = mean(dataset_figure8_10.SIF_PSII_TOT_PPFD_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se2(i) = std(dataset_figure8_10.SIF_PSII_TOT_PPFD_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
end

g1 = errorbar(x_input,y_input,y_input_se,'-o','Color','b', 'MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
xlim([min(x_input)-1 max(x_input)+1])
ylim([0 2000])
set(gca,'YTick',[0: 500 : 2500],'ycolor','k')
ylabel(['PAR (\mumol m^-^2 s^-^1) '])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));

grid on

legend([g1],'PAR','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {['(b)']};
text(0.95*max(xlim),0.96*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)

subplot(3,1,2);
set(gca,'position', [0.15 0.38 0.7 0.28]);   % [x0 y0 width height]
[C,~,~] = unique(dataset_figure8_10.Hour);
x_input = C;
clear y_input y_input_se y_input2 y_input_se2

for i  = 1: size(C,1)
    y_input(i) = mean(dataset_figure8_10.qL_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se(i) = std(dataset_figure8_10.qL_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input2(i) = mean(dataset_figure8_10.qL_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se2(i) = std(dataset_figure8_10.qL_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
end

 g1 =  errorbar(x_input,y_input2,y_input_se2,'-s', 'Color','b', 'MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
xlim([min(x_input)-1 max(x_input)+1])
ylim([0 1])
set(gca,'YTick',[0: 0.2 : 1],'ycolor','k')
ylabel(['{\itq}_L'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
legend([g1],'{\it{q}}_{L}_\__P_A_R','Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {'(d)'};
text(0.95*max(xlim),0.96*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)


subplot(3,1,3);
set(gca,'position', [0.15 0.06 0.7 0.28]);   % [x0 y0 width height]
[C,~,~] = unique(dataset_figure8_10.Hour);
x_input = C;
clear y_input y_input_se y_input2 y_input_se2 y_input3 y_input_se3
for i  = 1: size(C,1)
    y_input(i) = mean(dataset_figure8_10.GPP_EC(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se(i) = std(dataset_figure8_10.GPP_EC(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input2(i) = mean(dataset_figure8_10.GPP_MLR_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se2(i) = std(dataset_figure8_10.GPP_MLR_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input3(i) = mean(dataset_figure8_10.GPP_MLR_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se3(i) = std(dataset_figure8_10.GPP_MLR_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
end

g3 = errorbar(x_input,y_input,y_input_se,'-d','Color','k', 'MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1);
hold on
%g4 = errorbar(x_input,y_input2,y_input_se2,'-^','Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
g4 = errorbar(x_input,y_input3,y_input_se3,'-v','Color','b', 'MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
hold off

xlim([min(x_input)-1 max(x_input)+1])
ylim([0 60])
set(gca,'YTick',[0: 10 : 70],'ycolor','k')
ylabel(['GPP (\mumol m^-^2 s^-^1)'],'Color','k')
xlabel(['Hour of day, 2022'],'Color','k')

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);

set(gca,'LooseInset',get(gca,'TightInset'));
grid on

legend([g3,g4],'GPP_E_C', 'GPP_{\it{q}}_{L}_\__P_A_R', 'Location','northwest','NumColumns',2,'FontSize',18)
legend('boxoff')
txt = {['(f)']};
text(0.95*max(xlim),0.96*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)
%print('Figure9-2','-djpeg','-r300')

%% Figure10 GPPSIF_GPPEC_halfhourly
 dataset_figure8_10(dataset_figure8_10.GPP_EC>39 & dataset_figure8_10.DoY<80,:)=table(nan);

 dataset_figure8_10(dataset_figure8_10.GPP_MLR_SIF_mean>20 & dataset_figure8_10.DoY==148,:)=table(nan);
 dataset_figure8_10(dataset_figure8_10.GPP_MLR_SIF_mean>30 & dataset_figure8_10.DoY==72,:)=table(nan);
 dataset_figure8_10(dataset_figure8_10.GPP_MLR_SIF_mean>49 & dataset_figure8_10.GPP_EC<23,:)=table(nan);


figure('visible','on');
set(gcf,'Position',[50 50 1600 600]);

subplot(1,2,1);
set(gca,'position', [0.07 0.18 0.38 0.80]);   % [x0 y0 width height]

x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean;
y_input_temp =  dataset_figure8_10.GPP_EC;
color_input = dataset_figure8_10.PAR;

xy_temp = [x_input_temp y_input_temp color_input];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input =xy_temp(:,3);
x_input = t1;
y_input = t2;

ydata=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
xdata=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(xdata,ydata,250,color_bar,'.');
hold on
h(1) = plot(nan,nan,'k.','MarkerSize',20);
h1= colorbar('eastoutside');
 colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
h1.Limits = [0 	2000];
h1.Ticks = [0:300:1800];
xlim([0 80])
set(gca,'XTick',[0: 20 :80],'ycolor','k')
ylim([0 80])
set(gca,'YTick',[0: 20 :80],'ycolor','k')
xlabel(['GPP_{\it{q}}_{L}_\__S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [0 : (55-0)/1000 : 55]';

h(2) = plot(xrange, xrange,'r-.','LineWidth',2); % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(3) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
hold off
lgd = legend(h(1:3),'Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',18);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
box on
grid on
txt = {'(a)'};
text(0.9*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)


subplot(1,2,2);
set(gca,'position', [0.57 0.18 0.38 0.80]);   % [x0 y0 width height]

x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean;
y_input_temp =  dataset_figure8_10.GPP_EC;
color_input = dataset_figure8_10.PAR;

xy_temp = [x_input_temp y_input_temp color_input];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input =xy_temp(:,3);
x_input = t1;
y_input = t2;

ydata=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
xdata=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(xdata,ydata,250,color_bar,'.');
hold on
h(1) = plot(nan,nan,'k.','MarkerSize',20);

h1= colorbar('eastoutside');
 colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
h1.Limits = [0 	2000];
h1.Ticks = [0:300:1800];
xlim([0 80])
set(gca,'XTick',[0: 20 :80],'ycolor','k')
ylim([0 80])
set(gca,'YTick',[0: 20 :80],'ycolor','k')
xlabel(['GPP_{\it{q}}_{L}_\__P_A_R (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [0 : (55-0)/1000 : 55]';

h(2) = plot(xrange, xrange,'r-.','LineWidth',2); % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(3) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
hold off
lgd = legend(h(1:3),'Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',18);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);

rrmse = 100.*RMSE2./(max(y_input)-min(y_input));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
box on
grid on
txt = {'(b)'};
text(0.9*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)
%print('Figure10','-djpeg','-r300')


%% Figure10 GPPSIF_GPPEC_halfhourly
 dataset_figure8_10(dataset_figure8_10.GPP_EC>39 & dataset_figure8_10.DoY<80,:)=table(nan);

 dataset_figure8_10(dataset_figure8_10.GPP_MLR_SIF_mean>20 & dataset_figure8_10.DoY==148,:)=table(nan);
 dataset_figure8_10(dataset_figure8_10.GPP_MLR_SIF_mean>30 & dataset_figure8_10.DoY==72,:)=table(nan);
 dataset_figure8_10(dataset_figure8_10.GPP_MLR_SIF_mean>49 & dataset_figure8_10.GPP_EC<23,:)=table(nan);



% Assuming x and y are your input data
% Example data:
x1 = dataset_figure8_10.SIF_TOC_mean;
y = dataset_figure8_10.GPP_EC;
z = dataset_figure8_10.PAR;
x2= 3.14.*dataset_figure8_10.SIF_TOT_mean;
% Combine x and y into a single matrix
data = [x1, y, x2, z];

% Specify the percentage for training (e.g., 70%)
trainingPercentage = 70;

% Calculate the number of data points for training
numTrainingPoints = round(trainingPercentage / 100 * length(data));

% Randomly shuffle the data
rng('default') % For reproducibility

shuffledData = data(randperm(length(data)), :);

% Select 70% of the shuffled data for training
trainingData = shuffledData(1:numTrainingPoints, :);

% Select the remaining 30% for validation
validationData = shuffledData(numTrainingPoints+1:end, :);

% Extract x and y from training and validation data
xTraining = trainingData(:, 1);
yTraining = trainingData(:, 2);

yValidation = validationData(:, 2);

% Perform linear fitting on the training data
mdl1 = fitlm(trainingData(:, 1),  trainingData(:, 2));
mdl2 = fitlm(trainingData(:, 3),  trainingData(:, 2));

% Predict y values using the fitted model on the validation data

yPredictedValidation1 = predict(mdl1,validationData(:, 1));
yPredictedValidation2 = predict(mdl2,validationData(:, 3));



%% Figure10_comment50
figure('visible','on');
set(gcf,'Position',[50 50 1600 600]);

subplot(1,2,1);
set(gca,'position', [0.07 0.18 0.38 0.80]);   % [x0 y0 width height]

x_input_temp = yPredictedValidation1;
y_input_temp =  yValidation;
color_input = validationData(:,4);

xy_temp = [x_input_temp y_input_temp color_input];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input =xy_temp(:,3);
x_input = t1;
y_input = t2;

ydata=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
xdata=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(xdata,ydata,250,color_bar,'.');
hold on
h(1) = plot(nan,nan,'k.','MarkerSize',20);

h1= colorbar('eastoutside');
 colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
h1.Limits = [0 	2000];
h1.Ticks = [0:300:1800];
xlim([0 80])
set(gca,'XTick',[0: 20 :80],'ycolor','k')
ylim([0 80])
set(gca,'YTick',[0: 20 :80],'ycolor','k')
xlabel(['GPP_S_I_F_\__T_O_C_\__L_I_N (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [0 : (55-0)/1000 : 55]';

h(2) = plot(xrange, xrange,'r-.','LineWidth',2); % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(3) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
hold off
lgd = legend(h(1:3),'Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',24)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',28);
box on
grid on
txt = {'(a)'};
text(0.9*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)


subplot(1,2,2);
set(gca,'position', [0.57 0.18 0.38 0.80]);   % [x0 y0 width height]

x_input_temp = yPredictedValidation2;
y_input_temp =  yValidation;
color_input = validationData(:,4);

xy_temp = [x_input_temp y_input_temp color_input];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input =xy_temp(:,3);
x_input = t1;
y_input = t2;

ydata=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
xdata=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(xdata,ydata,250,color_bar,'.');
hold on
h(1) = plot(nan,nan,'k.','MarkerSize',20);

h1= colorbar('eastoutside');
 colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
h1.Limits = [0 	2000];
h1.Ticks = [0:300:1800];
xlim([0 80])
set(gca,'XTick',[0: 20 :80],'ycolor','k')
ylim([0 80])
set(gca,'YTick',[0: 20 :80],'ycolor','k')
xlabel(['GPP_S_I_F_\__T_O_T_\__L_I_N (\mumol m^-^2 s^-^1)'],'Color','k')
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));

xrange = [0 : (55-0)/1000 : 55]';

h(2) = plot(xrange, xrange,'r-.','LineWidth',2); % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(3) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
hold off
lgd = legend(h(1:3),'Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',24)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',28);
box on
grid on
txt = {'(b)'};
text(0.9*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)
%print('Figure10_comment50','-djpeg','-r300')



%% FigureS5 Seaonal_GPPSIF_GPPEC daily

dataset_figure8_10_daily = retime(dataset_figure8_10,'daily','mean');% Use 'hourly' timestep with 'mean' method to get the hourly average
dataset_figure8_10_daily_se = retime(dataset_figure8_10,'daily',@(x)std(x,'omitnan'));

figure('visible','on');
set(gcf,'Position',[50 50 1600 600]);

subplot(1,2,1);
set(gca,'position', [0.07 0.18 0.4 0.80]);   % [x0 y0 width height]
x_input =  dataset_figure8_10_daily.GPP_MLR_SIF_mean;
y_input = dataset_figure8_10_daily.GPP_EC;

clear color_input_temp;
h(1) = plot(x_input,y_input,'^','Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r');
hold on
xlim([0 50])
set(gca,'XTick',[0: 10 :60],'ycolor','k')
ylim([0 50])
set(gca,'YTick',[0: 10 :60],'ycolor','k')
xlabel(['GPP_{\it{q}}_{L}_\__S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
h(2) = plot(xrange, xrange,'r-.','LineWidth',2); % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(3) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
hold off
lgd = legend(h(1:3),'Daily','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',18);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
box on
grid on
txt = {'(a)'};
text(0.9*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)


subplot(1,2,2);
set(gca,'position', [0.56 0.18 0.4 0.80]);   % [x0 y0 width height]

x_input =  dataset_figure8_10_daily.GPP_MLR_PAR_mean;
y_input = dataset_figure8_10_daily.GPP_EC;

clear color_input_temp;
h(1) = plot(x_input,y_input,'v','Color','b', 'MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','b');
hold on
xlim([0 50])
set(gca,'XTick',[0: 10 :60],'ycolor','k')
ylim([0 50])
set(gca,'YTick',[0: 10 :60],'ycolor','k')
xlabel(['GPP_{\it{q}}_{L}_\__P_A_R (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
h(2) = plot(xrange, xrange,'r-.','LineWidth',2); % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(3) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
hold off
lgd = legend(h(1:3),'Daily','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',18);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
box on
grid on
txt = {'(b)'};
text(0.9*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)
%print('FigureS5','-djpeg','-r300')

%% FigureS4 Diurnal GPPSIF_GPPEC_halfhourly

[C,~,~] = unique(dataset_figure8_10.Hour);
x_input = C;
clear y_input y_input_se y_input2 y_input_se2 y_input3 y_input_se3
for i  = 1: size(C,1)
    y_input(i) = mean(dataset_figure8_10.GPP_EC(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se(i) = std(dataset_figure8_10.GPP_EC(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input2(i) = mean(dataset_figure8_10.GPP_MLR_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se2(i) = std(dataset_figure8_10.GPP_MLR_SIF_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
    
    y_input3(i) = mean(dataset_figure8_10.GPP_MLR_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan');
    y_input_se3(i) = std(dataset_figure8_10.GPP_MLR_PAR_mean(dataset_figure8_10.Hour == C(i)),'omitnan')';
end

figure('visible','on');
set(gcf,'Position',[50 50 1600 600]);

subplot(1,2,1);
set(gca,'position', [0.07 0.18 0.4 0.80]);   % [x0 y0 width height]

x_input1 =  y_input2(4:20);
y_input1 = y_input(4:20);

clear color_input_temp;
h(1) = plot(x_input1,y_input1,'^','Color','r', 'MarkerSize',9,'MarkerEdgeColor','r','MarkerFaceColor','r');
hold on
xlim([0 40])
set(gca,'XTick',[0: 10 :60],'ycolor','k')
ylim([0 40])
set(gca,'YTick',[0: 10 :60],'ycolor','k')
xlabel(['GPP_{\it{q}}_{L}_\__S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
model_temp = fitlm(x_input1,y_input1);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input1):(max(x_input1)-min(x_input1))./1000:max(x_input1)]';
h(2) = plot(xrange, xrange,'r-.','LineWidth',2); % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(3) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
hold off
lgd = legend(h(1:3),'Hour of day','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',18);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input1,y_input1),-4);
%rrmse =RMSE2./mean(y_input,'omitnan')
rrmse = 100.*RMSE2./(max(y_input1)-min(y_input1));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.55*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);

box on
grid on
txt = {'(a)'};
text(0.9*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)

subplot(1,2,2);
set(gca,'position', [0.56 0.18 0.4 0.80]);   % [x0 y0 width height]
% 
% axes(ha(1));
x_input1 =  y_input3(4:20);
y_input1 = y_input(4:20);

clear color_input_temp;
h(1) = plot(x_input1,y_input1,'v','Color','b', 'MarkerSize',9,'MarkerEdgeColor','b','MarkerFaceColor','b');
hold on
xlim([0 40])
set(gca,'XTick',[0: 10 :60],'ycolor','k')
ylim([0 40])
set(gca,'YTick',[0: 10 :60],'ycolor','k')
xlabel(['GPP_{\it{q}}_{L}_\__P_A_R (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
model_temp = fitlm(x_input1,y_input1);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input1):(max(x_input1)-min(x_input1))./1000:max(x_input1)]';
h(2) = plot(xrange, xrange,'r-.','LineWidth',2); % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(3) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
hold off
lgd = legend(h(1:3),'Hour of day','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',18);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input1,y_input1),-4);
%rrmse =RMSE2./mean(y_input,'omitnan')
rrmse = 100.*RMSE2./(max(y_input1)-min(y_input1));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope^* = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.55*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);

box on
grid on
txt = {'(b)'};
text(0.9*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)
%print('FigureS4','-djpeg','-r300')

 
%% Table 1
median(dataset_figure8_10.PAR,'omitnan')


R2_LAI = nan(4, 12); % all low mid high_LAI


LAI_index_mid = find(dataset_figure8_10.LAI_all >= 0.5 & dataset_figure8_10.LAI_all <= 4 & dataset_figure8_10.PAR >30);
x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean(LAI_index_mid);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_mid);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,1)  = [R22; RMSE2; rrmse; N];


LAI_index_mid = find(dataset_figure8_10.LAI_all >= 0.5 & dataset_figure8_10.LAI_all <= 4 & dataset_figure8_10.PAR >30 & dataset_figure8_10.PAR < 500);
x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean(LAI_index_mid);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_mid);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,2)  = [R22; RMSE2; rrmse; N];



LAI_index_mid = find(dataset_figure8_10.LAI_all >= 0.5 & dataset_figure8_10.LAI_all <= 4 & dataset_figure8_10.PAR >600);
x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean(LAI_index_mid);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_mid);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,3)  = [R22; RMSE2; rrmse; N];



LAI_index_high = find(dataset_figure8_10.LAI_all > 4 & dataset_figure8_10.PAR >30);
x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean(LAI_index_high);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_high);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,4)  = [R22; RMSE2; rrmse; N];


LAI_index_high = find(dataset_figure8_10.LAI_all > 4 & dataset_figure8_10.PAR >30 & dataset_figure8_10.PAR < 600);
x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean(LAI_index_high);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_high);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,5)  = [R22; RMSE2; rrmse; N];

LAI_index_high = find(dataset_figure8_10.LAI_all > 4 & dataset_figure8_10.PAR >600);
x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean(LAI_index_high);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_high);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,6)  = [R22; RMSE2; rrmse; N];

% PAR

LAI_index_mid = find(dataset_figure8_10.LAI_all >= 0.5 & dataset_figure8_10.LAI_all <= 4 & dataset_figure8_10.PAR >30);
x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean(LAI_index_mid);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_mid);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,7)  = [R22; RMSE2; rrmse; N];


LAI_index_mid = find(dataset_figure8_10.LAI_all >= 0.5 & dataset_figure8_10.LAI_all <= 4 & dataset_figure8_10.PAR >30 & dataset_figure8_10.PAR < 600);
x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean(LAI_index_mid);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_mid);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,8)  = [R22; RMSE2; rrmse; N];


LAI_index_mid = find(dataset_figure8_10.LAI_all >= 0.5 & dataset_figure8_10.LAI_all <= 4 & dataset_figure8_10.PAR >600);
x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean(LAI_index_mid);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_mid);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,9)  = [R22; RMSE2; rrmse; N];

LAI_index_high = find(dataset_figure8_10.LAI_all > 4 & dataset_figure8_10.PAR >30);
x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean(LAI_index_high);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_high);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,10)  = [R22; RMSE2; rrmse; N];


LAI_index_high = find(dataset_figure8_10.LAI_all > 4 & dataset_figure8_10.PAR >30 & dataset_figure8_10.PAR < 600);
x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean(LAI_index_high);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_high);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,11)  = [R22; RMSE2; rrmse; N];

LAI_index_high = find(dataset_figure8_10.LAI_all > 4 & dataset_figure8_10.PAR >600);
x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean(LAI_index_high);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_high);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_LAI(:,12)  = [R22; RMSE2; rrmse; N];


%% Table2

mean(dataset_figure8_10.Tair,'omitnan')

R2_Tair = nan(4,4); % all low mid high_LAI

LAI_index_low = find(dataset_figure8_10.Tair <= 17.4 &  dataset_figure8_10.PAR >30 );
x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean(LAI_index_low);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_low);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_Tair(:,1)  = [R22; RMSE2; rrmse; N];

LAI_index_high = find(dataset_figure8_10.Tair > 17.4 & dataset_figure8_10.PAR >30);
x_input_temp = dataset_figure8_10.GPP_MLR_SIF_mean(LAI_index_high);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_high);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_Tair(:,2)  = [R22; RMSE2; rrmse; N];


LAI_index_low = find(dataset_figure8_10.Tair <= 17.4  & dataset_figure8_10.PAR >30);
x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean(LAI_index_low);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_low);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_Tair(:,3)  = [R22; RMSE2; rrmse; N];

LAI_index_high = find(dataset_figure8_10.Tair > 17.4 & dataset_figure8_10.PAR >30);
x_input_temp = dataset_figure8_10.GPP_MLR_PAR_mean(LAI_index_high);
y_input_temp =  dataset_figure8_10.GPP_EC(LAI_index_high);
xy_temp = [x_input_temp y_input_temp];
xy_temp(isnan(xy_temp(:,1)),:) = [];
xy_temp(isnan(xy_temp(:,2)),:)=[];
x_input_temp = xy_temp(:,1);
y_input_temp =  xy_temp(:,2);
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);
color_input = dataset_figure8_10.PAR;
x_input =t1;
y_input = t2;
model_temp = fitlm(x_input,y_input);
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(y_input)-min(y_input));
N = size(x_input,1);
R2_Tair(:,4)  = [R22; RMSE2; rrmse; N];

% clearvars -except  R2_SAT  R2_CI R2_Tair  R2_LAI R2_LAI_PAR





