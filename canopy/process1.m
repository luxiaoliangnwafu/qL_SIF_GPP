clear;
clc;
close all;

load  dataset_canopy.mat

%% Figure 1a DOY Tair
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.16 .02],[0.12 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 12]);
axes(ha(1));
x_input = GPP_temp.DoY;
y_input = GPP_temp.Tair;

Hour_group_dataset0 = [GPP_temp table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'DateTime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [GPP_temp table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'DateTime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')
% xlim([min(x_input)-days(5) max(x_input)+days(5)])
% set(gca,'XTick',[min(x_input)+days(8): calmonths(1): max(x_input)],'xcolor','k')
% xtickformat('MM/yyyy')
%ylim([0 1500])
xlabel(['Doy of year (Doy, 2022)'],'Color','k')
ylabel(['Air temperature (^oC)'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',19);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

%print('Figure1a_YL_DOY_Tair','-dtiff','-r300')

%% Figure1c DOY Rred+Rnir
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.16 .02],[0.12 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 12]);
axes(ha(1));

x_input = GPP_temp.DoY;
y_input = mean([SIF_dataset1.NIR_v_O2B SIF_dataset2.NIR_v_O2B SIF_dataset3.NIR_v_O2B...
    SIF_dataset4.NIR_v_O2B SIF_dataset5.NIR_v_O2B SIF_dataset6.NIR_v_O2B],2,'omitnan');
Hour_group_dataset0 = [SIF_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;
Hour_group_dataset1 = [SIF_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

yyaxis  left
g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[255 177 177],...
    'MarkerFaceColor',1/255.*[255 177 177]);

hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
ylim([0 0.82])
set(gca,'YTick',[0: 0.2 :1],'ycolor','r')
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

xtickformat
xlabel(['Doy of year (Doy, 2022)'],'Color','k')
ylabel(['{\itR}_6_8_0'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
hold off

yyaxis right

x_input = GPP_temp.DoY;
y_input = mean([SIF_dataset1.NIR_v_O2A SIF_dataset2.NIR_v_O2A SIF_dataset3.NIR_v_O2A...
    SIF_dataset4.NIR_v_O2A SIF_dataset5.NIR_v_O2A SIF_dataset6.NIR_v_O2A],2,'omitnan');
Hour_group_dataset0 = [SIF_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;
Hour_group_dataset1 = [SIF_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g3 = plot(x_input,y_input,'*',...
    'MarkerSize',4,...
    'MarkerEdgeColor',1/255.*[255 177 255],...
    'MarkerFaceColor',1/255.*[255 177 255]);
hold on
g4 =plot(x_input2,y_input2,'s',...
    'MarkerSize',7,...
    'MarkerEdgeColor','m',...
    'MarkerFaceColor','m');
hold off
ylim([0 0.82])
set(gca,'YTick',[0: 0.2 :1],'ycolor','m')
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

ylabel(['{\itR}_7_5_5'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

legend([g1,g2,g3,g4], '{\itR}_6_8_0 Half hourly','{\itR}_6_8_0 Daily','{\itR}_7_5_5 Half hourly','{\itR}_7_5_5 Daily','Location','north','NumColumns',2,'FontSize',20)
legend('boxoff')

%print('Figure1c_YL_DOY_Rred_Rnir','-dtiff','-r300')


%% Figure1d DOY NDVI
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.16 .02],[0.12 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 12]);
axes(ha(1));

x_input = GPP_temp.DoY;
y_input = mean([SIF_dataset1.NDVI SIF_dataset2.NDVI SIF_dataset3.NDVI...
    SIF_dataset4.NDVI SIF_dataset5.NDVI SIF_dataset6.NDVI],2,'omitnan');
Hour_group_dataset0 = [SIF_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;
Hour_group_dataset1 = [SIF_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

xlabel(['Doy of year (Doy, 2022)'],'Color','k')
ylabel(['NDVI'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',19);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
%print('Figure1d_YL_DOY_NDVI','-dtiff','-r300')


%% Figure1e DOY NIRv
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.16 .02],[0.12 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 12]);

axes(ha(1));

x_input = GPP_temp.DoY;
y_input =mean([SIF_dataset1.NDVI.*SIF_dataset1.NIR_v_O2A	 SIF_dataset2.NDVI.*SIF_dataset2.NIR_v_O2A SIF_dataset3.NDVI.*SIF_dataset3.NIR_v_O2A...
    SIF_dataset4.NDVI.*SIF_dataset4.NIR_v_O2A	 SIF_dataset5.NDVI.*SIF_dataset5.NIR_v_O2A SIF_dataset6.NDVI.*SIF_dataset6.NIR_v_O2A],2,'omitnan');
Hour_group_dataset0 = [SIF_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;
Hour_group_dataset1 = [SIF_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;


g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

ylim([0 1])
xlabel(['Doy of year (Doy, 2022)'],'Color','k')
ylabel(['NIR_v'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',19);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
%print('Figure1e_YL_DOY_NIRv','-dtiff','-r300')


%% Figure1f DOY fAPAR
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.16 .02],[0.12 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 12]);

axes(ha(1));

x_input = GPP_temp.DoY;
y_input = mean([SIF_dataset1.fPAR SIF_dataset2.fPAR SIF_dataset3.fPAR...
    SIF_dataset4.fPAR SIF_dataset5.fPAR SIF_dataset6.fPAR],2,'omitnan');
Hour_group_dataset0 = [SIF_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;
Hour_group_dataset1 = [SIF_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;


g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

xlabel(['Doy of year (Doy, 2022)'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',19);
legend('boxoff')
%set(gca,'Xticklabel',[])
ylabel(['{\itf}_A_P_A_R '])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

%print('Figure1f_YL_DOY_fPAR','-dtiff','-r300')




%% Figure1f DOY fesc
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.16 .02],[0.12 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 12]);
axes(ha(1));
x_input = GPP_temp.DoY;
y_input = 0.9.*mean( ...
    [SIF_dataset1.NDVI.*SIF_dataset1.NIR_v_O2A./SIF_dataset1.fPAR ...
    SIF_dataset2.NDVI.*SIF_dataset2.NIR_v_O2A./SIF_dataset2.fPAR ...
    SIF_dataset3.NDVI.*SIF_dataset3.NIR_v_O2A./SIF_dataset3.fPAR ...
    SIF_dataset4.NDVI.*SIF_dataset4.NIR_v_O2A./SIF_dataset4.fPAR ...
    SIF_dataset5.NDVI.*SIF_dataset5.NIR_v_O2A./SIF_dataset5.fPAR ...
    SIF_dataset6.NDVI.*SIF_dataset6.NIR_v_O2A./SIF_dataset6.fPAR] ...
    ,2,'omitnan');

Hour_group_dataset0 = [SIF_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;
Hour_group_dataset1 = [SIF_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

ylim([0 0.5])
xlabel(['Doy of year (Doy, 2022)'],'Color','k')
ylabel(['{\itf}_e_s_c_\__P_-_C'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',19);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

%print('Figure1f_YL_DOY_fesc','-dtiff','-r300')



%% Figure1g DOY SIF760toc
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.16 .02],[0.12 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 12]);

axes(ha(1));

x_input = GPP_temp.DoY;
y_input = 10.*mean([SIF_dataset1.SFM_O2A SIF_dataset2.SFM_O2A SIF_dataset3.SFM_O2A...
    SIF_dataset4.SFM_O2A SIF_dataset5.SFM_O2A SIF_dataset6.SFM_O2A],2,'omitnan');
Hour_group_dataset0 = [SIF_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;
Hour_group_dataset1 = [SIF_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')


ylim([0 2.6])
xlabel(['Doy of year (Doy, 2022)'],'Color','k')
ylabel(['SIF_T_O_C_\__7_6_0  (mW m^-^2 nm^-^1 sr^-^1)'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',19);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

%print('Figure1f_YL_DOY_SIF760','-dtiff','-r300')

%% Figure1h DOY GPPec
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.16 .02],[0.12 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 12]);
axes(ha(1));

x_input = GPP_temp.DoY;
y_input = GPP_temp.GPP_DT_uStar;

Hour_group_dataset0 = [GPP_temp table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'DateTime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [GPP_temp table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'DateTime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([50 160])
set(gca,'XTick',[50: 20 : 160],'xcolor','k')

ylim([0 max(y_input)])
xlabel(['Doy of year (Doy, 2022)'],'Color','k')
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',19);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',19);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

%print('Figure1h_YL_DOY_GPP','-dtiff','-r300')
