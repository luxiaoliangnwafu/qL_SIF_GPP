%% Fig4-1 corn qL_PAM vs ChlF_PSII

% Leaf level results
clear;
clc;
close all;

load leaf_dataset.mat

color_temp =jet(7);

ax1= figure('visible','on');
set(gcf,'Position',[50 50 1600 840]);

  h1 = subplot(1,1,1);
  set(gca,'position', [0.06 0.17 0.9 0.75]);   % [x0 y0 width height]

qL_mean_maize       = mean(qL_cal(1:12,1:28),2,'omitnan');
qL_se_maize            = std(qL_cal(1:12,1:28),0,2,'omitnan');
TChlF_mean_maize = mean(TChlF(1:12,1:28),2,'omitnan');
TChlF_se_maize      = std(TChlF(1:12,1:28),0,2,'omitnan');

% x =TChlF_mean_maize;
% y = qL_mean_maize;
% x_std =TChlF_se_maize;
% y_std = qL_se_maize;

x = [mean(TChlF(1:12,1:4),2)
       mean(TChlF(1:12,5:8),2)
       mean(TChlF(1:12,9:12),2)
       mean(TChlF(1:12,13:16),2)
       mean(TChlF(1:12,17:20),2)
       mean(TChlF(1:12,21:24),2)
       mean(TChlF(1:12,25:28),2)];

y = [mean(qL_cal(1:12,1:4),2)
       mean(qL_cal(1:12,5:8),2)
       mean(qL_cal(1:12,9:12),2)
       mean(qL_cal(1:12,13:16),2)
       mean(qL_cal(1:12,17:20),2)
       mean(qL_cal(1:12,21:24),2)
       mean(qL_cal(1:12,25:28),2)];
%plot(TChlF(1:12,1:28),qL_cal(1:12,1:28),'.','MarkerSize',12,'Color',1/255.*[230 230 230],'MarkerFaceColor',1/255.*[230 230 230],'LineWidth',1);
i =1;
h(2) = plot(mean(TChlF(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(1,:));
hold on
i =5;
h(3) = plot(mean(TChlF(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(2,:));
i =9;
h(4) = plot(mean(TChlF(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(3,:));
i =13;
h(5) = plot(mean(TChlF(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(4,:));
i =17;
h(6) = plot(mean(TChlF(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(5,:));
i =21;
h(7) = plot(mean(TChlF(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(6,:));
i =25;
h(8) = plot(mean(TChlF(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(7,:));

% errorbar(x,y,y_std,y_std,x_std,x_std,'d','MarkerSize',14,'Color','k','MarkerFaceColor','k','LineWidth',0.8);
% g(3) = plot(x,y,'d','MarkerSize',14,'Color','k','MarkerFaceColor','k');

[xData, yData] = prepareCurveData(x, y);
ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.2;
%Fit model to data.
[fitresult, gof] = fit(xData, yData, ft, opts );
m_PAR_QL_maize = fitresult.a;
r2_PAR_QL_maize= gof.rsquare;
rmse_PAR_QL_maize= gof.rmse;
rrmse = 100.*rmse_PAR_QL_maize./(max(yData)-min(yData));
xData = (min(xData):(max(xData)-min(xData))/1000:max(xData))';
Yhat = feval(fitresult, xData);
g(4)=plot(xData,Yhat,'Color','k','LineWidth',2);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);

%xlim([0 32])
set(gca,'XTick',0: 5 :35,'ycolor','k')
ylim([0.0 1])
set(gca,'YTick',0: 0.2 :1,'ycolor','k')

g(1) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(2) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(3) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(5) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(6) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(7) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(8) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(9) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');

hold off
xlabel(['SIF_T_O_T_\__F_U_L_L_\__P_S_I_I (\mumol m^-^2 s^-^1)'],'Color','k')% ChlFTOT_FULL_PSII
ylabel('{\itq}_L_\__P_A_M','Color','k')

% all pvalue > 0.01
legend(g(1:9),'','', '','Fit line',...
    '{\it{y}} = m/({\it{x}}^{1/m}+m)',...
    ['{\it{m}} = ',num2str(m_PAR_QL_maize,'%6.2f')],...
    ['{\it{R}}^2 = ',num2str(r2_PAR_QL_maize,'%6.2f')],...
    ['rRMSE = ',num2str(rrmse,'%6.2f'),'%'],...
    'Location','northeast','NumColumns',1,'FontSize',22);
legend('boxoff')

txt = {'(a)'};
text(0.93*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)

set(h1,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(h1,'LooseInset',get(gca,'TightInset'));

grid(h1,'on')
box(h1,'on')

% Copy the axes and plot the second legned
ah1=axes('position',get(gca,'position'),'visible','off');
lh2 = legend(ah1, h(2:8),'15 ^{o}C','20 ^{o}C','25 ^{o}C','30 ^{o}C','35 ^{o}C','40 ^{o}C','45 ^{o}C','Location','southoutside','NumColumns',7,'FontSize',22,'Orientation','horizontal');
legend('boxoff')
lh2.Position = [0.23 0.95 0.7583 0.0328];% [0.25 0.005 0.7583 0.0328];
set(lh2,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
%set(lgd,'Fontname','times new Roman','FontWeight','bold','FontSize',20);


%print('Fig4-1 qL_TChlF','-djpeg','-r300')


%% Fig4-2 Tair-qL-PAR_wheat


ax1= figure('visible','on');
set(gcf,'Position',[50 50 1600 500]);

h1 = subplot(1,3,1);
set(gca,'position', [0.065 0.21 0.26 0.70]);   % [x0 y0 width height]

x_input = [15 20 25 30 35 40 45];
i=3;
j=1; y1 =mean(qL_cal(i+1,j:j+3),2); y1_se = std(qL_cal(i+1,j:j+3));
j=5; y2 =mean(qL_cal(i+1,j:j+3),2); y2_se = std(qL_cal(i+1,j:j+3));
j=9; y3 =mean(qL_cal(i+1,j:j+3),2); y3_se = std(qL_cal(i+1,j:j+3));
j=13; y4 =mean(qL_cal(i+1,j:j+3),2); y4_se = std(qL_cal(i+1,j:j+3));
j=17; y5 =mean(qL_cal(i+1,j:j+3),2); y5_se = std(qL_cal(i+1,j:j+3));
j=21; y6 =mean(qL_cal(i+1,j:j+3),2); y6_se = std(qL_cal(i+1,j:j+3));
j=25; y7 =mean(qL_cal(i+1,j:j+3),2); y7_se = std(qL_cal(i+1,j:j+3));
y_input = [y1 y2 y3 y4 y5 y6 y7];

clear h;
    errorbar(x_input(1),y1,y1_se,y1_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
hold on
h(1) = plot(x_input(1),y1,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(1,:));

    errorbar(x_input(2),y2,y2_se,y2_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(2) = plot(x_input(2),y2,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(2,:));

    errorbar(x_input(3),y3,y3_se,y3_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(3) = plot(x_input(3),y3,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(3,:));

    errorbar(x_input(4),y4,y4_se,y4_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(4) = plot(x_input(4),y4,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(4,:));

    errorbar(x_input(5),y5,y5_se,y5_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(5) = plot(x_input(5),y5,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(5,:));

    errorbar(x_input(6),y6,y6_se,y6_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(6) = plot(x_input(6),y6,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(6,:));

    errorbar(x_input(7),y7,y7_se,y7_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(7) = plot(x_input(7),y7,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(7,:));

h(8) = plot(x_input,y_input,'o-','MarkerSize',12,'Color','k','LineWidth',2);

hold off
xlim([14 46])
ylim([0 0.9])
set(gca,'YTick',0: 0.2 :0.9,'ycolor','k')
set(gca,'XTick',15: 5 : 45,'ycolor','k')
xlabel('{\itT}_L_e_a_f (^{o}C)','Color','k');
ylabel('{\itq}_L_\__P_A_M','Color','k')


% lgd = legend(h(1:7),'15 ^{o}C','20 ^{o}C','25 ^{o}C','30 ^{o}C','35 ^{o}C','40 ^{o}C','45 ^{o}C','Location','south','NumColumns',2,'FontSize',18);
% legend('boxoff')

txt = {'(b)'};
text(0.32*max(xlim),0.93*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)

txt = {'Low PAR (100 \mumol m^-^2 s^-^1)'};
text(0.36*max(xlim),1.08*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
box on
grid on

% ah1=axes('position',get(gca,'position'),'visible','off');
% lh2 = legend(ah1, h(1:7),'15 ^{o}C','20 ^{o}C','25 ^{o}C','30 ^{o}C','35 ^{o}C','40 ^{o}C','45 ^{o}C','Location','southoutside','NumColumns',7,'FontSize',22,'Orientation','horizontal');
% legend('boxoff')
% lh2.Position = [0.24 0.97 0.7583 0.0328];% [0.25 0.005 0.7583 0.0328];
% set(lh2,'Fontname','times new Roman','FontWeight','bold','FontSize',22);


    
h2 = subplot(1,3,2);
set(gca,'position', [0.385 0.21 0.26 0.7]);   % [x0 y0 width height]

x_input = [15 20 25 30 35 40 45];
i=7;
j=1; y1 =mean(qL_cal(i+1,j:j+3),2); y1_se = std(qL_cal(i+1,j:j+3));
j=5; y2 =mean(qL_cal(i+1,j:j+3),2); y2_se = std(qL_cal(i+1,j:j+3));
j=9; y3 =mean(qL_cal(i+1,j:j+3),2); y3_se = std(qL_cal(i+1,j:j+3));
j=13; y4 =mean(qL_cal(i+1,j:j+3),2); y4_se = std(qL_cal(i+1,j:j+3));
j=17; y5 =mean(qL_cal(i+1,j:j+3),2); y5_se = std(qL_cal(i+1,j:j+3));
j=21; y6 =mean(qL_cal(i+1,j:j+3),2); y6_se = std(qL_cal(i+1,j:j+3));
j=25; y7 =mean(qL_cal(i+1,j:j+3),2); y7_se = std(qL_cal(i+1,j:j+3));
y_input = [y1 y2 y3 y4 y5 y6 y7];

clear h;
    errorbar(x_input(1),y1,y1_se,y1_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
hold on
h(1) = plot(x_input(1),y1,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(1,:));

    errorbar(x_input(2),y2,y2_se,y2_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(2) = plot(x_input(2),y2,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(2,:));

    errorbar(x_input(3),y3,y3_se,y3_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(3) = plot(x_input(3),y3,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(3,:));

    errorbar(x_input(4),y4,y4_se,y4_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(4) = plot(x_input(4),y4,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(4,:));

    errorbar(x_input(5),y5,y5_se,y5_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(5) = plot(x_input(5),y5,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(5,:));

    errorbar(x_input(6),y6,y6_se,y6_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(6) = plot(x_input(6),y6,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(6,:));

    errorbar(x_input(7),y7,y7_se,y7_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(7) = plot(x_input(7),y7,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(7,:));

h(8) = plot(x_input,y_input,'o-','MarkerSize',12,'Color','k','LineWidth',2);

hold off
xlim([14 46])
ylim([0 0.9])
set(gca,'YTick',0: 0.2 :0.9,'ycolor','k')
set(gca,'XTick',15: 5 : 45,'ycolor','k')
xlabel('{\itT}_L_e_a_f (^{o}C)','Color','k');
ylabel('{\itq}_L_\__P_A_M','Color','k')

txt = {'(c) '};
text(0.32*max(xlim),0.93*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)

txt = {'Intermediate PAR (900 \mumol m^-^2 s^-^1)'};
text(0.28*max(xlim),1.08*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
box on
grid on




h3 = subplot(1,3,3);
set(gca,'position', [0.71 0.21 0.26 0.7]);   % [x0 y0 width height]


x_input = [15 20 25 30 35 40 45];
i=10;
j=1; y1 =mean(qL_cal(i+1,j:j+3),2); y1_se = std(qL_cal(i+1,j:j+3));
j=5; y2 =mean(qL_cal(i+1,j:j+3),2); y2_se = std(qL_cal(i+1,j:j+3));
j=9; y3 =mean(qL_cal(i+1,j:j+3),2); y3_se = std(qL_cal(i+1,j:j+3));
j=13; y4 =mean(qL_cal(i+1,j:j+3),2); y4_se = std(qL_cal(i+1,j:j+3));
j=17; y5 =mean(qL_cal(i+1,j:j+3),2); y5_se = std(qL_cal(i+1,j:j+3));
j=21; y6 =mean(qL_cal(i+1,j:j+3),2); y6_se = std(qL_cal(i+1,j:j+3));
j=25; y7 =mean(qL_cal(i+1,j:j+3),2); y7_se = std(qL_cal(i+1,j:j+3));
y_input = [y1 y2 y3 y4 y5 y6 y7];

clear h;
    errorbar(x_input(1),y1,y1_se,y1_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
hold on
h(1) = plot(x_input(1),y1,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(1,:));

    errorbar(x_input(2),y2,y2_se,y2_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(2) = plot(x_input(2),y2,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(2,:));

    errorbar(x_input(3),y3,y3_se,y3_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(3) = plot(x_input(3),y3,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(3,:));

    errorbar(x_input(4),y4,y4_se,y4_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(4) = plot(x_input(4),y4,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(4,:));

    errorbar(x_input(5),y5,y5_se,y5_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(5) = plot(x_input(5),y5,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(5,:));

    errorbar(x_input(6),y6,y6_se,y6_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(6) = plot(x_input(6),y6,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(6,:));

    errorbar(x_input(7),y7,y7_se,y7_se,0,0,'o','MarkerSize',12,'Color','k','MarkerFaceColor','k','LineWidth',1);
h(7) = plot(x_input(7),y7,'o','MarkerSize',12,'Color','k','MarkerFaceColor',color_temp(7,:));

h(8) = plot(x_input,y_input,'o-','MarkerSize',12,'Color','k','LineWidth',2);

hold off
xlim([14 46])
ylim([0 0.9])
set(gca,'YTick',0: 0.2 :0.9,'ycolor','k')
set(gca,'XTick',15: 5 : 45,'ycolor','k')
xlabel('{\itT}_L_e_a_f (^{o}C)','Color','k');
ylabel('{\itq}_L_\__P_A_M','Color','k')

txt = {'(d)'};
text(0.32*max(xlim),0.93*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)

txt = {'High PAR  (1800 \mumol m^-^2 s^-^1)'};
text(0.33*max(xlim),1.08*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',22)

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(gca,'LooseInset',get(gca,'TightInset'));
box on
grid on

%print('Fig4-2 qL_TChlF','-djpeg','-r300')



%% Fig5-1  qL_PAM vs T_leaf

clear; clc;

load leaf_dataset.mat

color_temp =jet(7);

% using mean value of qL to fit with PAR and TChlF, then compare their R2 with each other
PAR(1,:)=0;



ax1= figure('visible','on');
set(gcf,'Position',[50 50 800 800]);

% PAR vs qL
%h1 = subplot(1,2,1);
set(gca,'position', [0.12 0.14 0.84 0.84]);   % [x0 y0 width height]


temp_i = 0;
qL_mean_wheat_15      = mean(qL_cal(1:12,temp_i+1:temp_i+4),2,'omitnan');
qL_se_wheat_15           = std(qL_cal(1:12,temp_i+1:temp_i+4),0,2,'omitnan');
chlf_mean_wheat_15      = mean(TChlF(1:12,temp_i+1:temp_i+4),2,'omitnan');
chlf_se_wheat_15           = std(TChlF(1:12,temp_i+1:temp_i+4),0,2,'omitnan');

x =chlf_mean_wheat_15;
x_std =chlf_se_wheat_15;
y = qL_mean_wheat_15;
y_std = qL_se_wheat_15;

%plot(TChlF(1:12,temp_i+1:temp_i+4),qL_cal(1:12,temp_i+1:temp_i+4),'.','MarkerSize',12,'Color',1/255.*[163 200 255],'MarkerFaceColor','#0072BD','LineWidth',1);
errorbar(x,y,y_std,y_std,x_std,x_std,'o','MarkerSize',8,'Color',color_temp(1,:),'MarkerFaceColor',color_temp(1,:));

hold on

h(1) = plot(x,y,'o','MarkerSize',8,'Color','k','MarkerFaceColor',color_temp(1,:));
[xData, yData] = prepareCurveData(x, y);
%Set up fittype and options.
%ft = fittype( '1/(a.*x.^a+1)', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.2;
%Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
m_wheat_15 = fitresult.a;
r2_wheat_15= gof.rsquare;
rmse_wheat_15= 100.*gof.rmse./(max(yData)-min(yData)); % rrmse
%rrmse_wheat_15 = 100.*gof.rmse./(max(yData)-min(yData));

xData = [min(xData):(max(xData)-min(xData))./1000:max(xData)]';
Yhat = feval(fitresult, xData);
plot(xData,Yhat,'Color',color_temp(1,:),'LineWidth',1.5);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);
qL_mod_wheat_15 = feval(fitresult, x);


temp_i = 4;
qL_mean_wheat_20      = mean(qL_cal(1:12,temp_i+1:temp_i+4),2,'omitnan');
qL_se_wheat_20           = std(qL_cal(1:12,temp_i+1:temp_i+4),0,2,'omitnan');
chlf_mean_wheat_20      = mean(TChlF(1:12,temp_i+1:temp_i+4),2,'omitnan');
chlf_se_wheat_20           = std(TChlF(1:12,temp_i+1:temp_i+4),0,2,'omitnan');

x =chlf_mean_wheat_20;
x_std =chlf_se_wheat_20;
y = qL_mean_wheat_20;
y_std = qL_se_wheat_20;

%plot(TChlF(1:12,temp_i+1:temp_i+4),qL_cal(1:12,temp_i+1:temp_i+4),'.','MarkerSize',12,'Color',1/255.*[255 215 155],'MarkerFaceColor','#EDB120','LineWidth',1);

errorbar(x,y,y_std,y_std,x_std,x_std,'s','MarkerSize',8,'Color',color_temp(2,:),'MarkerFaceColor',color_temp(2,:));

h(2) = plot(x,y,'s','MarkerSize',8,'Color','k','MarkerFaceColor',color_temp(2,:));

[xData, yData] = prepareCurveData(x, y);
%Set up fittype and options.
ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.2;
%Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
m_wheat_20 = fitresult.a;
r2_wheat_20= gof.rsquare;
rmse_wheat_20= 100.*gof.rmse./(max(yData)-min(yData)); % rrmse
xData = [min(xData):(max(xData)-min(xData))./1000:max(xData)]';
Yhat = feval(fitresult, xData);
plot(xData,Yhat,'Color',color_temp(2,:),'LineWidth',1.5);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);
qL_mod_wheat_20 = feval(fitresult, x);


temp_i = 8;
qL_mean_wheat_25      = mean(qL_cal(1:12,temp_i+1:temp_i+4),2,'omitnan');
qL_se_wheat_25           = std(qL_cal(1:12,temp_i+1:temp_i+4),0,2,'omitnan');
chlf_mean_wheat_25      = mean(TChlF(1:12,temp_i+1:temp_i+4),2,'omitnan');
chlf_se_wheat_25           = std(TChlF(1:12,temp_i+1:temp_i+4),0,2,'omitnan');

x =chlf_mean_wheat_25;
x_std =chlf_se_wheat_25;
y = qL_mean_wheat_25;
y_std = qL_se_wheat_25;

%plot(TChlF(1:12,temp_i+1:temp_i+4),qL_cal(1:12,temp_i+1:temp_i+4),'.','MarkerSize',12,'Color',1/255.*[255 170 170],'MarkerFaceColor','#D95319','LineWidth',1);

errorbar(x,y,y_std,y_std,x_std,x_std,'d','MarkerSize',8,'Color',color_temp(3,:),'MarkerFaceColor',color_temp(3,:));
h(3) = plot(x,y,'d','MarkerSize',8,'Color','k','MarkerFaceColor',color_temp(3,:));

[xData, yData] = prepareCurveData(x, y);
%Set up fittype and options.
ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.2;
%Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
m_wheat_25 = fitresult.a;
r2_wheat_25= gof.rsquare;
rmse_wheat_25= 100.*gof.rmse./(max(yData)-min(yData)); % rrmse
xData = [min(xData):(max(xData)-min(xData))./1000:max(xData)]';
Yhat = feval(fitresult, xData);
plot(xData,Yhat,'Color',color_temp(3,:),'LineWidth',1.5);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);
qL_mod_wheat_25 = feval(fitresult, x);


temp_i = 12;
qL_mean_wheat_30      = mean(qL_cal(1:12,temp_i+1:temp_i+4),2,'omitnan');
qL_se_wheat_30           = std(qL_cal(1:12,temp_i+1:temp_i+4),0,2,'omitnan');
chlf_mean_wheat_30      = mean(TChlF(1:12,temp_i+1:temp_i+4),2,'omitnan');
chlf_se_wheat_30           = std(TChlF(1:12,temp_i+1:temp_i+4),0,2,'omitnan');

x =chlf_mean_wheat_30;
x_std =chlf_se_wheat_30;
y = qL_mean_wheat_30;
y_std = qL_se_wheat_30;

%plot(TChlF(1:12,temp_i+1:temp_i+4),qL_cal(1:12,temp_i+1:temp_i+4),'.','MarkerSize',12,'Color',1/255.*[255 170 170],'MarkerFaceColor','#D95319','LineWidth',1);

errorbar(x,y,y_std,y_std,x_std,x_std,'^','MarkerSize',8,'Color',color_temp(4,:),'MarkerFaceColor',color_temp(4,:));
h(4) = plot(x,y,'^','MarkerSize',8,'Color','k','MarkerFaceColor',color_temp(4,:));

[xData, yData] = prepareCurveData(x, y);
%Set up fittype and options.
ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.2;
%Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
m_wheat_30 = fitresult.a;
r2_wheat_30= gof.rsquare;
rmse_wheat_30= 100.*gof.rmse./(max(yData)-min(yData)); % rrmse
xData = [min(xData):(max(xData)-min(xData))./1000:max(xData)]';
Yhat = feval(fitresult, xData);
plot(xData,Yhat,'Color',color_temp(4,:),'LineWidth',1.5);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);
qL_mod_wheat_30 = feval(fitresult, x);


temp_i = 16;
qL_mean_wheat_35      = mean(qL_cal(1:12,temp_i+1:temp_i+4),2,'omitnan');
qL_se_wheat_35           = std(qL_cal(1:12,temp_i+1:temp_i+4),0,2,'omitnan');
chlf_mean_wheat_35      = mean(TChlF(1:12,temp_i+1:temp_i+4),2,'omitnan');
chlf_se_wheat_35           = std(TChlF(1:12,temp_i+1:temp_i+4),0,2,'omitnan');

x =chlf_mean_wheat_35;
x_std =chlf_se_wheat_35;
y = qL_mean_wheat_35;
y_std = qL_se_wheat_35;

%plot(TChlF(1:12,temp_i+1:temp_i+4),qL_cal(1:12,temp_i+1:temp_i+4),'.','MarkerSize',12,'Color',1/255.*[255 170 170],'MarkerFaceColor','#D95319','LineWidth',1);

errorbar(x,y,y_std,y_std,x_std,x_std,'>','MarkerSize',8,'Color',color_temp(5,:),'MarkerFaceColor',color_temp(5,:));
h(5) = plot(x,y,'>','MarkerSize',8,'Color','k','MarkerFaceColor',color_temp(5,:));

[xData, yData] = prepareCurveData(x, y);
%Set up fittype and options.
ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.2;
%Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
m_wheat_35 = fitresult.a;
r2_wheat_35= gof.rsquare;
rmse_wheat_35= 100.*gof.rmse./(max(yData)-min(yData)); % rrmse
xData = [min(xData):(max(xData)-min(xData))./1000:1.1*max(xData)]';
Yhat = feval(fitresult, xData);
plot(xData,Yhat,'Color',color_temp(5,:),'LineWidth',2);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);
%plot(xData,Yhat,'--k','LineWidth',0.5);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);
qL_mod_wheat_35 = feval(fitresult, x);


temp_i = 20;
qL_mean_wheat_40      = mean(qL_cal(1:12,temp_i+1:temp_i+4),2,'omitnan');
qL_se_wheat_40           = std(qL_cal(1:12,temp_i+1:temp_i+4),0,2,'omitnan');
chlf_mean_wheat_40      = mean(TChlF(1:12,temp_i+1:temp_i+4),2,'omitnan');
chlf_se_wheat_40           = std(TChlF(1:12,temp_i+1:temp_i+4),0,2,'omitnan');

x =chlf_mean_wheat_40;
x_std =chlf_se_wheat_40;
y = qL_mean_wheat_40;
y_std = qL_se_wheat_40;

%plot(TChlF(1:12,temp_i+1:temp_i+4),qL_cal(1:12,temp_i+1:temp_i+4),'.','MarkerSize',12,'Color',1/255.*[255 170 170],'MarkerFaceColor','#D95319','LineWidth',1);

errorbar(x,y,y_std,y_std,x_std,x_std,'p','MarkerSize',8,'Color',color_temp(6,:),'MarkerFaceColor',color_temp(6,:));
h(6) = plot(x,y,'p','MarkerSize',8,'Color','k','MarkerFaceColor',color_temp(6,:));

[xData, yData] = prepareCurveData(x, y);
%Set up fittype and options.
ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.2;
%Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
m_wheat_40 = fitresult.a;
r2_wheat_40= gof.rsquare;
rmse_wheat_40= 100.*gof.rmse./(max(yData)-min(yData)); % rrmse
xData = [min(xData):(max(xData)-min(xData))./1000:max(xData)]';
Yhat = feval(fitresult, xData);
plot(xData,Yhat,'Color',color_temp(6,:),'LineWidth',1.5);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);
qL_mod_wheat_40=feval(fitresult, x);


temp_i = 24;
qL_mean_wheat_45      = mean(qL_cal(1:12,temp_i+1:temp_i+4),2,'omitnan');
qL_se_wheat_45           = std(qL_cal(1:12,temp_i+1:temp_i+4),0,2,'omitnan');
chlf_mean_wheat_45      = mean(TChlF(1:12,temp_i+1:temp_i+4),2,'omitnan');
chlf_se_wheat_45           = std(TChlF(1:12,temp_i+1:temp_i+4),0,2,'omitnan');

x =chlf_mean_wheat_45;
x_std =chlf_se_wheat_45;
y = qL_mean_wheat_45;
y_std = qL_se_wheat_45;

%plot(TChlF(1:12,temp_i+1:temp_i+4),qL_cal(1:12,temp_i+1:temp_i+4),'.','MarkerSize',12,'Color',1/255.*[255 170 170],'MarkerFaceColor','#D95319','LineWidth',1);

errorbar(x,y,y_std,y_std,x_std,x_std,'h','MarkerSize',8,'Color',color_temp(7,:),'MarkerFaceColor',color_temp(7,:));
h(7) = plot(x,y,'h','MarkerSize',8,'Color','k','MarkerFaceColor',color_temp(7,:));

[xData, yData] = prepareCurveData(x, y);
%Set up fittype and options.
ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.2;
%Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
m_wheat_45 = fitresult.a;
r2_wheat_45= gof.rsquare;
rmse_wheat_45= 100.*gof.rmse./(max(yData)-min(yData)); % rrmse
xData = [min(xData):(max(xData)-min(xData))./1000:max(xData)]';
Yhat = feval(fitresult, xData);
plot(xData,Yhat,'Color',color_temp(7,:),'LineWidth',1.5);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);
qL_mod_wheat_45=feval(fitresult, x);


xlim([0 45])
set(gca,'XTick',[0: 10 :60],'ycolor','k')
ylim([0.0 1])
set(gca,'YTick',[0: 0.2 :1],'ycolor','k')

hold off

xlabel(['SIF_F_U_L_L_\__P_S_I_I_\__L (\mumol m^-^2 s^-^1)'],'Color','k')% ChlFTOT_FULL_PSII
ylabel(['{\itq}_L_\__P_A_M'],'Color','k')


lgd = legend(h(1:7),...
    ['  15 ^{o}C ', '    ',num2str(m_wheat_15,'%6.2f'), '   ',num2str(r2_wheat_15,'%6.2f'), '   ',num2str(rmse_wheat_15,'%6.2f'),'%  '],...
    ['  20 ^{o}C ', '    ',num2str(m_wheat_20,'%6.2f'), '   ',num2str(r2_wheat_20,'%6.2f'), '   ',num2str(rmse_wheat_20,'%6.2f'),'%  '],...
    ['  25 ^{o}C ', '    ',num2str(m_wheat_25,'%6.2f'), '   ',num2str(r2_wheat_25,'%6.2f'), '   ',num2str(rmse_wheat_25,'%6.2f'),'%  '],...
    ['  30 ^{o}C ', '    ',num2str(m_wheat_30,'%6.2f'), '   ',num2str(r2_wheat_30,'%6.2f'), '   ',num2str(rmse_wheat_30,'%6.2f'),'%  '],...
    ['  35 ^{o}C ', '    ',num2str(m_wheat_35,'%6.2f'), '   ',num2str(r2_wheat_35,'%6.2f'), '   ',num2str(rmse_wheat_35,'%6.2f'),'%  '],...
    ['  40 ^{o}C ', '    ',num2str(m_wheat_40,'%6.2f'), '   ',num2str(r2_wheat_40,'%6.2f'), '   ',num2str(rmse_wheat_40,'%6.2f'),'%  '],...
    ['  45 ^{o}C ', '    ',num2str(m_wheat_45,'%6.2f'), '   ',num2str(r2_wheat_45,'%6.2f'), '   ',num2str(rmse_wheat_45,'%6.2f'),'%  '],...
    'Location','northeast','NumColumns',1,'FontSize',30);
title(lgd,{['{\it{y}} = m/({\it{x}}^{1/m}+m)  m     {\it{R}}^2  rRMSE']},'FontSize',30);

legend('boxoff')

txt = {'(a)'};
text(0.04*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',34)

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',32);
set(gca,'LooseInset',get(gca,'TightInset'));
box on
grid   on

%print('Fig5-1','-djpeg','-r300')


qL_mean_wheat = [qL_mean_wheat_15 qL_mean_wheat_20 qL_mean_wheat_25 qL_mean_wheat_30 qL_mean_wheat_35 qL_mean_wheat_40 qL_mean_wheat_45];
qL_mod_wheat = [qL_mod_wheat_15  qL_mod_wheat_20 qL_mod_wheat_25 qL_mod_wheat_30 qL_mod_wheat_35 qL_mod_wheat_40 qL_mod_wheat_45];



%% Fig5-2 Soybean T_leaf vs 1/m
m_qL_CHLF = nan(3,size(qL_cal,2));

for jj = 1: size(qL_cal,2)

    x =TChlF(1:12,jj);  
    y = qL_cal(1:12,jj);

    [xData, yData] = prepareCurveData(x, y);
    %Set up fittype and options.
    ft = fittype( 'a/(x.^(1/a)+a)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';   
    opts.StartPoint = 0.2;
    %Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    m_temp = fitresult.a;
    r2_temp= gof.rsquare;
    rmse_temp= 100.*gof.rmse./(max(yData)-min(yData)); % rrmse
    Yhat = feval(fitresult, xData);
    %plot(xData,Yhat,'Color','#A2142F','LineWidth',1.5);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);

    m_qL_CHLF(1,jj)=m_temp;
    m_qL_CHLF(2,jj)=r2_temp;
    m_qL_CHLF(3,jj)=rmse_temp;

end

% ax1= figure('visible','on');
% set(gcf,'Position',[50 50 1400 650]);



ax1= figure('visible','on');
set(gcf,'Position',[50 50 800 800]);

% PAR vs qL
%h1 = subplot(1,2,1);
set(gca,'position', [0.12 0.14 0.84 0.84]);   % [x0 y0 width height]




Temperature_Leaf_single =mean(Temperature_Leaf(1:12,1:end));
Temperature_Leaf_single1 = reshape(Temperature_Leaf_single,[4, 7])';
x = mean(Temperature_Leaf_single1(1:7,:),2,'omitnan');
x_std     = std(Temperature_Leaf_single1(1:7,:),0,2,'omitnan');

m_qL_CHLF1 =  reshape(m_qL_CHLF(1,:),[4, 7])';
y = mean(m_qL_CHLF1(1:7,:),2,'omitnan');
y_std     = std(m_qL_CHLF1(1:7,:),0,2,'omitnan');

%plot(Temperature_Leaf_single(:,29:56),m_qL_CHLF(1,29:56),'.','MarkerSize',12,'Color',1/255.*[200 200 200]);

errorbar(x,y,y_std,y_std,x_std,x_std,'o','MarkerSize',10,'Color','k','MarkerFaceColor','k','LineWidth',1.5);
hold on
h(1) = plot(x,y,'o','MarkerSize',10,'Color','k','MarkerFaceColor','k');

% x_input = Temperature_Leaf_single(29:56);
% y_input =m_qL_CHLF(1,29:56);

xlabel(['{\itT}_L_e_a_f (^{o}C)'],'Color','k')
ylabel(['{\itm}'],'Color','k')

xlim([14 46])
% ylim([1 2.3])
% set(gca,'YTick',[1: 0.3 :2.3],'ycolor','k')
set(gca,'XTick',[15: 5 : 45],'ycolor','k')

set(gca,'ycolor','k');


x_input = x+273.15;
y_input =y;

modelfun = @(b,x)b(1).*(b(2).*exp((b(3).*(x-b(4))./(x.*8.314.*b(4)))))./(b(2)-b(3).*(1-exp((b(2).*(x-b(4))./(x.*8.314.*b(4))))));
beta0 = [50 200000 43360 273.15+35.24];

% modelfun = @(b,x)b(1).*(200000.*exp((b(2).*(x-b(3))./(x.*8.314.*b(3)))))./(200000-b(2).*(1-exp((200000.*(x-b(3))./(x.*8.314.*b(3))))));
% beta0 = [50 43360 273.15+35.24];


model_fit = fitnlm(x_input,y_input, modelfun, beta0);
SE = diag(sqrt(model_fit.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
[ypred,delta_t] = predict(model_fit,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange-273.15,ypred,'k-','LineWidth',2)
boundedline(xrange-273.15, ypred, delta,'alpha');
%plot(xrange-273.15,delta_t,'r--','LineWidth',1)
plot(xrange-273.15,ypred,'k-','LineWidth',2)
hold off

b_wheat = model_fit.Coefficients.Estimate;

r2_temp= model_fit.Rsquared.Ordinary;
%rmse_temp= model_fit.RMSE;
rmse_temp= 100.*model_fit.RMSE./(max(y_input)-min(y_input)); % rrmse

%Yhat = feval(fitresult, xData);

dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
lgd = legend(dummyh, [''],'Location','southeast','NumColumns',1,'FontSize',14);

legend('boxoff')

R22=roundn(r2_temp,-4);
RMSE2=roundn(rmse_temp ,-4);
%title(lgd,{['Soybean'],['{\itq}_L_\__m_o_d = 1/({\itm}ChlF_P_S_I_I ^{\itm}+1)'],['{\it{R}}^2 = ',num2str(R22,'%6.2f')],['RMSE = ',num2str(RMSE2,'%6.2f')]},'FontSize',18);
% title(lgd,{[''],['{\it{R}}^2 = ',num2str(R22,'%6.2f')],['RMSE = ',num2str(RMSE2,'%6.4f')]},'FontSize',18);

% lgd.Title.NodeChildren.Position =  [0.5 0.6 0];


txt = {[' {\it{m}}_o_p_t = ',num2str(model_fit.Coefficients.Estimate(1),'%.2f'),...
    ';       {\it{H}}_d = ',num2str(model_fit.Coefficients.Estimate(2),'%.2e'),' J mol^-^1'];...
    ['{\it{T}}_o_p_t = ',num2str(model_fit.Coefficients.Estimate(4)-273.15,'%.2f'),' ^{o}C;',...
    ' {\it{H}}_a = ',num2str(model_fit.Coefficients.Estimate(3),'%.2e'),' J mol^-^1'...
    ]};
text(0.4*max(xlim),0.25*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',30)

txt = {'',['{\it{R}}^2 = ',num2str(R22,'%6.2f')],['rRMSE = ',num2str(RMSE2,'%6.2f'),'%']};
text(0.35*max(xlim),0.90*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',30)

b = model_fit.Coefficients.Estimate;
save qL_4parameters.mat  b

txt = {'(b) '};
text(0.35*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',34)
box on
grid   on

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',32);
set(gca,'LooseInset',get(gca,'TightInset'));

%print('Fig5-2','-djpeg','-r300')



%% FigureS1 PAR vs f_PSII_752


ax1= figure('visible','on');
set(gcf,'Position',[50 50 800 700]);
h1 = subplot(1,2,1);
set(gca,'position', [0.15 0.17 0.8 0.8]);   % [x0 y0 width height]
clear h

% corn

i =1;
h(1) = plot(reshape(PAR(1:12,i:i+3),[],1), reshape(fPSII_760_set(1:12,i:i+3),[],1),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(1,:),'LineWidth',1.5);
hold on
i =5;
h(2) = plot(reshape(PAR(1:12,i:i+3),[],1), reshape(fPSII_760_set(1:12,i:i+3),[],1),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(2,:),'LineWidth',1.5);
i =9;
h(3) = plot(reshape(PAR(1:12,i:i+3),[],1), reshape(fPSII_760_set(1:12,i:i+3),[],1),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(3,:),'LineWidth',1.5);
i =13;
h(4) = plot(reshape(PAR(1:12,i:i+3),[],1), reshape(fPSII_760_set(1:12,i:i+3),[],1),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(4,:),'LineWidth',1.5);
i =17;
h(5) = plot(reshape(PAR(1:12,i:i+3),[],1), reshape(fPSII_760_set(1:12,i:i+3),[],1),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(5,:),'LineWidth',1.5);
i =21;
h(6) = plot(reshape(PAR(1:12,i:i+3),[],1), reshape(fPSII_760_set(1:12,i:i+3),[],1),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(6,:),'LineWidth',1.5);
i =25;
h(7) = plot(reshape(PAR(1:12,i:i+3),[],1), reshape(fPSII_760_set(1:12,i:i+3),[],1),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(7,:),'LineWidth',1.5);



i =1;
h(1) = plot(mean(PAR(1:12,i:i+3),2), mean(fPSII_760_set(1:12,i:i+3),2),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(1,:),'LineWidth',1.5);
hold on
i =5;
h(1) = plot(mean(PAR(1:12,i:i+3),2), mean(fPSII_760_set(1:12,i:i+3),2),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(2,:),'LineWidth',1.5);
i =9;
h(1) = plot(mean(PAR(1:12,i:i+3),2), mean(fPSII_760_set(1:12,i:i+3),2),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(3,:),'LineWidth',1.5);
i =13;
h(1) = plot(mean(PAR(1:12,i:i+3),2), mean(fPSII_760_set(1:12,i:i+3),2),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(4,:),'LineWidth',1.5);
i =17;
h(1) = plot(mean(PAR(1:12,i:i+3),2), mean(fPSII_760_set(1:12,i:i+3),2),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(5,:),'LineWidth',1.5);
i =21;
h(1) = plot(mean(PAR(1:12,i:i+3),2), mean(fPSII_760_set(1:12,i:i+3),2),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(6,:),'LineWidth',1.5);
i =25;
h(1) = plot(mean(PAR(1:12,i:i+3),2), mean(fPSII_760_set(1:12,i:i+3),2),'o','MarkerSize',6,'Color','k','MarkerFaceColor',color_temp(7,:),'LineWidth',1.5);




% xData = reshape(TChlF(1:12,1:28),[],1);
% yData=reshape(ChlF_PAM(1:12,1:28),[],1);
% fitresult = fitlm(xData,yData);
% a_maize = fitresult.Coefficients.Estimate(2);
% b_maize = fitresult.Coefficients.Estimate(1);
% r2_maize= fitresult.Rsquared.Ordinary;
% rmse_maize= fitresult.RMSE;
% xData = [min(xData):(max(xData)-min(xData))./1000:max(xData)]';
% Yhat = feval(fitresult, xData);
% plot(xData,Yhat,'Color','k','LineWidth',2);


fPSII_mean = mean(reshape(fPSII_760_set(2:12,1:i+3),[],1));
fPSII_std = std(reshape(fPSII_760_set(2:12,1:i+3),[],1));
plot([0 2100],[fPSII_mean+fPSII_std fPSII_mean+fPSII_std],'r--','LineWidth',1)
plot([0 2100],[fPSII_mean fPSII_mean],'k--','LineWidth',2)
plot([0 2100],[fPSII_mean-fPSII_std fPSII_mean-fPSII_std],'r--','LineWidth',1)


hold off

ylabel(['{\it{f}}_7_6_0_\__P_S_I_I'],'Color','k')% ChlFTOT_FULL_PSII
xlabel(['PAR (\mumol m^-^2 s^-^1)'],'Color','k')

xlim([0 2200])
set(gca,'XTick',[0: 300 :2400],'ycolor','k')
ylim([0.4 0.9])
set(gca,'YTick',[0.0 : 0.1 :1],'ycolor','k')


lgd = legend(h(1:7),'15 ^{o}C','20 ^{o}C','25 ^{o}C','30 ^{o}C','35 ^{o}C','40 ^{o}C','45 ^{o}C','Location','northeast','NumColumns',4,'FontSize',18);
legend('boxoff')

% if fitresult.Coefficients.pValue(1)<0.05
%     txt = {...
%     ['y-intercept = ',num2str(fitresult.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(fitresult.Coefficients.SE(1),'%10.2g')],...
%         ['slope = ',num2str(fitresult.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(fitresult.Coefficients.SE(2),'%10.2g')],...
%         };
% else
%     txt = {...
%         ['y-intercept^* = ',num2str(fitresult.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(fitresult.Coefficients.SE(1),'%10.2g')],...
%         ['slope = ',num2str(fitresult.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(fitresult.Coefficients.SE(2),'%10.2g')],...
%         };
% end
% text(0.5*max(xlim),0.1*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)


% txt = {...
%     ['{\it{R}}^2 = ',num2str(r2_maize,'%6.2f')], ['RMSE = ',num2str(rmse_maize,'%6.2f')],...
%     };
% text(0.7*max(xlim),0.3*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
%
% txt = {...
%     ['(a) Corn'],...
%     };
% text(0.05*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',24)


set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
box on
grid   on

%print('FigureS1','-djpeg','-r300')


%% Figure S2 temperature_qL_chlF_qL_PAM

ax1= figure('visible','on');
set(gcf,'Position',[50 50 800 800]);

% PAR vs qL
%h1 = subplot(1,2,1);
set(gca,'position', [0.12 0.14 0.84 0.84]);   % [x0 y0 width height]

h(2) = plot(qL_mod_wheat(:,1), qL_mean_wheat(:,1),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(1,:));
hold on

h(3) = plot(qL_mod_wheat(:,2), qL_mean_wheat(:,2),'s','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(2,:));

h(4) = plot(qL_mod_wheat(:,3), qL_mean_wheat(:,3),'d','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(3,:));

h(5) = plot(qL_mod_wheat(:,4), qL_mean_wheat(:,4),'^','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(4,:));

h(6) = plot(qL_mod_wheat(:,5), qL_mean_wheat(:,5),'>','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(5,:));

h(7) = plot(qL_mod_wheat(:,6), qL_mean_wheat(:,6),'p','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(6,:));

h(8) = plot(qL_mod_wheat(:,7), qL_mean_wheat(:,7),'h','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(7,:));

xData = reshape(qL_mod_wheat(1:12,1:7),[],1);
yData=reshape(qL_mean_wheat(1:12,1:7),[],1);
fitresult = fitlm(xData,yData);
a_maize = fitresult.Coefficients.Estimate(2);
b_maize = fitresult.Coefficients.Estimate(1);
r2_maize= fitresult.Rsquared.Ordinary;
rmse_maize= 100.*fitresult.RMSE./(max(yData)-min(yData)); % rrmse

xData = [min(xData):(max(xData)-min(xData))./1000:max(xData)]';
Yhat = feval(fitresult, xData);
plot(xData,Yhat,'Color','k','LineWidth',2);
plot([0:0.01:1],[0:0.01:1],'r--','LineWidth',2)



h(1)  = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');

hold off

xlabel('{\itq}_L_\__S_I_F','Color','k')
ylabel('{\itq}_L_\__P_A_M','Color','k')

xlim([0 1])
set(gca,'XTick',0: 0.2 :1,'ycolor','k')
ylim([0.0 1])
set(gca,'YTick',0: 0.2 :1,'ycolor','k')


% lgd = legend(h(1:8),'Corn','15 ^{o}C','20 ^{o}C','25 ^{o}C','30 ^{o}C','35 ^{o}C','40 ^{o}C','45 ^{o}C','Location','northoutside','NumColumns',4,'FontSize',18,'Orientation','horizontal');
% legend('boxoff')
% lgd.Title.NodeChildren.Position =  [0.5 0.85 0];

if fitresult.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(fitresult.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(fitresult.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(fitresult.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(fitresult.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(fitresult.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(fitresult.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(fitresult.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(fitresult.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.4*max(xlim),0.1*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',30)

txt = {['{\it{R}}^2 = ',num2str(r2_maize,'%6.2f')], ['rRMSE = ',num2str(rmse_maize,'%6.2f'),'%']};
text(0.08*max(xlim),0.86*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',30)

txt = {'(c) '};
text(0.08*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',34)
box on
grid   on

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',32);
box on
grid   on




% Copy the axes and plot the second legned
% ah1=axes('position',get(gca,'position'),'visible','off');
% lh2 = legend(ah1, h(2:8),'15 ^{o}C','20 ^{o}C','25 ^{o}C','30 ^{o}C','35 ^{o}C','40 ^{o}C','45 ^{o}C','Location','southoutside','NumColumns',4,'FontSize',18,'Orientation','horizontal');
% legend('boxoff')
% lh2.Position = [0.15 0.93 0.7583 0.0328];% [0.25 0.005 0.7583 0.0328];
% set(lh2,'Fontname','times new Roman','FontWeight','bold','FontSize',20);
% %set(lgd,'Fontname','times new Roman','FontWeight','bold','FontSize',20);




%print(['FigureS2'],'-djpeg','-r300')


%% Figure Supporting qL-PAR

clc;
clear;
load leaf_dataset.mat
color_temp = jet(7);

figure('visible','on');
set(gcf,'Position',[50 50 800 600]);

h1 = subplot(1,1,1);

set(gca,'position', [0.11 0.15 0.85 0.8]);   % [x0 y0 width height]

qL_mean_maize       = mean(qL_cal(1:12,1:28),2,'omitnan');
qL_se_maize            = std(qL_cal(1:12,1:28),0,2,'omitnan');
TChlF_mean_maize = mean(TChlF(1:12,1:28),2,'omitnan');
TChlF_se_maize      = std(TChlF(1:12,1:28),0,2,'omitnan');

% x =TChlF_mean_maize;
% y = qL_mean_maize;
% x_std =TChlF_se_maize;
% y_std = qL_se_maize;

x = [mean(PAR(1:12,1:4),2)
       mean(PAR(1:12,5:8),2)
       mean(PAR(1:12,9:12),2)
       mean(PAR(1:12,13:16),2)
       mean(PAR(1:12,17:20),2)
       mean(PAR(1:12,21:24),2)
       mean(PAR(1:12,25:28),2)];

y = [mean(qL_cal(1:12,1:4),2)
       mean(qL_cal(1:12,5:8),2)
       mean(qL_cal(1:12,9:12),2)
       mean(qL_cal(1:12,13:16),2)
       mean(qL_cal(1:12,17:20),2)
       mean(qL_cal(1:12,21:24),2)
       mean(qL_cal(1:12,25:28),2)];


%plot(TChlF(1:12,1:28),qL_cal(1:12,1:28),'.','MarkerSize',12,'Color',1/255.*[230 230 230],'MarkerFaceColor',1/255.*[230 230 230],'LineWidth',1);

i =1;
h(2) = plot(mean(PAR(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(1,:));
hold on
i =5;
h(3) = plot(mean(PAR(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(2,:));
i =9;
h(4) = plot(mean(PAR(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(3,:));
i =13;
h(5) = plot(mean(PAR(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(4,:));
i =17;
h(6) = plot(mean(PAR(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(5,:));
i =21;
h(7) = plot(mean(PAR(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(6,:));
i =25;
h(8) = plot(mean(PAR(1:12,i:i+3),2), mean(qL_cal(1:12,i:i+3),2),'o','MarkerSize',10,'Color','k','MarkerFaceColor',color_temp(7,:));

% errorbar(x,y,y_std,y_std,x_std,x_std,'d','MarkerSize',14,'Color','k','MarkerFaceColor','k','LineWidth',0.8);
% g(3) = plot(x,y,'d','MarkerSize',14,'Color','k','MarkerFaceColor','k');

[xData, yData] = prepareCurveData(x, y);
ft = fittype( 'a.*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0];
%Fit model to data.
[fitresult, gof] = fit(xData, yData, ft, opts );
a_PAR_QL_maize = fitresult.a;
b_PAR_QL_maize = fitresult.b;

r2_PAR_QL_maize= gof.rsquare;
rmse_PAR_QL_maize= gof.rmse;
rrmse = 100.*rmse_PAR_QL_maize./(max(yData)-min(yData));
xData = (min(xData):(max(xData)-min(xData))/1000:max(xData))';
Yhat = feval(fitresult, xData);
g(4)=plot(xData,Yhat,'Color','k','LineWidth',2);      %     Yhat = feval(fitresult, xData)';  %     fit_15=fitlm(y,Yhat);

xlim([0 2200])
set(gca,'XTick',0: 300 :2100,'ycolor','k')
ylim([0.0 1])
set(gca,'YTick',0: 0.2 :1,'ycolor','k')

g(1) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(2) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(3) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(5) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(6) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(7) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(8) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
g(9) = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');

hold off
xlabel(['PAR (\mumol m^-^2 s^-^1)'],'Color','k')% ChlFTOT_FULL_PSII
ylabel('{\itq}_L_\__P_A_M','Color','k')

% all pvalue > 0.01
legend(g(1:9),'','', '','Fit line',...
    ['{\it{y}} = {\it{a}}{\it{e}}^{-{\it{b}}{\it{x}}}'],...
    ['{\it{a}} = ',num2str(a_PAR_QL_maize,'%6.2f'),'; {\it{b}} = ',num2str(b_PAR_QL_maize,'%6.2e')],...
    ['{\it{R}}^2 = ',num2str(r2_PAR_QL_maize,'%6.2f')],...
    ['rRMSE = ',num2str(rrmse,'%6.2f'),'%'],...
    'Location','northeast','NumColumns',1,'FontSize',18);
legend('boxoff')

% txt = {'(a) Corn'};
% text(0.8*max(xlim),0.93*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',28)

set(h1,'Fontname','times new Roman','FontWeight','bold','FontSize',22);
set(h1,'LooseInset',get(gca,'TightInset'));


% Copy the axes and plot the second legned
ah1=axes('position',get(gca,'position'),'visible','off');
lh2 = legend(ah1, h(2:8),'15 ^{o}C','20 ^{o}C','25 ^{o}C','30 ^{o}C','35 ^{o}C','40 ^{o}C','45 ^{o}C','Location','southoutside','NumColumns',4,'FontSize',18,'Orientation','horizontal');
legend('boxoff')
lh2.Position = [0.15 0.87 0.7583 0.0328];% [0.25 0.005 0.7583 0.0328];
set(lh2,'Fontname','times new Roman','FontWeight','bold','FontSize',20);
%set(lgd,'Fontname','times new Roman','FontWeight','bold','FontSize',20);


grid(h1,'on')
box(h1,'on')
%print('FigureSX1 qL_PAR','-djpeg','-r300')
