%% RF model
clc; clear; close all;

 load CB6F_simulated_dataset % produce by the CB6F model with the following inputs, see the Johnson et al., 2021 for code
% %          1PAR      2Tleaf  3Cm       4fAPAR  5Vmax_CB6F  6Vmax_RUB  7Ku-50%  8fAPAR_PSII-(Pedrós2010 PSI/PSII = 1-2.4)
% VarMin=[1           10           1             0.1              87.5             1                       1                  0.5] ;
% VarMax=[2400   50          1000      0.9              262.5          200                   3                0.705];

input = [YA(:,1), A(:,2:3), A(:,5:end)];
output  =YA(:,2);
nTree = 50;%树的tree number
Factor = TreeBagger(nTree, input(1:0.7*size(input,1),:), output(1:0.7*size(input,1),:),'Method','regression','NumPredictorsToSample',122,'OOBpredictorImportance','on','MinLeafSize',5);%

[Predict_label,~] = predict(Factor, input(0.7*size(input,1):end,:));
cct=corrcoef(output(0.7*size(input,1):end,:),Predict_label');
cct=cct(2,1);
save Factor_SIF2.mat Factor
%% set the range for parameters


myFunction=@(x) randomforest(Factor, x);% target function, or a black box
D=7;% the number of parameter
nPop=1e06;% sample size
%          1SIF           2Tleaf   3Cm       5Vmax_CB6F  6Vmax_RUB  7Ku-50%  8fAPAR_PSII-(Pedrós2010 PSI/PSII = 1-2.4)
VarMin=[0.001      10          1               87.5                 1                         1              0.1] ;%min
VarMax=[30          50         1000         262.5              200                   3              0.9];%max of parameter

% GSA
numnPop=numel(nPop);
SAll=zeros(numnPop,D+1);%
STAll=zeros(numnPop,D+1);
i=numnPop;

[S,ST]=sobol(D,nPop(i),VarMin,VarMax,myFunction);

SAll(i,1:D)=S';
SAll(i,D+1)=sum(SAll(i,1:D));
STAll(i,1:D)=ST';
STAll(i,D+1)=sum(STAll(i,1:D));

load train
sound(y,Fs)

%% figure - bar

%load SIF_STi.mat

figure('visible','on');
set(gcf,'Position',[50 50 1000 600]);

x = 1:7;
y = ST./sum(ST);
x_label = {'SIF_F_U_L_L_\__P_S_I_I ', '{\itT}_L_e_a_f','{\itC}_m','{\itV}_M_A_X_\__C_B_6_f','{\itV}_M_A_X_\__R_u_b','{\itK}_u','{\it\beta}_2'};

bar(x, 100*y,0.5)
set(gca,'xticklabel',x_label,'Fontname','times new Roman','FontWeight','bold','FontSize',15)
set(gca, 'XTickLabel',x_label)
ylabel('{\itS}_T_i (%)','Color','k','Fontname','times new Roman','FontWeight','bold','FontSize',18)
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',20);
set(gca,'LooseInset',get(gca,'TightInset'));
box on
grid  on
print('GSA_SlF_others','-djpeg','-r300')

