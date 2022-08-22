clear;
ANM_BN={'Beta';'Gamma'};ANM_Ctrl={'Tien'};
CsvFiles = dir(['B_' '*.csv']);
Data_AllANM=[];
for i = 1 : length(CsvFiles)
    temp = readtable(CsvFiles(i).name);
    if find(strcmp(ANM_BN,temp.Subject(i)))
        temp.Stn=repelem(1,height(temp))';
    elseif find(strcmp(ANM_Ctrl,temp.Subject(i)))
        temp.Stn=repelem(0,height(temp))';
    else
        temp.Stn=repelem(2,height(temp))';
    end
    Data_AllANM = [Data_AllANM;temp];
end
clear temp

ANM_all = unique(Data_AllANM.Subject);
FPtype = [0.5 1 1.5];%unique(Data_AllANM.FP(Data_AllANM.Stage==1));
%% calculate
tablenames = {'ANM','Strain', 'Task', 'Treatment',...
    'nTrial','Accuracy','RT_median','WaitDur_median','ST_median',...
    'PerfTotal','PerfPort1','PerfPort2','Perf_shortFP','Perf_midFP','Perf_longFP',...
    'RT_Port','RT_FP','WaitDur_Port','WaitDur_FP'};
ATF_ANM=[];ATF_Stn=[];ATF_Task=[];ATF_Treatment=[];ATF_nTrial=[];ATF_Acc=[];ATF_RT=[];ATF_WaitDur=[];ATF_ST=[];
PerfTotal=[]; PerfPort1=[];PerfPort2=[];Perf_shortFP=[];Perf_midFP=[];Perf_longFP=[];
RT_Port=[];RT_FP=[];WaitDur_Port=[];WaitDur_FP=[];

for i = 1: length(ANM_all)  
    Temp=Data_AllANM(strcmp(Data_AllANM.Subject(:),ANM_all(i)) & Data_AllANM.Stage==1 & Data_AllANM.Treatment>-1,:);
    ATF_ANM=[ATF_ANM; Temp.Subject(1);Temp.Subject(1)];
    if  find(strcmp(ANM_BN, ANM_all(i)))
        ATF_Stn=[ATF_Stn;1;1]; % BN =1;
    elseif find(strcmp(ANM_Ctrl, ANM_all(i)))
        ATF_Stn=[ATF_Stn;0;0];
    else
        ATF_Stn=[ATF_Stn;2;2]; % LE=2;
    end
    ATF_Task=[ATF_Task; Temp.Task(1); Temp.Task(1)];
    ATF_Treatment=[ATF_Treatment; 0;1];% Saline;DCZ
    ATF_nTrial=[ATF_nTrial; sum(Temp.Treatment==0); sum(Temp.Treatment>0);];
    
    ATF_Acc=[ATF_Acc; sum(Temp.OutcomeCode==1 & Temp.Treatment==0)/sum(Temp.Treatment==0); sum(Temp.OutcomeCode==1 & Temp.Treatment>0)/sum(Temp.Treatment>0)];
    
    ATF_RT=[ATF_RT; compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.Treatment==0),'quartile')'; compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.Treatment==1),'quartile')']; % no warmup; only correct
    
    ATF_WaitDur=[ATF_WaitDur;median(Temp.WaitDur(Temp.Treatment==0),'omitnan'); median(Temp.WaitDur(Temp.Treatment>0),'omitnan')];
    ATF_ST=[ATF_ST; median(rmoutliers(Temp.ST(Temp.Treatment==0)),'omitnan'); median(rmoutliers(Temp.ST(Temp.Treatment>0)),'omitnan')];

    PerfTotal=[PerfTotal; [sum(Temp.Treatment==0), sum(Temp.OutcomeCode==1 & Temp.Treatment==0),sum(Temp.OutcomeCode<0 & Temp.Treatment==0),sum(Temp.OutcomeCode==0 & Temp.Treatment==0)];
                          [sum(Temp.Treatment>0), sum(Temp.OutcomeCode==1 & Temp.Treatment>0),sum(Temp.OutcomeCode<0 & Temp.Treatment>0),sum(Temp.OutcomeCode==0 & Temp.Treatment>0)]  ]; 
    PerfPort1=[PerfPort1; [sum(Temp.CorPort==1 & Temp.Treatment==0),sum(Temp.OutcomeCode(Temp.CorPort==1 & Temp.Treatment==0)==1),sum(Temp.OutcomeCode(Temp.CorPort==1 & Temp.Treatment==0)<0),sum(Temp.OutcomeCode(Temp.CorPort==1 & Temp.Treatment==0))==0];
                          [sum(Temp.CorPort==1 & Temp.Treatment>0),sum(Temp.OutcomeCode(Temp.CorPort==1 & Temp.Treatment>0)==1),sum(Temp.OutcomeCode(Temp.CorPort==1 & Temp.Treatment>0)<0),sum(Temp.OutcomeCode(Temp.CorPort==1 & Temp.Treatment>0)==0)]];
    PerfPort2=[PerfPort2; [sum(Temp.CorPort==2 & Temp.Treatment==0),sum(Temp.OutcomeCode(Temp.CorPort==2 & Temp.Treatment==0)==1),sum(Temp.OutcomeCode(Temp.CorPort==2 & Temp.Treatment==0)<0),sum(Temp.OutcomeCode(Temp.CorPort==2 & Temp.Treatment==0)==0)];
                          [sum(Temp.CorPort==2 & Temp.Treatment>0),sum(Temp.OutcomeCode(Temp.CorPort==2 & Temp.Treatment>0)==1),sum(Temp.OutcomeCode(Temp.CorPort==2 & Temp.Treatment>0)<0),sum(Temp.OutcomeCode(Temp.CorPort==2 & Temp.Treatment>0)==0)]];
    
    Perf_shortFP=[Perf_shortFP; [sum(Temp.FP==0.5 & Temp.Treatment==0),sum(Temp.OutcomeCode(Temp.FP==0.5 & Temp.Treatment==0)==1),sum(Temp.OutcomeCode(Temp.FP==0.5 & Temp.Treatment==0)<0),sum(Temp.OutcomeCode(Temp.FP==0.5 & Temp.Treatment==0)==0)];
                                [sum(Temp.FP==0.5 & Temp.Treatment>0),sum(Temp.OutcomeCode(Temp.FP==0.5 & Temp.Treatment>0)==1),sum(Temp.OutcomeCode(Temp.FP==0.5 & Temp.Treatment>0)<0),sum(Temp.OutcomeCode(Temp.FP==0.5 & Temp.Treatment>0)==0)]];
    Perf_midFP=[Perf_midFP; [sum(Temp.FP==1 & Temp.Treatment==0),sum(Temp.OutcomeCode(Temp.FP==1 & Temp.Treatment==0)==1),sum(Temp.OutcomeCode(Temp.FP==1 & Temp.Treatment==0)<0),sum(Temp.OutcomeCode(Temp.FP==1 & Temp.Treatment==0)==0)];
                            [sum(Temp.FP==1 & Temp.Treatment>0),sum(Temp.OutcomeCode(Temp.FP==1 & Temp.Treatment>0)==1),sum(Temp.OutcomeCode(Temp.FP==1 & Temp.Treatment>0)<0),sum(Temp.OutcomeCode(Temp.FP==1 & Temp.Treatment>0)==0)]];
    Perf_longFP=[Perf_longFP; [sum(Temp.FP==1.5 & Temp.Treatment==0),sum(Temp.OutcomeCode(Temp.FP==1.5 & Temp.Treatment==0)==1),sum(Temp.OutcomeCode(Temp.FP==1.5 & Temp.Treatment==0)<0),sum(Temp.OutcomeCode(Temp.FP==1.5 & Temp.Treatment==0)==0)];
                              [sum(Temp.FP==1.5 & Temp.Treatment>0),sum(Temp.OutcomeCode(Temp.FP==1.5 & Temp.Treatment>0)==1),sum(Temp.OutcomeCode(Temp.FP==1.5 & Temp.Treatment>0)<0),sum(Temp.OutcomeCode(Temp.FP==1.5 & Temp.Treatment>0)==0)]];
    
    RT_Port=[RT_Port;compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.CorPort==1 & Temp.Treatment==0),'quartile')', compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.CorPort==2 & Temp.Treatment==0),'quartile')';
                    compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.CorPort==1 & Temp.Treatment>0),'quartile')', compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.CorPort==2 & Temp.Treatment>0),'quartile')'];
    RT_FP=[RT_FP;compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.FP==0.5 & Temp.Treatment==0),'quartile')', compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.FP==1 & Temp.Treatment==0),'quartile')', compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.FP==1.5 & Temp.Treatment==0),'quartile')';
                 compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.FP==0.5 & Temp.Treatment>0),'quartile')', compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.FP==1 & Temp.Treatment>0),'quartile')', compute_stat_summary(Temp.RT(Temp.OutcomeCode==1 & Temp.FP==1.5 & Temp.Treatment>0),'quartile')'];
    
    WaitDur_Port=[WaitDur_Port;median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.CorPort==1 & Temp.Treatment==0),'omitnan'), median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.CorPort==2 & Temp.Treatment==0),'omitnan');
                               median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.CorPort==1 & Temp.Treatment>0),'omitnan'), median(Temp.RT(Temp.OutcomeCode==1 & Temp.CorPort==2 & Temp.Treatment>0),'omitnan')];
    WaitDur_FP=[WaitDur_FP;median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.FP==0.5 & Temp.Treatment==0),'omitnan'), median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.FP==1 & Temp.Treatment==0),'omitnan'), median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.FP==1.5 & Temp.Treatment==0),'omitnan');
                           median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.FP==0.5 & Temp.Treatment>0),'omitnan'), median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.FP==1 & Temp.Treatment>0),'omitnan'), median(Temp.WaitDur(Temp.OutcomeCode==1 & Temp.FP==1.5 & Temp.Treatment>0),'omitnan')];

end
clear Temp
ATF=table(ATF_ANM,ATF_Stn,ATF_Task,ATF_Treatment, ...
    ATF_nTrial,ATF_Acc,ATF_RT,ATF_WaitDur,ATF_ST, ...
    PerfTotal, PerfPort1,PerfPort2,Perf_shortFP,Perf_midFP,Perf_longFP,...
    RT_Port,RT_FP,WaitDur_Port,WaitDur_FP,...
    'VariableNames', tablenames);%Perf_shortFP,Perf_midFP,Perf_longFP,ATF_WT,
savename = ['GPS_',cell2mat(ATF_Task(1)), '_Training_Information'];
save(savename, 'ATF');
writetable(ATF, [savename '.csv']);
clear ATF_ANM ATF_Stn ATF_Date ATF_Task ATF_WaitDur ATF_Treatment ATF_nTrial ATF_ToC ATF_Acc ATF_Pre ATF_RT ATF_ST ATF_Bpodfiles ATF_Csvfiles PerfTotal PerfPort1 PerfPort2 Perf_shortFP Perf_midFP Perf_longFP RT_Port RT_FP WaitDur_Port WaitDur_FP

%% Add color
Color = [colororder;[0 1 0]*0.75];cSFP= [0.5 0.5 0.5];cLFP=[0.929,0.49,0.192];cMFP=mean([cSFP;cLFP]);
cAcc=[85 225 0]/255;cPre=[1 0 0];cLate=[140 140 140]/255;cErr=[70 130 180]/255;cDCZ=[253, 166, 0]/255;cBlue=[0 184 255]/255;
cOutcome=[cAcc;cPre;cErr];
cANM=[Color(1,:);Color(3,:);Color(5,:)];cBN=Color(1,:);cLE=Color(3,:);cCtrl=Color(5,:);

%% Group Performance
figure()
set(gcf,'unit', 'centimeters', 'position',[1 1 13 8], 'paperpositionmode', 'auto' );
% Performance
% FP=0.5s
axes('Units', 'centimeter', 'Position', [1 1 2.5 2.5],'FontSize',6); hold on;
ylabel('Performace(%)');ylim([0 100]);
xlim([0 3]);xticks([0:3]);xticklabels({'','Saline','DCZ'});xtickangle(0);
perf_FP=[ATF.Perf_shortFP(ATF.Treatment==0,2)./ATF.Perf_shortFP(ATF.Treatment==0,1),...
        ATF.Perf_shortFP(ATF.Treatment==0,3)./ATF.Perf_shortFP(ATF.Treatment==0,1),...
        ATF.Perf_shortFP(ATF.Treatment==0,4)./ATF.Perf_shortFP(ATF.Treatment==0,1),...        
        ATF.Perf_shortFP(ATF.Treatment==1,2)./ATF.Perf_shortFP(ATF.Treatment==1,1),... 
        ATF.Perf_shortFP(ATF.Treatment==1,3)./ATF.Perf_shortFP(ATF.Treatment==1,1),...
        ATF.Perf_shortFP(ATF.Treatment==1,4)./ATF.Perf_shortFP(ATF.Treatment==1,1)];
perf_FP_Ebar=[];
for  i =1 : length(perf_FP)
    perf_FP_Ebar(:,i)=compute_stat_summary(perf_FP(1:4,i),'mean-sem');
end
perf_FP_Ebar(2,:)=perf_FP_Ebar(1,:)-perf_FP_Ebar(2,:);
perf_FP_Ebar(3,:)=perf_FP_Ebar(3,:)-perf_FP_Ebar(1,:);
perf_FP_Ebar=perf_FP_Ebar.*100;
b=bar([perf_FP_Ebar(1,1:3);perf_FP_Ebar(1,4:6)],1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
for k = 1:3
    b(k).FaceColor = cOutcome(k,:); 
    b(k).EdgeColor = 'k';   
end
hold on;
errorbar([0.775 1.775], [perf_FP_Ebar(1,1);perf_FP_Ebar(1,4)], [perf_FP_Ebar(2,1);perf_FP_Ebar(2,4)],[perf_FP_Ebar(3,1);perf_FP_Ebar(3,4)], 'k', 'Linestyle', 'None');  hold on;
errorbar([1 2], [perf_FP_Ebar(1,2);perf_FP_Ebar(1,5)], [perf_FP_Ebar(2,2);perf_FP_Ebar(2,5)],[perf_FP_Ebar(3,2);perf_FP_Ebar(3,5)],'k', 'Linestyle', 'None'); hold on;
errorbar([1.225 2.225], [perf_FP_Ebar(1,3);perf_FP_Ebar(1,6)], [perf_FP_Ebar(2,3);perf_FP_Ebar(2,6)],[perf_FP_Ebar(3,3);perf_FP_Ebar(3,6)], 'k', 'Linestyle', 'None'); 
text(-1,100,'B','FontSize',10)
hold off;

% FP=1.0s
axes('Units', 'centimeter', 'Position', [4.5 1 2.5 2.5],'FontSize',6); hold on;
ylabel('Performace(%)');ylim([0 100]);
xlim([0 3]);xticks([0:3]);xticklabels({'','Saline','DCZ'});xtickangle(0);
perf_FP=[ATF.Perf_midFP(ATF.Treatment==0,2)./ATF.Perf_midFP(ATF.Treatment==0,1),...
        ATF.Perf_midFP(ATF.Treatment==0,3)./ATF.Perf_midFP(ATF.Treatment==0,1),...
        ATF.Perf_midFP(ATF.Treatment==0,4)./ATF.Perf_midFP(ATF.Treatment==0,1),...        
        ATF.Perf_midFP(ATF.Treatment==1,2)./ATF.Perf_midFP(ATF.Treatment==1,1),... 
        ATF.Perf_midFP(ATF.Treatment==1,3)./ATF.Perf_midFP(ATF.Treatment==1,1),...
        ATF.Perf_midFP(ATF.Treatment==1,4)./ATF.Perf_midFP(ATF.Treatment==1,1)];
perf_FP_Ebar=[];
for  i =1 : length(perf_FP)
    perf_FP_Ebar(:,i)=compute_stat_summary(perf_FP(1:4,i),'mean-sem');
end
perf_FP_Ebar(2,:)=perf_FP_Ebar(1,:)-perf_FP_Ebar(2,:);
perf_FP_Ebar(3,:)=perf_FP_Ebar(3,:)-perf_FP_Ebar(1,:);
perf_FP_Ebar=perf_FP_Ebar.*100;
b=bar([perf_FP_Ebar(1,1:3);perf_FP_Ebar(1,4:6)],1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
for k = 1:3
    b(k).FaceColor = cOutcome(k,:); 
    b(k).EdgeColor = 'k';   
end
hold on;
errorbar([0.775 1.775], [perf_FP_Ebar(1,1);perf_FP_Ebar(1,4)], [perf_FP_Ebar(2,1);perf_FP_Ebar(2,4)],[perf_FP_Ebar(3,1);perf_FP_Ebar(3,4)], 'k', 'Linestyle', 'None');  hold on;
errorbar([1 2], [perf_FP_Ebar(1,2);perf_FP_Ebar(1,5)], [perf_FP_Ebar(2,2);perf_FP_Ebar(2,5)],[perf_FP_Ebar(3,2);perf_FP_Ebar(3,5)],'k', 'Linestyle', 'None'); hold on;
errorbar([1.225 2.225], [perf_FP_Ebar(1,3);perf_FP_Ebar(1,6)], [perf_FP_Ebar(2,3);perf_FP_Ebar(2,6)],[perf_FP_Ebar(3,3);perf_FP_Ebar(3,6)], 'k', 'Linestyle', 'None'); 
hold off;

% FP=1.5s
axes('Units', 'centimeter', 'Position', [8 1 2.5 2.5],'FontSize',6); hold on;
ylabel('Performace(%)');ylim([0 100]);
xlim([0 3]);xticks([0:3]);xticklabels({'','Saline','DCZ'});xtickangle(0);
perf_FP=[ATF.Perf_longFP(ATF.Treatment==0,2)./ATF.Perf_longFP(ATF.Treatment==0,1),...
        ATF.Perf_longFP(ATF.Treatment==0,3)./ATF.Perf_longFP(ATF.Treatment==0,1),...
        ATF.Perf_longFP(ATF.Treatment==0,4)./ATF.Perf_longFP(ATF.Treatment==0,1),...        
        ATF.Perf_longFP(ATF.Treatment==1,2)./ATF.Perf_longFP(ATF.Treatment==1,1),... 
        ATF.Perf_longFP(ATF.Treatment==1,3)./ATF.Perf_longFP(ATF.Treatment==1,1),...
        ATF.Perf_longFP(ATF.Treatment==1,4)./ATF.Perf_longFP(ATF.Treatment==1,1)];
perf_FP_Ebar=[];
for  i =1 : length(perf_FP)
    perf_FP_Ebar(:,i)=compute_stat_summary(perf_FP(1:4,i),'mean-sem');
end
perf_FP_Ebar(2,:)=perf_FP_Ebar(1,:)-perf_FP_Ebar(2,:);
perf_FP_Ebar(3,:)=perf_FP_Ebar(3,:)-perf_FP_Ebar(1,:);
perf_FP_Ebar=perf_FP_Ebar.*100;
b=bar([perf_FP_Ebar(1,1:3);perf_FP_Ebar(1,4:6)],1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
for k = 1:3
    b(k).FaceColor = cOutcome(k,:); 
    b(k).EdgeColor = 'k';   
end
hold on;
errorbar([0.775 1.775], [perf_FP_Ebar(1,1);perf_FP_Ebar(1,4)], [perf_FP_Ebar(2,1);perf_FP_Ebar(2,4)],[perf_FP_Ebar(3,1);perf_FP_Ebar(3,4)], 'k', 'Linestyle', 'None');  hold on;
errorbar([1 2], [perf_FP_Ebar(1,2);perf_FP_Ebar(1,5)], [perf_FP_Ebar(2,2);perf_FP_Ebar(2,5)],[perf_FP_Ebar(3,2);perf_FP_Ebar(3,5)],'k', 'Linestyle', 'None'); hold on;
errorbar([1.225 2.225], [perf_FP_Ebar(1,3);perf_FP_Ebar(1,6)], [perf_FP_Ebar(2,3);perf_FP_Ebar(2,6)],[perf_FP_Ebar(3,3);perf_FP_Ebar(3,6)], 'k', 'Linestyle', 'None'); 
hold off;

axes('Units', 'centimeter', 'Position', [11 1 2.5 2.5],'FontSize',6,'XLim',[0 10],'ylim',[0,10]); hold on;
plot([0,2],[9 9],'Color',cOutcome(1,:),'LineWidth',4);
plot([0,2],[8 8],'Color',cOutcome(2,:),'LineWidth',4);
plot([0,2],[7 7],'Color',cOutcome(3,:),'LineWidth',4);
text(3,9,'Correct','FontSize',6)
text(3,8,'Premature','FontSize',6)
text(3,7,'ChoiceError','FontSize',6)
axis off;

%% Group WaitDur
timebins = 0:0.05:4;
DataUsed=Data_AllANM(Data_AllANM.Stage==1 & Data_AllANM.Stn>0,:);
WaitDur=DataUsed.WaitDur;
% figure()
% set(gcf,'unit', 'centimeters', 'position',[1 1 13 8], 'paperpositionmode', 'auto' );
for j=1:length(FPtype)
    axes('Units', 'centimeter', 'Position', [1+3.5*(j-1) 4.5 2.5 2.5],'FontSize',6); hold on;
    xlim([0,3]);xlabel('Wait duration (s)');xticks([0:0.5:4]);xtickangle(0);
    ylim([0,1]);ylabel('CDF');
    title(['FP=',sprintf('%.1f',FPtype(j)),'s']);
    Acdf_Saline=[];Acdf_DCZ=[];cdf_Saline=[];cdf_DCZ=[];
    for i =1:length(ANM_all)-1
    
        Acdf_Saline(i,:)= ksdensity([WaitDur(strcmp(DataUsed.Subject,ANM_all(i)) & DataUsed.FP==FPtype(j) & DataUsed.Treatment==0) ], timebins,'Support','positive',...
            'Function','cdf','Bandwidth',0.05);
        Acdf_DCZ(i,:) = ksdensity([WaitDur(strcmp(DataUsed.Subject,ANM_all(i)) & DataUsed.FP==FPtype(j) & DataUsed.Treatment==1)], timebins,'Support','positive',...
            'Function','cdf','Bandwidth',0.05);
    end
    for i =1 :size(Acdf_DCZ,2)
        cdf_Saline(:,i)=compute_stat_summary(Acdf_Saline(:,i),'mean-sem');
        cdf_DCZ(:,i)=compute_stat_summary(Acdf_DCZ(:,i),'mean-sem');
    end
    plot(timebins, cdf_Saline(1,:), '-', 'Color',  [0 0 0], 'linewidth', 1);hold on;
    plot(timebins, cdf_DCZ(1,:), '-', 'Color',  cDCZ, 'linewidth', 1);hold on;
    plotshaded(timebins,cdf_Saline(2:3,:),[0 0 0]);hold on;
    plotshaded(timebins,cdf_DCZ(2:3,:),cDCZ);hold on;
    plot([FPtype(j),FPtype(j)],[0 1],'--','Color','k','LineWidth',.7);hold on;
    if j==1
        text(-1,1,'A','FontSize',10)
    end
    hold off;
end

axes('Units', 'centimeter', 'Position', [11 4.5 2.5 2.5],'FontSize',6,'XLim',[0 10],'ylim',[0,10]); hold on;
plot([0,2],[9 9],'Color','k','LineWidth',1);
plot([0,2],[8 8],'Color',cDCZ,'LineWidth',1);

text(3,9,'Saline','FontSize',6)
text(3,8,'DCZ','FontSize',6)

axis off;


if ~isfolder('plot')
    mkdir('plot');
end
savename = ['GPS_',cell2mat(ATF.Task(1)),'_GroupPerf'];
savename = fullfile(pwd, 'plot', savename);
saveas(gcf, savename, 'png')
saveas(gcf, savename, 'fig')
%% Sample Performance+Waitdur
ANM_Sample={'Beta'};

figure()
set(gcf,'unit', 'centimeters', 'position',[1 1 13 8], 'paperpositionmode', 'auto' );
% Performance
% FP=0.5s
axes('Units', 'centimeter', 'Position', [1 1 2.5 2.5],'FontSize',6); hold on;
ylabel('Performace(%)');ylim([0 100]);
xlim([0 3]);xticks([0:3]);xticklabels({'','Saline','DCZ'});xtickangle(0);

perf_FP=[ATF.Perf_shortFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==0,2:4)./ATF.Perf_shortFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==0,1);    
        ATF.Perf_shortFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==1,2:4)./ATF.Perf_shortFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==1,1)];

perf_FP=perf_FP.*100;
b=bar(perf_FP,1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
for k = 1:3
    b(k).FaceColor = cOutcome(k,:); 
    b(k).EdgeColor = 'k';   
end

text(-1,100,'B','FontSize',10)

hold off;

% FP=1.0s
axes('Units', 'centimeter', 'Position', [4.5 1 2.5 2.5],'FontSize',6); hold on;
ylabel('Performace(%)');ylim([0 100]);
xlim([0 3]);xticks([0:3]);xticklabels({'','Saline','DCZ'});xtickangle(0);
perf_FP=[ATF.Perf_midFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==0,2:4)./ATF.Perf_midFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==0,1);    
        ATF.Perf_midFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==1,2:4)./ATF.Perf_midFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==1,1)];

perf_FP=perf_FP.*100;
b=bar(perf_FP,1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
for k = 1:3
    b(k).FaceColor = cOutcome(k,:); 
    b(k).EdgeColor = 'k';   
end
hold off;

% FP=1.5s
axes('Units', 'centimeter', 'Position', [8 1 2.5 2.5],'FontSize',6); hold on;
ylabel('Performace(%)');ylim([0 100]);
xlim([0 3]);xticks([0:3]);xticklabels({'','Saline','DCZ'});xtickangle(0);
perf_FP=[ATF.Perf_longFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==0,2:4)./ATF.Perf_longFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==0,1);    
        ATF.Perf_longFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==1,2:4)./ATF.Perf_longFP(strcmp(ATF.ANM,ANM_Sample) & ATF.Treatment==1,1)];

perf_FP=perf_FP.*100;
b=bar(perf_FP,1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
for k = 1:3
    b(k).FaceColor = cOutcome(k,:); 
    b(k).EdgeColor = 'k';   
end
hold off;

axes('Units', 'centimeter', 'Position', [11 1 2.5 2.5],'FontSize',6,'XLim',[0 10],'ylim',[0,10]); hold on;
plot([0,2],[9 9],'Color',cOutcome(1,:),'LineWidth',4);
plot([0,2],[8 8],'Color',cOutcome(2,:),'LineWidth',4);
plot([0,2],[7 7],'Color',cOutcome(3,:),'LineWidth',4);
text(3,9,'Correct','FontSize',6)
text(3,8,'Premature','FontSize',6)
text(3,7,'ChoiceError','FontSize',6)
axis off;

% WaitDur
timebins = 0:0.05:4;ylabel('CDF');
Data_Sample=Data_AllANM(Data_AllANM.Stage==1 & strcmp(Data_AllANM.Subject,ANM_Sample),:);
WaitDur=Data_Sample.WaitDur;
% figure()
% set(gcf,'unit', 'centimeters', 'position',[1 1 13 8], 'paperpositionmode', 'auto' );
for j=1:length(FPtype)
    axes('Units', 'centimeter', 'Position', [1+3.5*(j-1) 4.5 2.5 2.5],'FontSize',6); hold on;
    xlim([0,3]);xlabel('Wait duration (s)');xticks([0:0.5:4]);xtickangle(0);
    ylim([0,1]);ylabel('CDF');
    title(['FP=',sprintf('%.1f',FPtype(j)),'s']);    
    cdf_Saline= ksdensity([WaitDur(Data_Sample.FP==FPtype(j) & Data_Sample.Treatment==0) ], timebins,'Support','positive',...
        'Function','cdf','Bandwidth',0.05);
    cdf_DCZ= ksdensity([WaitDur(Data_Sample.FP==FPtype(j) & Data_Sample.Treatment==1)], timebins,'Support','positive',...
        'Function','cdf','Bandwidth',0.05);
    plot(timebins, cdf_Saline(1,:), '-', 'Color',  [0 0 0], 'linewidth', 1);hold on;
    plot(timebins, cdf_DCZ(1,:), '-', 'Color',  cDCZ, 'linewidth', 1);hold on;
    plot([FPtype(j),FPtype(j)],[0 1],'--','Color','k','LineWidth',.7);hold on;
    hold off;
    if j==1
        text(-1,1,'A','FontSize',10)
    end
end

axes('Units', 'centimeter', 'Position', [11 4.5 2.5 2.5],'FontSize',6,'XLim',[0 10],'ylim',[0,10]); hold on;
plot([0,2],[9 9],'Color','k','LineWidth',1);
plot([0,2],[8 8],'Color',cDCZ,'LineWidth',1);

text(3,9,'Saline','FontSize',6)
text(3,8,'DCZ','FontSize',6)

axis off;


if ~isfolder('plot')
    mkdir('plot');
end
savename = ['GPS_',cell2mat(ATF.Task(1)),'_SamplePerf'];
savename = fullfile(pwd, 'plot', savename);
saveas(gcf, savename, 'png')
saveas(gcf, savename, 'fig')
%% ST
% Sample ST
Data_Sample=Data_AllANM(Data_AllANM.Stage==1 & strcmp(Data_AllANM.Subject,ANM_Sample),:);
DataST=Data_Sample(Data_Sample.Treatment>-1,:);
DataST = rmoutliers(DataST, 'median', 'DataVariables', 'ST');
TreatmentST = DataST.Treatment;TreatmentST(TreatmentST>0)=1;
ST = DataST.ST;

figure()
set(gcf, 'unit', 'centimeters', 'position',[18 2 17 5], 'paperpositionmode', 'auto' );
axes('Units', 'centimeter', 'Position', [1 1 3 3],'FontSize', 6, 'YScale', 'log'); hold on;
violinplot(ST, cellstr(string(TreatmentST)), 'ViolinColor', [.5 .5 .5; .5 .5 .5], ...
    'ShowMedian', false, 'ShowWhiskers', false, 'ShowBox', false, 'ShowData', false);
scatter(1, median(ST(TreatmentST==0), 'omitnan'), 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.3 .3 .3], 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.2);
scatter(2, median(ST(TreatmentST>0), 'omitnan'), 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.3 .3 .3], 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.2);
ylim([0.5 15]); xlim([.5 2.5]);
ylabel('ShuttleTime / s', 'FontSize', 8);
xticklabels({'Saline', 'DCZ'}); xlabel('Treatment', 'FontSize', 8);
p = anovan(ST(TreatmentST>-1), TreatmentST(TreatmentST>-1));
box off;
pval = sprintf('%.6f', p);
text(1,12,['p = ',pval],'FontSize', 6);
text(0,15,'A','FontSize',10)
title('Sample','FontSize', 8)

DataSTG=Data_AllANM(Data_AllANM.Stage==1 & Data_AllANM.Treatment>-1 & Data_AllANM.Stn>0,:);
DataSTG = rmoutliers(DataSTG, 'median', 'DataVariables', 'ST');
TreatmentST = DataSTG.Treatment;TreatmentST(TreatmentST>0)=1;
ST = DataSTG.ST;
axes('Units', 'centimeter', 'Position', [5 1 3 3],'FontSize', 6, 'YScale', 'log'); hold on;
violinplot(ST, cellstr(string(TreatmentST)), 'ViolinColor', [.5 .5 .5; .5 .5 .5], ...
    'ShowMedian', false, 'ShowWhiskers', false, 'ShowBox', false, 'ShowData', false);
for i =1: length(ANM_all)-1

    scatter(1, median(ST(strcmp(DataSTG.Subject,ANM_all(i)) & TreatmentST==0), 'omitnan'), 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.3 .3 .3], 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.2);
    scatter(2, median(ST(strcmp(DataSTG.Subject,ANM_all(i)) & TreatmentST>0), 'omitnan'), 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.3 .3 .3], 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.2);
end
ylim([0.5 15]); xlim([.5 2.5]);
ylabel('ShuttleTime / s', 'FontSize', 8);
xticklabels({'Saline', 'DCZ'}); xlabel('Treatment', 'FontSize', 8);
p = anovan(ST(TreatmentST>-1), TreatmentST(TreatmentST>-1));
box off;
pval = sprintf('%.6f', p);
text(1,12,['p = ',pval],'FontSize', 6);
text(0,15,'B','FontSize',10)
title('Group','FontSize', 8)

axes('Units', 'centimeter', 'Position', [9 1 3 3],'FontSize', 6, 'YScale', 'log'); hold on;
violinplot(ST(DataSTG.Stn==1), cellstr(string(TreatmentST(DataSTG.Stn==1))), 'ViolinColor', [.5 .5 .5; .5 .5 .5], ...
    'ShowMedian', false, 'ShowWhiskers', false, 'ShowBox', false, 'ShowData', false);
for i =1: 2

    scatter(1, median(ST(strcmp(DataSTG.Subject,ANM_all(i)) & TreatmentST==0), 'omitnan'), 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.3 .3 .3], 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.2);
    scatter(2, median(ST(strcmp(DataSTG.Subject,ANM_all(i)) & TreatmentST>0), 'omitnan'), 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.3 .3 .3], 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.2);
end
ylim([0.5 15]); xlim([.5 2.5]);
ylabel('ShuttleTime / s', 'FontSize', 8);
xticklabels({'Saline', 'DCZ'}); xlabel('Treatment', 'FontSize', 8);
p = anovan(ST(TreatmentST>-1), TreatmentST(TreatmentST>-1));
box off;
pval = sprintf('%.6f', p);
text(1,12,['p = ',pval],'FontSize', 6);
text(0,15,'C','FontSize',10)
title('BN','FontSize', 8)

axes('Units', 'centimeter', 'Position', [13 1 3 3],'FontSize', 6, 'YScale', 'log'); hold on;
violinplot(ST(DataSTG.Stn==2), cellstr(string(TreatmentST(DataSTG.Stn==2))), 'ViolinColor', [.5 .5 .5; .5 .5 .5], ...
    'ShowMedian', false, 'ShowWhiskers', false, 'ShowBox', false, 'ShowData', false);
for i =3:4

    scatter(1, median(ST(strcmp(DataSTG.Subject,ANM_all(i)) & TreatmentST==0), 'omitnan'), 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.3 .3 .3], 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.2);
    scatter(2, median(ST(strcmp(DataSTG.Subject,ANM_all(i)) & TreatmentST>0), 'omitnan'), 20, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.3 .3 .3], 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.2);
end
ylim([0.5 15]); xlim([.5 2.5]);
ylabel('ShuttleTime / s', 'FontSize', 8);
xticklabels({'Saline', 'DCZ'}); xlabel('Treatment', 'FontSize', 8);
p = anovan(ST(TreatmentST>-1), TreatmentST(TreatmentST>-1));
box off;
pval = sprintf('%.6f', p);
text(1,12,['p = ',pval],'FontSize', 6);
text(0,15,'D','FontSize',10)
title('LE','FontSize', 8)

if ~isfolder('plot')
    mkdir('plot');
end
savename = ['GPS_',cell2mat(ATF.Task(1)),'_ST'];
savename = fullfile(pwd, 'plot', savename);
saveas(gcf, savename, 'png')
saveas(gcf, savename, 'fig')
