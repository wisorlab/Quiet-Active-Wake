% script to make Active vs. Quiet wake lactate figure for Janne's paper
%
% UPDATE 6.26.15: 
% You can run this script after loading beta_paper_dataALL.mat 
% and comment out the lines that call PROCESSLBATCHMODE.m 
% 
% Or just run the script as is and just load beta_paper_data.mat
load beta_paper_data.mat


% File (and directory) to use:
directory='\\FS1\WisorData\Gronli\Stat_Analyses\NEW_incl_beta_gamma_Sleep_Report_Txt_Output_112514\Txt_files\Michael\Used_in_figure\';

% file = strcat(directory,'2012_S.txt');  % or whichever file you have in "directory"
% As of 3.5.2015, this script uses file 2012_S.txt (not 2116_T.txt)
%


% I'm assuming that PROCESSLBATCHMODE.m has already been run and state_data and signal_data have been set up for both the
% 3-state version of the model and the 5-state version. 

% Alternatively you can run them here:
% the 3-state model:
cd C:\Users\wisorlab\Documents\MATLAB\mrempe\Epoch-Based-Processing\Timed_Intervals\internal\used_in_beta_paper\
[signal_data,state_data,residual,scaled_L_3state,UppA,LowA,dynamic_range,Timer,Taui,TauD]=PROCESSLBATCHMODE(directory,'lactate','NelderMead');
%
% the 5-state model: (5 states, but only 2 tau values: one for rising during wake, REMS or AW, and one for falling during SWS or QW)
cd C:\Users\wisorlab\Documents\MATLAB\mrempe\Todays_code_testing\Quiet-Active-Wake\Only2Taus\Quiet-Active-Wake\
[signal_data,state_data,residual,best_S,UppA,LowA,dynamic_range,Timer,Taui,TauD]=PROCESSLBATCHMODE(directory,'lactate','NelderMead');

% axis_label_fontsize sets the size (in points) of the font labeling all axes
axis_label_fontsize = 12;


if iscell(signal_data)
	signal_data=cell2mat(signal_data);
end
if iscell(scaled_L_3state)
  scaled_L_3state=cell2mat(scaled_L_3state);
end
if iscell(state_data)
	state_data=cell2mat(state_data);
end
if iscell(best_S)
	best_S=cell2mat(best_S);
end
if iscell(UppA)
	UA=cell2mat(UppA);
end
if iscell(LowA)
	LA=cell2mat(LowA);
end

% Since the output from PROCESSLBATCHMODE.m is a scaled version of S, (the best fit of the model), unscale it.
scaled_L = best_S;
best_fit = (UA-LA).*scaled_L+LA;



dt=1/(60*60/10);
t=0:dt:dt*(length(signal_data)-1);
tS=t((4/2)*(60*60/10)+1:end-(4/2)*(60*60/10));

figure
% --- TOP PANEL ---------------------
subplot(3,6,1:6)
L_indices = (4/2)*(60*60/10)+1:length(t)-(4/2)*(60*60/10);

only_sleep_indices = find(state_data==1);
    only_wake_indices  = find(state_data==0);
    only_rem_indices   = find(state_data==2);
    only_quiet_wake_indices  = find(state_data==3);
    only_active_wake_indices = find(state_data==4);

only_sleep_indices_L = find(state_data(L_indices)==1);
    only_wake_indices_L  = find(state_data(L_indices)==0);
    only_rem_indices_L   = find(state_data(L_indices)==2);
    only_quiet_wake_indices_L  = find(state_data(L_indices)==3);
    only_active_wake_indices_L = find(state_data(L_indices)==4);

    scaled_lactate_data = ((UA-LA)-(UA-signal_data(L_indices)'))./(UA-LA);

    sleep_lactate       = signal_data(only_sleep_indices);
    wake_lactate        = signal_data(only_wake_indices);
    rem_lactate         = signal_data(only_rem_indices);
    quiet_wake_lactate  = signal_data(only_quiet_wake_indices);
    active_wake_lactate = signal_data(only_active_wake_indices);

    sleep_lactate_scaled = scaled_lactate_data(only_sleep_indices_L);
    wake_lactate_scaled  = scaled_lactate_data(only_wake_indices_L);
    rem_lactate_scaled   = scaled_lactate_data(only_rem_indices_L);
    quiet_wake_lactate_scaled  = scaled_lactate_data(only_quiet_wake_indices_L);
    active_wake_lactate_scaled = scaled_lactate_data(only_active_wake_indices_L);

    scatter(tS(only_wake_indices_L),wake_lactate_scaled,25,'d','MarkerEdgeColor',[0.5 0.5 0.5])
      
    hold on
    scatter(tS(only_sleep_indices_L),sleep_lactate_scaled,25,'k+')
    scatter(tS(only_rem_indices_L),rem_lactate_scaled,25,'x','MarkerEdgeColor',[0.75 0.75 0.75])
    % scatter(tS(only_quiet_wake_indices_L),quiet_wake_lactate_scaled,25,'bs')
    % scatter(tS(only_active_wake_indices_L),active_wake_lactate_scaled,25,'ro')
    
   % %tS=t(361:end-(60*60/epoch_length));
   %  %tS=t((4/2)*(60*60/10)+1:end-(4/2)*(60*60/10));
    
    plot(tS,scaled_L_3state,'k')
    % plot(tS,LA,'k--')
    % plot(tS,UA,'k--')
    hold off
    set(gca,'FontSize',axis_label_fontsize)
    axis([28 46 -0.1 1.2])
    ylabel('Lactate','FontSize',axis_label_fontsize)
    xlabel('TIME (H)','FontSize',axis_label_fontsize)
    %title(['Best fit of model to lactate data for file ' filename 'using ' num2str(epoch_length) '-second epochs'])
    one_legend=legend('Wake','SWS','REMS','Model','Location','EastOutside');
    legend boxoff
    set(one_legend,'FontSize',14);
    panel_label_text = text(28,1.15,'A');
    set(panel_label_text,'FontSize',18);

% --- SECOND PANEL -----------------
	subplot(3,6,7:12)
    %scatter(t(only_wake_indices),wake_lactate,25,'rd')
      
    % hold on
    %scatter(t(only_sleep_indices),sleep_lactate,25,'b+')
    %scatter(t(only_rem_indices),rem_lactate,25,'rx')
    scatter(tS(only_quiet_wake_indices_L),quiet_wake_lactate_scaled,25,'ks')
    hold on
    scatter(tS(only_active_wake_indices_L),active_wake_lactate_scaled,25,'o','MarkerEdgeColor',[0.5 0.5 0.5])
    
   %tS=t(361:end-(60*60/epoch_length));
    tS=t((4/2)*(60*60/10)+1:end-(4/2)*(60*60/10));
    
    plot(tS,scaled_L,'k')
    % plot(tS,LA,'k--')
    % plot(tS,UA,'k--')
    hold off
    set(gca,'FontSize',axis_label_fontsize)
    axis([28 46 -0.1 1.2])
    ylabel('Lactate','FontSize',axis_label_fontsize)
    xlabel('TIME (H)','FontSize',axis_label_fontsize)
    two_legend=legend('QW','AW','Model','Location','EastOutside');
    legend boxoff
    set(two_legend,'FontSize',14);
    panel_label_text = text(28,1.15,'B');
    set(panel_label_text,'FontSize',18);

% --- THIRD PANEL (BOTTOM LEFT)-----------------
subplot(3,6,13:14)
X = [ones(1,14); 2*ones(1,14)];
Y=[residuals_old_model_lactate; residuals_active_quiet_lactate];
H=line(X,Y);
axis([0.5 2.5 0 3.5])
set(H(:),'Color',[0 0 0])
set(H(:),'Marker','.')
set(H(:),'MarkerSize',10)
ylabel('Residuals','FontSize',axis_label_fontsize)
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'3-State\newlineModel';'5-State\newlineModel';})
panel_label_text = text(0.5,3.4,'C');
set(panel_label_text,'FontSize',18);


% --- FOURTH PANEL ----------------
subplot(3,6,16:17)
bar([0 1],[mean(residuals_old_model_lactate) mean(residuals_active_quiet_lactate)],0.3)
colormap(gray)
set(gca,'box','off')
axis([-0.8 1.9 0 1.4])
pval_text=text(0.8,1.2,'P=0.0273');
set(pval_text,'FontSize',12);
hold on 
E=[std(residuals_old_model_lactate)/sqrt(length(residuals_old_model_lactate)) std(residuals_active_quiet_lactate)/sqrt(length(residuals_active_quiet_lactate))];
errorbar([0 1],[mean(residuals_old_model_lactate) mean(residuals_active_quiet_lactate)],E,'Color','k','LineStyle','none');
hold off
ylabel('Residuals')
set(gca,'XTickLabel',{'3-State\newlineModel'; '5-State\newlineModel'})
panel_label_text = text(-0.8,1.35,'D');
set(panel_label_text,'FontSize',18);