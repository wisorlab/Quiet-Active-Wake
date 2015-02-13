% dynamic range script

% directories
BA10secdir = 'D:\mrempe\autoscore_and_epoch_length_study_data\10-Second_epochs\BA\';
BL10secdir = 'D:\mrempe\autoscore_and_epoch_length_study_data\10-Second_epochs\BL\';


kappa10secALL       = [kappaBA10sec kappaBL10sec];
globalagree10secALL = [globalagreeBA10sec globalagreeBL10sec];




% run PROCESSLBATCHMODE on both strains using delta to compute the dynamic range of delta power
[signal_data,state_data,residual,best_S,UppA,LowA,dynamic_rangeBAdelta,Timer,Taui,Taud]=PROCESSLBATCHMODE(BA10secdir,'delta2','NelderMead');
[signal_data,state_data,residual,best_S,UppA,LowA,dynamic_rangeBLdelta,Timer,Taui,Taud]=PROCESSLBATCHMODE(BL10secdir,'delta2','NelderMead');

% run PROCESSLBATCHMODE on both strain using lactate to compute the dynamic range of lactate 
[signal_data,state_data,residual,best_S,UppA,LowA,dynamic_rangeBAlactate,Timer,Taui,Taud]=PROCESSLBATCHMODE(BA10secdir,'lactate','NelderMead');
[signal_data,state_data,residual,best_S,UppA,LowA,dynamic_rangeBLlactate,Timer,Taui,Taud]=PROCESSLBATCHMODE(BL10secdir,'lactate','NelderMead');

% now concatenate dynamic range vectors
dynamic_range_delta   = [dynamic_rangeBAdelta dynamic_rangeBLdelta];
dynamic_range_lactate = [dynamic_rangeBAlactate dynamic_rangeBLlactate]

save dynamic_range_figure_data.mat kappa10secALL globalagree10secALL dynamic_range_delta dynamic_range_lactate 

figure 
subplot(1,2,1)
plot(dynamic_range_delta,kappa10secALL,'ko')
hold on 
plot(dynamic_range_delta,globalagree10secALL,'ro')
hold off
subplot(1,2,2)
plot(dynamic_range_lactate,kappa10secALL,'ko')
hold on
plot(dynamic_range_lactate,globalagree10secALL,'ro')
hold off
