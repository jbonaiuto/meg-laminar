function contrasts=create_contrast_structure()

contrasts=[];
% Resp-locked Beta rebound:
contrasts(1).zero_event='resp';
contrasts(1).foi=[15 30];
contrasts(1).invwoi=[-750 1500];
contrasts(1).comparison_name='resp_beta_rebound';
contrasts(1).comparison_woi=[500 1000];
contrasts(1).baseline_woi=[-250 250];
contrasts(1).hemisphere='left';
contrasts(1).region='lh_motor';
contrasts(1).direction='positive';

% Resp-locked gamma MRGS
contrasts(2).zero_event='resp';
contrasts(2).foi=[60 90];
contrasts(2).invwoi=[-1500 500];
contrasts(2).comparison_name='resp_mrgs';
contrasts(2).comparison_woi=[-100 200];
contrasts(2).baseline_woi=[-1500 -1000];
contrasts(2).hemisphere='left';
contrasts(2).region='lh_motor';
contrasts(2).direction='positive';

% Dots-aligned visual alpha
contrasts(3).zero_event='instr';
contrasts(3).foi=[7 13];
contrasts(3).invwoi=[-3500 0];
contrasts(3).comparison_name='dots_alpha';
contrasts(3).comparison_woi=[-2500 -500];
%contrasts(3).baseline_woi=[-3000 -2500];
contrasts(3).baseline_woi=[-3500 -3000];
contrasts(3).hemisphere='';
contrasts(3).region='visual';
contrasts(3).direction='negative';

% Dots-aligned visual gamma
contrasts(4).zero_event='instr';
contrasts(4).foi=[60 90];
contrasts(4).invwoi=[-3500 -1500];
contrasts(4).comparison_name='dots_gamma';
contrasts(4).comparison_woi=[-2250 -2000];
contrasts(4).baseline_woi=[-3000 -2750];
contrasts(4).hemisphere='';
contrasts(4).region='visual';
contrasts(4).direction='positive';

% Dots-aligned beta ERD
contrasts(5).zero_event='instr';
contrasts(5).foi=[15 30];
contrasts(5).invwoi=[-3500 0];
contrasts(5).comparison_name='dots_beta_erd';
contrasts(5).comparison_woi=[-2500 -500];
contrasts(5).baseline_woi=[-3000 -2500];
contrasts(5).hemisphere='left';
contrasts(5).region='lh_motor';
contrasts(5).direction='negative';

% Instructed-aligned visual gamma
contrasts(6).zero_event='instr';
contrasts(6).foi=[60 90];
contrasts(6).invwoi=[-750 500];
contrasts(6).comparison_name='instr_gamma';
contrasts(6).comparison_woi=[100 400];
contrasts(6).baseline_woi=[-500 -100];
contrasts(6).hemisphere='';
contrasts(6).region='visual';
contrasts(6).direction='positive';