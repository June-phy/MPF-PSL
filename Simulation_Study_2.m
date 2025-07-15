%------------------------ Header Comment -------------------------------%
% Authors:      Qianxue Shan, Ziqiang Yu, Weitian Chen
% Email:        qxshan@link.cuhk.edu.hk
% Version:      1.0
% Date:         2025-07-13
%
% Copyright (c) 2025 Qianxue Shan, Ziqiang Yu, Weitian Chen. All rights reserved.
%
% License and Usage Notice:
%   This code is provided strictly for academic and research purposes only.
%   Any commercial use, including but not limited to sale, redistribution,
%   or integration into proprietary software, is strictly prohibited without
%   explicit written permission from the authors.
%
%   Modification of this code, its header comments, or removal of this notice,
%   in whole or in part, is EXPRESSLY FORBIDDEN without prior written consent
%   from the authors.
%
%   By using, copying, or referencing this code, you agree to abide by these terms.
%   For any inquiries or requests, please contact the authors at the email above.
%
% Description:  
%   Demo script for Simulation Study 2 in the MPF-PSL paper.
%
%   This script investigates how the Relative Measurement Precision (RMP)
%   changes with increasing spin-lock duration.
%
% Functions:
%   - Given the spin-lock duration, calculate the RMP of Rmpfsl,pul under
%     different noise levels. 
%
% History: 
%   - First built on 2025-07-13.
%
% ----------------------------------------------------------------------%
%%
clc; clear; close all;

% load adibatic pulse
load('AFP_waveform.mat');
down_afp_rfa = down_afp_rfa.*42.58;
down_afp_time = down_afp_time.*1e-3;
diff_afp_time = diff(down_afp_time);
np_afp = length(diff_afp_time);

% set adibatic parameters
adi_set = 1;
dictionary_train_transient_pool_parameters;
amscale_fx = 400;
fmscale_fx = 4000;
adi_freq_fx = 4000;
hsdur_fx = 5*1.e-3;

MxA = zeros(2,nfo, 10);
MyA = MxA;
MzA = MxA;
MxMT = MxA;
MyMT = MxA;
MzMT = MxA;

MxA_conv = MxA;                
MyA_conv = MxA;                
MzA_conv = MxA;                
MzMT_conv = MxA;
MxMT_conv = MxA;
MyMT_conv = MxA;

use_Mini_t1rec = 1;
t1_recovery_time = 1.441;
TI = 0;
Mini=[0 0 0 0 0 0 0 0 0];

ntrain = 10;
% The duration of single spin-lock unit
tsl = 10*1.e-3;

for itog = 1:2
for ifo = fost:foed 
    b1_fx = sla(ifo); 
    fo_fx = freq_offset(ifo); 
    b1 = 1;
    b0 = 0;     
    wait_time2 = 50 *1.e-3;  
       
    %% initialize simulation parameters
    M0b=pb; 
    M0c=pc; 
    M0a=1-M0b-M0c; 
    kab=(M0b/M0a)*kba; 

    kac = kca*(M0c/M0a); 
    SL=sla(ifo)*2*pi; 

    Mall = Mini;                      
    start = 0;
    fin   = t1_recovery_time - TI; 
    M_t1rec = []; 
    t_t1rec = []; 
    w=b0;
    w1=0;
    w2=0;
    [t_t1rec,M_t1rec]=ode45('superL_transient',[start fin],Mall);
    Mini_t1rec = M_t1rec(end,:);
    Mini_t1rec(1) = 0;                     
    Mini_t1rec(2) = 0; 
    Mini_t1rec(3) = 0;
    Mini_t1rec(4) = 0; 
    Mini_t1rec(8) = 0;
    Mini_t1rec(9) = 0;    
    %%  Generate ahp-rahp                
    amscale = sla(ifo);
    [am_ahp, fm_ahp, am_rahp, fm_rahp] = gen_ahp_rahp(amscale_fx, fmscale_fx, hsdur_fx, hsn, adi_freq_fx, tp, beta, disprf);
    np_ahp = length(am_ahp);
    
    am_ahp=am_ahp*2*pi*b1; 
    fm_ahp=fm_ahp*2*pi;
    am_rahp=am_rahp*2*pi*b1;
    fm_rahp=fm_rahp*2*pi;
    Binh=Bexp*b1;
    a_rf_inh=r*Binh*2*pi;
    %%
    M=[];
    t=[];
    tf2=[];
    Mf2=[];
    tw=[];
    Mw=[];
    tw2=[];
    Mw2=[];
    Ma=0;
    Mb=0;
    Mc=0;
    ta=0;
    tb=0;
    tc=0;
   
    M1ini  = Mini_t1rec; 
    Mall = M1ini;
    Ma=M1ini;
    start = 0;
    fin   = 0;
    
    %% 180 flip
    if itog==1 
        for i = 1:np_afp 

            start = fin;
            fin   = start + diff_afp_time(i); 
            w=down_afp_rff(i)*2*pi+b0+AFP_cf*2*pi; 
            w1=down_afp_rfa(i)*2*pi;
            w2=0;

            [tf2,Mf2]=ode45('superL_transient',[start fin],Mall);
            Mall = Mf2(end,:);
        end
    else
        start = fin;
        fin   = start + down_afp_time(end);

        w=b0;
        w1=0;
        w2=0;

        [tf2,Mf2]=ode45('superL_transient',[start fin],Mall);
        Mall = Mf2(end,:);
    end

    %% crusher wait time 1
    if wait_time1~=0
        Mall(1)=0; 
        Mall(3)=0;
        
        start=fin;
        fin=start+wait_time1;
        w=b0;
        w1=0;
        w2=0;
        [tw,Mw]=ode45('superL_transient',[start fin],Mall); 
        Mall = Mw(end,:);
    end
    
    clc;
    fprintf('  Progress Report:\n');
    fprintf('  TOG Acquisition:     %d / 2\n', itog);
    fprintf('  FO:                  %d / %d\n', ifo, nfo);
    fprintf('  Training Set:         1 / %d\n', ntrain);
 
    %% AHP
    for i=0:(np_ahp-1)
        start_ahp=fin+i*(tp); 
        fin_ahp=fin+(i+1)*(tp);
        w=fm_ahp(i+1)+b0; 
        w1=am_ahp(i+1);
        w2=0;
        [t1,M1]=ode45('superL_transient',[start_ahp fin_ahp],Mall); 
        Mall = M1(end,:);
    end
    
    fin = fin_ahp;
    fin_temp = fin;
    Mall_temp = Mall;
    
    M2=[];
    t2=[];

    %%   
    if tsl ~=0
        start=fin;
        fin=start+ntrain*tsl; 
        w=freq_offset(ifo)*2*pi+b0; 
        w1=SL*b1; 
        w2=0;
        [t2,M2]=ode45('superL_transient',[start fin],Mall);
        Mall = M2(end,:);
    end
    %% RAHP
    for k=0:(np_ahp-1)
        start_rahp=fin+k*(tp);
        fin_rahp=fin+(k+1)*(tp);
        w=fm_rahp(k+1)+b0;
        w1=am_rahp(k+1);
        w2=0;
        [t3,M3]=ode45('superL_transient',[start_rahp fin_rahp],Mall);
        Mall_conv = M3(end,:);
    end

    %% save magnetization for orignial MPFSL
    MxA_conv(itog, ifo) = Mall_conv(1);
    MyA_conv(itog, ifo) = Mall_conv(3); 
    MzA_conv(itog, ifo) = Mall_conv(5);
    MzMT_conv(itog,ifo) = Mall_conv(7);
    MxMT_conv(itog,ifo) = Mall_conv(8);
    MyMT_conv(itog,ifo) = Mall_conv(9);
    %%
    if tsl ~=0                  
        start=fin_temp;
        fin=start+tsl; 
        w=freq_offset(ifo)*2*pi+b0; 
        w1=SL*b1; 
        w2=0;
        [t2,M2]=ode45('superL_transient',[start fin],Mall_temp);
        Mall = M2(end,:);
    end
    %% RAHP
    for k=0:(np_ahp-1)
        start_rahp=fin+k*(tp);
        fin_rahp=fin+(k+1)*(tp);
        w=fm_rahp(k+1)+b0;
        w1=am_rahp(k+1);
        w2=0;
        [t3,M3]=ode45('superL_transient',[start_rahp fin_rahp],Mall);
        Mall = M3(end,:);
    end
          
    MxA(itog, ifo,  1) = Mall(1);            
    MyA(itog, ifo,  1) = Mall(3); 
    MzA(itog, ifo,  1) = Mall(5);             
    MzMT(itog, ifo, 1) = Mall(7);
    MxMT(itog, ifo, 1) = Mall(8);
    MyMT(itog, ifo, 1) = Mall(9);
    %% train loop equal to number of tsl
    if ntrain>0                    
        for itrain = 1:ntrain-1 
            fin = fin_rahp; 
            if wait_time2~=0 
                Mall(1)=0;
                Mall(3)=0;
                start=fin;
                fin=start+wait_time2;
                w=b0;
                w1=0;
                w2=0;
                [tw2,Mw2]=ode45('superL_transient',[start fin],Mall);
                Mall = Mw2(end,:);
            end
            %% AHP
  
            clc;
            fprintf('  Progress Report:\n');
            fprintf('  TOG Acquisition:     %d / 2\n', itog);
            fprintf('  FO:                  %d / %d\n', ifo, nfo);
            fprintf('  Training Set:        %d / %d\n', itrain+1, ntrain);
            M4=[];
            t4=[];
            for i=0:(np_ahp-1)
                start_ahp=fin+i*(tp);
                fin_ahp=fin+(i+1)*(tp);
                w=fm_ahp(i+1)+b0;
                w1=am_ahp(i+1);
                w2=0;
                [t4,M4]=ode45('superL_transient',[start_ahp fin_ahp],Mall);
                Mall = M4(end,:);
            end
            
            fin = fin_ahp;            
            M5=[];
            t5=[];
            %%
            if tsl ~=0             
                start=fin;
                fin=start+tsl;
                w=freq_offset(ifo)*2*pi+b0;
                w1=SL*b1;
                w2=0;
                [t5,M5]=ode45('superL_transient',[start fin],Mall);
                Mall = M5(end,:);
            end
            %% RAHP
            M6=[];
            t6=[];
            for k=0:(np_ahp-1)
                start_rahp=fin+k*(tp);
                fin_rahp=fin+(k+1)*(tp);
                w=fm_rahp(k+1)+b0;
                w1=am_rahp(k+1);
                w2=0;
                [t6,M6]=ode45('superL_transient',[start_rahp fin_rahp],Mall);
                Mall = M6(end,:);

            end
            MxA(itog, ifo,  itrain+1) = Mall(1);             
            MyA(itog, ifo,  itrain+1) = Mall(3); 
            MzA(itog, ifo,  itrain+1) = Mall(5);             
            MzMT(itog, ifo, itrain+1) = Mall(7);
            MxMT(itog, ifo, itrain+1) = Mall(8);
            MyMT(itog, ifo, itrain+1) = Mall(9);   
        end
    end                             
end
end
toc

%% plot image
tsl_arr = linspace(1,10,10)*tsl;

rmpfsl_conv_gt = abs(-log((MzA_conv(2,2)-MzA_conv(1,2))./(MzA_conv(2,1)-MzA_conv(1,1)))./ntrain./tsl); 

rmpfsl_train_gt = zeros(1,ntrain);
for i = 1:ntrain
    rmpfsl_train_gt(i) = abs(-log((MzA(2,2,i)-MzA(1,2,i))./(MzA(2,1,i)-MzA(1,1,i)))./i./tsl); 
end

signal_level = MzA_conv(2,2);

db_array = [30, 40, 50, 60];
ndb = length(db_array);

snr_curve = zeros(ndb,ntrain);
snr_train = zeros(1, ntrain);

for m = 1:ndb
    snr_db = db_array(m);
    npoints = 500; 
    
    signal_power = signal_level.^2;
    noise_power = signal_power / (10^(snr_db/10));

    for i = 1:ntrain
        MzA_signal = repmat(MzA(:,:,i), [1,1,npoints]);
        noise1 = sqrt(noise_power)*randn(2,2,npoints);
        MzA_signal = MzA_signal+noise1;
        
        data = MzA_signal; 
        rmpfsl_train = abs(-log((data(2,2,:)-data(1,2,:))./(data(2,1,:)-data(1,1,:)))./tsl./i); 
        rmpfsl_train= squeeze(rmpfsl_train);
    
        train_s = rmpfsl_train_gt(i)*ones(npoints,1);
        snr_train(i) = snr(train_s, train_s-rmpfsl_train);
    end
    snr_curve(m,:) = snr_train/2/10;
end

figure
plot(snr_curve(1,:), '-o', 'LineWidth',2,'MarkerSize',4,'Color',[0, 0.4470, 0.7410]); 
hold on;
plot(snr_curve(2,:), '-s', 'LineWidth',2,'MarkerSize',4,'Color',[0.8500, 0.3250, 0.0980]); 
hold on;
plot(snr_curve(3,:),'-d', 'LineWidth',2,'MarkerSize',4,'Color',[0.9290, 0.6940, 0.1250]);
hold on;
plot(snr_curve(4,:), '-*', 'LineWidth',2,'MarkerSize',4,'Color',[0.4940, 0.1840, 0.5560]); 
legend('30 dB','40 dB','50 dB','60 dB');
title(['Tp = ' num2str(tsl*1e3) ' ms'],'FontSize', 15,'FontWeight', 'bold');
ylabel('log_{10}(RMP)','FontSize', 15,'FontWeight', 'bold');
xlabel('Number of spin-lock pulse (n)','FontSize', 15,'FontWeight', 'bold');
