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
%   Demo script for Simulation Study 1 in the MPF-PSL paper.
%
%   This script investigates the sensitivity of R_{mpfsl} and R_{mpfsl,pul}
%   to various tissue parameters, including R1a, R2a, T2b, Kca, and fb.
%
% Functions:
%   - Evaluates how R_{mpfsl} and R_{mpfsl,pul} respond to changes in key 
%     tissue parameters.
%
% History: 
%   - First built on 2025-07-13.
%
% ----------------------------------------------------------------------%
%%
clc; clear; close all;

% select tissues parameters we want to investigate
para = 'fb';

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

n_lines = 6;

% parameter select
switch para
    case 'fb'
        arr = mt*linspace(0,17,10)/100;
        nx = length(arr);
    case 'kca'
        mtr_arr = linspace(20,70,10);
        arr = mtr_arr;
        nx = length(arr);
    case 'R1a'
        arr = linspace(0.5,3,10);
        nx = length(arr);
    case 'R1c'
        arr = linspace(0.5,3,10);
        nx = length(arr);
    case 'R2a'
        arr=linspace(5,75,10);
        nx = length(arr);
    case 't2c'
        arr = linspace(6,12,10)*1e-6;
        nx = length(arr);
    otherwise  
        disp('Wrong parameter！');  
end

MxA = zeros(3,2, nfo, nx);
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

t1_recovery_time = 1.441;
TI = 0;
% Mag reset
Mini=[0 0 0 0 0 0 0 0 0];

% For the MPF-PSL method, we consider three different settings for the MPF-PSL parameters:
%   1) 10-10-50: spin-lock unit duration Tp = 10 ms, number of spin-lock units n = 10,
%                free precession time Tf = 50 ms
%   2) 20-5-50:  spin-lock unit duration Tp = 20 ms, number of spin-lock units n = 5,
%                free precession time Tf = 50 ms
%   3) 5-20-50:  spin-lock unit duration Tp = 5 ms, number of spin-lock units n = 20,
%                free precession time Tf = 50 ms
for scan = 1:3
    switch scan
        case 1
            Tf = 50;
            tsl = 10*1.e-3;
            ntrain = 10;
        case 2
            Tf = 50;
            tsl = 20*1.e-3;
            ntrain = 5;
        case 3
            Tf = 50;
            tsl = 5*1.e-3;
            ntrain = 20;
        otherwise  
            disp('Wrong scan scheme！');  
    end

for itog = 1:2 
for ifo = fost:foed 
    b1_fx = sla(ifo); 
    fo_fx = freq_offset(ifo); 
    b1=1;
    b0=0;
    wait_time2 = Tf *1.e-3;                
    for ix = 1:nx 
        %% initialize simulation parameters               
        M0b=pb; 
        M0c=0.069; 
        M0a=1-M0b-M0c; 
        kba=k_b; 
        kab=(M0b/M0a)*kba; 

        kca = 51; 
        kac = kca*(M0c/M0a); 
        SL=sla(ifo)*2*pi; 
        
        switch para
            case 'fb'                      
                M0b=pb; 
                M0c=arr(ix); 
                M0a=1-M0b-M0c; 
                kba=k_b; 
                kab=(M0b/M0a)*kba; 
                kca = 51; 
                kac = kca*(M0c/M0a);
            case 'kca' 
                kca = arr(ix);
                kac = kca*(M0c/M0a);
            case 'R1a'
                R1a = arr(ix);
                R1c = R1a;
            case 'R1c'
                R1c = arr(ix);
            case 'R2a'
                R2a = arr(ix);
            case 't2c'            
                t2c = arr(ix);                 
                T2c = t2c;
                R2c=mt*1/t2c;
            otherwise  
                disp('Wrong parameter！');  
        end
        
        % T1 recovery
        Mall = Mini;                      
        start = 0;
        fin   = t1_recovery_time - TI; 
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
                        
        %% tog pulse       
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

        %% crusher
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
        fprintf('  TOG Acquisition:      %d / 2\n', itog);
        fprintf('  FO:                  %d / %d\n', ifo, nfo);
        fprintf('  Parameter Index:      %d / %d\n', ix, nx);
        fprintf('  Training Set:         1 / %d\n', ntrain);
        fprintf('  Scan Mode:            %d / 3\n', scan);

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

        %% MPF-SL
        if tsl ~=0
            start=fin;
            fin=start+ntrain*tsl; 
            w=freq_offset(ifo)*2*pi+b0; 
            w1=SL*b1;
            w2=0;
            [t2,M2]=ode45('superL_transient',[start fin],Mall);
            Mall = M2(end,:);
        end
        %%  MPF-SL RAHP
        for k=0:(np_ahp-1)
            start_rahp=fin+k*(tp);
            fin_rahp=fin+(k+1)*(tp);
            w=fm_rahp(k+1)+b0;
            w1=am_rahp(k+1);
            w2=0;
            [t3,M3]=ode45('superL_transient',[start_rahp fin_rahp],Mall);
            Mall_conv = M3(end,:);
        end

        %% MPF-PSL
        if tsl ~=0 
            start=fin_temp;
            fin=start+tsl; 
            w=freq_offset(ifo)*2*pi+b0; 
            w1=SL*b1;
            w2=0;
            [t2,M2]=ode45('superL_transient',[start fin],Mall_temp);
            Mall = M2(end,:);
        end
        %% MPF-PSL RAHP
        for k=0:(np_ahp-1)
            start_rahp=fin+k*(tp);
            fin_rahp=fin+(k+1)*(tp);
            w=fm_rahp(k+1)+b0;
            w1=am_rahp(k+1);
            w2=0;
            [t3,M3]=ode45('superL_transient',[start_rahp fin_rahp],Mall);
            Mall = M3(end,:);
        end
        
        %% save magnetization for orignial MPFSL
        MxA_conv(scan, itog, ifo, ix) = Mall_conv(1);
        MyA_conv(scan, itog, ifo, ix) = Mall_conv(3); 
        MzA_conv(scan, itog, ifo, ix) = Mall_conv(5);
      
        MzMT_conv(scan, itog, ifo, ix) = Mall_conv(7);
        MxMT_conv(scan, itog, ifo, ix) = Mall_conv(8);
        MyMT_conv(scan, itog, ifo, ix) = Mall_conv(9);

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

                clc;
                fprintf('  Progress Report:\n');
                fprintf('  TOG Acquisition:      %d / 2\n', itog);
                fprintf('  FO:                  %d / %d\n', ifo, nfo);
                fprintf('  Parameter Index:      %d / %d\n', ix, nx);
                fprintf('  Training Set:         %d / %d\n', itrain+1, ntrain);
                fprintf('  Scan Mode:            %d / 3\n', scan);

                %% MPF-PSL AHP
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
                                       
            end
        end

        MxA(scan, itog, ifo, ix) = Mall(1);            
        MyA(scan, itog, ifo, ix) = Mall(3);              
        MzA(scan, itog, ifo, ix) = Mall(5);              
        MzMT(scan, itog, ifo, ix) = Mall(7);
        MxMT(scan, itog, ifo, ix) = Mall(8);
        MyMT(scan, itog, ifo, ix) = Mall(9);
    end
end          
end
end

toc

%%
sumtsl=ntrain*tsl;
for scan = 1:3
    for ix = 1:nx
        data = squeeze(MzA(scan,:,:,ix));
        rmpfsl_train(scan,ix) = abs(-log((data(2,2)-data(1,2))/(data(2,1)-data(1,1)))/sumtsl); 

        data_conv = squeeze(MzA_conv(scan,:,:,ix));
        rmpfsl_conv(scan,ix)= abs(-log((data_conv(2,2)-data_conv(1,2))/(data_conv(2,1)-data_conv(1,1)))/tsl/ntrain);
    end
end

figure;  
plot(arr, rmpfsl_conv(1,:), '-o', 'LineWidth',2,'MarkerSize',4,'Color',[0, 0.4470, 0.7410]); 
hold on;
plot(arr, rmpfsl_train(1,:), '-s', 'LineWidth',2,'MarkerSize',4,'Color',[0.8500, 0.3250, 0.0980]); 
hold on; 
plot(arr, rmpfsl_train(2,:),'-d', 'LineWidth',2,'MarkerSize',4,'Color',[0.9290, 0.6940, 0.1250]);
hold on;
plot(arr, rmpfsl_train(3,:), '-*', 'LineWidth',2,'MarkerSize',4,'Color',[0.4940, 0.1840, 0.5560]);
grid on;
legend('CW MPF-SL','MPF-PSL with 10-10-50','MPF-PSL with 20-5-50','MPF-PSL with 5-20-50');
xlabel(para)
ylabel('R_{mpfsl}[Hz]')