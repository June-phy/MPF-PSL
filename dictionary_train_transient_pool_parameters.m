%-------------------------------------------------------------------------%
% Script:      bm_simulation_params.m
% Authors:     Qianxue Shan, Ziqiang Yu, Weitian Chen
% Email:       qxshan@link.cuhk.edu.hk
% Version:     1.0
% Date:        2025-07-13
%
% Copyright (c) 2025 Qianxue Shan, Ziqiang Yu, Weitian Chen. All rights reserved.
%
% License and Usage Notice:
%   This script is provided strictly for academic and research purposes only.
%   Any commercial use, including but not limited to sale, redistribution,
%   or integration into proprietary software, is strictly prohibited without
%   explicit written permission from the authors.
%
%   Modification of this script, its header comments, or removal of this notice,
%   in whole or in part, is EXPRESSLY FORBIDDEN without prior written consent
%   from the authors.
%
%   By using, copying, or referencing this script, you agree to abide by these terms.
%   For any inquiries or requests, please contact the authors at the email above.
%
% Description:
%   Defines tissue parameters and sequence parameters used in Bloch-McConnell (BM)
%   simulations for liver tissue.
%
%   DO NOT modify, remove, or redistribute this header or any part of the script
%   without explicit written permission from the authors.
%
%-------------------------------------------------------------------------%

global T2c R1a R2a R1b R2b R1c R1d R2d R1e R2e w wa wb wc wd we w1 w2 kab kba kad kda kae kea kac kca M0a M0b M0c M0d M0e Rfc_temp R2c
%% For toggling simulation
itoggle = 1;
repeat = 2;
%%  liver
T1=812*1.e-3; 
T2=42*1.e-3;
t2c=7.7*1.e-6; 

Td = 50;
Tp = 10;

mtr_arr = [51];
ntrain_arr= 10; 
%% crusher wait time
wait_time1=2*1.e-3;
wait_time2=Td*1.e-3;
%% MT and CE switch
cb=0; 
mt=1;
%% select different tissue parameter
tissue=1;
t2c=mt*[t2c];
mtr_arr=mt*mtr_arr;
pc = 0.069; 
npc = length(pc);
%%
R1=1./T1;
R2=1./T2;
T2c=t2c;
R2c=mt*1/t2c;
R1a=R1;
R1b=cb*R1;
R1c=mt*R1; 
R2a=R2; 
R2b=cb*R2;
%% chemical shift 
wa=0;
wb=cb*1*128*2*pi;%Hz*2pi
wc=0*128*2*pi;%Hz*2pi
%% pool population
pb=cb*linspace(0.01,0.01,1);
pb=cb*pb;
npb=length(pb);
pd=cd*0.01;
%% chemical exchange rate 
Cestr=1500;
k_b=[cb*Cestr];
kba=k_b;
kca_arr=mtr_arr;

%% B1 and B0 inhomogeneity
b0_arr=linspace(-100,100,10);

nb0=length(b0_arr);
b0st=1;
b0ed=nb0;

b1=1;
%% tsl setting
tsl = Tp;
tsl=tsl*1.e-3;

sla=[400 100];
nsl = length(sla);
slst=1;
sled=nsl;
%% frequency offset
freq_offset = [4000 1000];
nfo = length(freq_offset);
fost=1;
foed=nfo;
%% r gyromagnetic ragio
r=42.58e6;
Bexp=13.5*1.e-6;
flip_angl_180=pi;
a_rf_exp=r*Bexp*2*pi;
pw_rf_180=abs(flip_angl_180/a_rf_exp);

flip_angl_90=pi/2;
pw_rf_90=abs(flip_angl_90/a_rf_exp);
%% adiabatic pulse parameters
if adi_set == 1
    fmscale = 150;
    hsdur =[35];
    hsdur=hsdur*1.e-3;
    nhs=length(hsdur);
    hsst=1;
    hsed=1;
    hsn = 1;
    tp = 3.2*1.e-5;
    beta = 1;
    disprf = 0;
end
