%-------------------------------------------------------------------------%
% Function:    gen_ahp_rahp
% Authors:     Qianxue Shan, Ziqiang Yu, Weitian Chen
% Email:       qxshan@link.cuhk.edu.hk
% Version:     1.0
% Date:        2025-07-13
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
%   Generates Adiabatic Half-Passage (AHP) and reverse AHP (rAHP) RF pulses.
%   Supports on-resonance and off-resonance T1rho applications.
%
%   Inputs:
%     amscale     (Hz):   Maximum B1 amplitude of AHP and rAHP (default: 500)
%     fmscale     (Hz):   Maximum frequency sweep of AHP and rAHP (default: 500)
%     hsdur      (s):     Duration of AHP or rAHP pulse (default: 25e-3)
%     hsn              :  HSn, n=1,8... for on-resonance T1rho; HS1 for off-resonance T1rho (default: 1)
%     freq_offset (Hz):   Frequency offset for off-resonance (default: 0)
%     tp         (s):     Sampling rate (default: 6.4e-6)
%     beta      (unitless): Shape parameter for AHP and rAHP (default: 4.0)
%     disprf    (0 or 1): Flag to display RF pulse (default: 0)
%
%   Outputs:
%     am_ahp (Hz):   Amplitude modulation for AHP
%     fm_ahp (Hz):   Frequency modulation for AHP
%     am_rahp (Hz):  Amplitude modulation for rAHP
%     fm_rahp (Hz):  Frequency modulation for rAHP
%
%   Note:
%     Further calculation of FM waveform is needed for off-resonance T1rho 
%     with HSn (n>1).
%
%   Do not modify, remove, or redistribute this header or any part of the code 
%   without explicit written permission from the authors.
%
%-------------------------------------------------------------------------%
function [am_ahp fm_ahp am_rahp fm_rahp] = gen_ahp_rahp(amscale, fmscale, hsdur, hsn, freq_offset, tp, beta, disprf)

if nargin <8 
    disprf = 0; 
    if nargin <7 
        beta = 4.0; 
        if nargin < 6 
           tp = 6.4*1.e-6; % sampling rate in second 
           if nargin < 5
             freq_offset = 0; 
             if nargin < 4
                hsn = 1; 
                if nargin < 3
                   hsdur = 25*1.e-3; 
                   if nargin < 2
                        fmscale = 500; 
                        if nargin < 1
                            amscale = 500; 
                        end
                   end
                end
             end
           end
        end
    end
end

% freq_offset = abs(freq_offset); %% JBY
                           
if (freq_offset ~= 0 & hsn > 1)
    hsn = 1; 
    warning('Change HSn to 1 for off-resonance T1rho'); 
end

% In PPE, this is written outside the cal_ahp_rahp function
% if abs(freq_offset) > 0.0
%     fmscale = abs(freq_offset + sign(freq_offset)*fmscale); 
% end

theta = atan2(amscale, freq_offset ); 

N = round(hsdur/tp); 
 
amtd = zeros(N,1); 
fmtd = amtd; 
amtu = amtd; 
fmtu = amtd; 
tt   = amtd; 

for nn=1:N
    amtd(nn) =0.0; 
    fmtd(nn) = 0.0; 
    amtu(nn) = 0.0;
    fmtu(nn) = 0.0; 
end

npw = N; 


if (freq_offset ~= 0.0)
    if (freq_offset > 0.0)     
        for nn = 0:npw-1          
            t  = (N-1-nn)/(N-1); 
            amtd( nn+1) = 1./cosh(beta*t) - 1/cosh(beta*1.0) ; 
            fmtd( nn+1) = tanh(beta*t); 
        end
    else
        for nn=0:npw-1
            t = (nn - (N-1) ) / (N-1);                         
            amtd( nn+1) = 1./cosh(beta*t) - 1/cosh(-beta*1.0); 
            fmtd( nn+1) = tanh(beta*t);                        
        end
    end      
    %%JBY reactivate this line
    amtd = amscale*amtd/max(abs(amtd(:))); % in unit Hz % move this line
    %after chopping the length of the array to npw
    fmtd = fmscale*fmtd/max(abs(fmtd(:)));  % in unit Hz
    %%JBY add frequency offset
    fmtd = fmtd + freq_offset;
    
    tmpv2 = (1.0 +  freq_offset/fmscale ) / (1.0 - freq_offset/fmscale ); 

    %%JBY comment these
%     if (tmpv2 <= 0.0)
%         tmpv = 0.0; 
%     else
%         tmpv = (0.5/beta)*log(tmpv2);
%     end
%     if (freq_offset > 0.0)
%         npw = round(((N-1)-(N-1)*tmpv)) + 2; 
%     else
%         npw = round((N-1)*tmpv + (N-1) ) + 2;                  
%     end
% 
%     amtd = amtd(1:npw); 
%     amtd = amscale*amtd/max(abs(amtd(:))); % in unit Hz
%     fmtd = fmtd(1:npw);    
    
    
    
    amtu = zeros(size(amtd)); 
    fmtu = amtu; 

    for nn=0:npw-1
        amtu(npw-nn) = amtd(nn+1); % in unit Hz
        fmtu(npw-nn) = fmtd(nn+1); % in unit Hz
    end

% 
%     if 1
%         fmtupc = 2*freq_offset-fmtu;
%     else
%         fmtupc = -fmtu; 
%     end  

else   

    for nn=0:N-1
        %t = (2*nn-(N-1)) / (N-1) ; 
        t = (nn-(N-1)) / (N-1) ; 
        tt(nn+1) = t; 
        %amtd(nn+1) = 1./(cosh(beta*t)) - 1./cosh(-beta*1.0); 
        %amtd(nn+1) = 1./(cosh(beta*t)); 
        
        % set the AM form to zero at the beginning
        if 0 
            amtd(nn+1) = sech(beta*t^hsn); 
        else
            amtd(nn+1) = sech(beta*t^hsn) - sech(beta*(-1.0)^hsn); 
        end
        tmpw(nn+1) = (sech(beta*t^hsn)).^2; 
        fmtd2(nn+1) = tanh(beta*t); 
    end
    fmtd = cumsum(tmpw); 
    fmtd = fmtd(:); 
    fmtd = fmtd/max(fmtd) ; 
    fmtd = fmtd - 1.0;
    amtd = amscale*amtd/max(abs(amtd(:))); 
    fmtd = fmscale*fmtd/max(abs(fmtd(:)));               
    for nn=0:N-1 
        amtu(N-nn) = amtd(nn+1); 
        fmtu(N-nn)    = fmtd(nn+1);       
    end
end

am_ahp  = amtd; 
am_rahp = amtu; 
fm_ahp  = fmtd; 
fm_rahp = fmtu; 

if disprf
    np = length(am_ahp); 
    np1 = length(am_rahp); 
    np2 = length(fm_ahp); 
    np3 = length(fm_rahp); 
    if ( (np~=np1) || (np ~= np2) || (np~=np3) )
        error('The duration of waveform does not match'); 
    end
    
    figure; 
    subplot(221); plot(am_ahp);  title(sprintf('AM of AHP, duration %d ms',  round(hsdur*1000))); 
    xlabel('sample point'); ylabel('Amplitude (Hz)'); axis([0 np+1 min(am_ahp)-5 max(am_ahp) + 5]); 
    
    subplot(222); plot(am_rahp); title('AM of reverse AHP');  
    xlabel('sample point'); ylabel('Amplitude (Hz)'); axis([0 np+1 min(am_rahp)-5 max(am_rahp) + 5]); 
    
    subplot(223); plot(fm_ahp);  title('FM of AHP');   
    xlabel('sample point'); ylabel('Frequency (Hz)'); axis([0 np+1 min(fm_ahp)-5 max(fm_ahp) + 5]); 
    subplot(224); plot(fm_rahp); title('FM of reverse AHP'); 
    xlabel('sample point'); ylabel('Frequency (Hz)'); axis([0 np+1 min(fm_rahp)-5 max(fm_rahp) + 5]);  
end
    
    
    
    
  
