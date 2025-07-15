%-------------------------------------------------------------------------%
% Function:    superL_transient
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
%   Computes the time derivative (dM/dt) for the Bloch-McConnell equations
%   including a Super-Lorentzian lineshape for the macromolecular pool (MT).
%   This function is intended for use with ODE solvers for simulating transient
%   magnetization dynamics in the presence of chemical exchange and MT effects.
%-------------------------------------------------------------------------%

function f= superL_transient(t,M)

global T2c wc R1a R2a R1b R2b w wa wb w1 w2 kba kac kca M0a M0b M0c kab R1c R2c
wx=sqrt(w1^2 + w2^2);
wy=w+wc;

f(1)=-(R2a+kab)*M(1)+kba*M(2)-(wa-w)*M(3)-w2*M(5);                    %MxA
f(2)=kab*M(1)-(R2b+kba)*M(2)-(wb-w)*M(4)-w2*M(6);                     %MxB
f(3)=(wa-w)*M(1)-(R2a+kab)*M(3)+kba*M(4)+w1*M(5);                     %MyA
f(4)=(wb-w)*M(2)+kab*M(3)-(R2b+kba)*M(4)+w1*M(6);                     %MyB
f(5)=w2*M(1)-w1*M(3)-(R1a+kab+kac)*M(5)+kba*M(6)+kca*M(7)+R1a*M0a;    %MzA
f(6)=w2*M(2)-w1*M(4)+kab*M(5)-(R1b+kba)*M(6)+R1b*M0b;                 %MzB
f(7)=w2*M(8)-w1*M(9)+kac*M(5)-(R1c+kca)*M(7)+R1c*M0c;                 %MzMT

f(8)=-R2c*M(8)-(wa-w)*M(9)-w2*M(7);                                   %MxMT
f(9)=(wa-w)*M(8)-R2c*M(9) +w1*M(7);                                   %MyMT
f = f(:);
end
