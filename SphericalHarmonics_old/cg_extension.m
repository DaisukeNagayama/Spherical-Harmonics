function [C] = cg_extension(coef1,coef2)
% CG_EXTENSION  クレブシュ・ゴルダン係数を使用して既約展開
%   [C] = CG_EXTENSION(coef1,coef2)
%       coef1: 球面調和係数１
%       coef2: 球面調和係数２
%       C    : 既約展開した結果の球面調和係数
%   
%   coef1 の球面場と coef2 の球面場を掛けた結果を求める
%   
%   Yl1m1 * Yl2m2 = sum[L,M]{
%                       √((2l1+1)(2l2+1)/(4π(2L+1)))
%                       〈l1 0 l2 0 | L 0〉
%                       〈l1 m1 l2 m2 | L M〉
%                       YLM
%                   }

%   1/√4π 倍になってる？
% 
% coef1 = [0 1 0; 0 0 0];
% coef2 = [0 0 1 0 0; 0 0 1 0 0; 0 1 1 0 0];
% C = cg_extension(coef1,coef2)
% SHnorm(SHBT(C))

%%
L1 = size(coef1,1)-1;  % coef1 の最大次数 l1
L2 = size(coef2,1)-1;  % coef2 の最大次数 l2

L = abs(L1 - L2):L1 + L2;  % C の取りうる次数 L
C = zeros(L(end)+1,2*L(end)+1);

for il1 = 0:L1
    for il2 = 0:L2
        for im1 = -il1:il1
            for im2 = -il2:il2
                iL = il1 + il2;
                term1 = sqrt((2*il1+1)*(2*il2+1)/(4*pi*(2*iL+1)));
                term2 = Wigner3j(il1,il2,iL,0,0,0);
                for iM = -iL:iL
                    C(iL+1,L(end)+iM+1) = C(iL+1,L(end)+iM+1)...
                        + term1 * term2 ...
                            * Wigner3j(il1,il2,iL,im1,im2,iM) ...
                            * coef1(il1+1,L1+im1+1) ...
                            * coef2(il2+1,L2+im2+1);
                end
            end
        end
    end
end



