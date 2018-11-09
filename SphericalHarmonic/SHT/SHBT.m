% Spherical Harmonic Back Transform
% 
% INPUTS:
% 
% coefficients  - Spherical harmonic coefficients, [L_max+1 x 2*L_max+1]
% RES           - Vector of # of points to use [#Theta x #Phi points],[1x2] or [2x1]
% 
% OUTPUTS:
% 
% field_reconstructed - Reconstructed spherical field, [RES(1) x RES(2)]
% 
% テスト実行するとき
%  plotSH(SHBT(SHT(SH(2,-1),3)),1)
% 
% Nagayama Daisuke, 2016


function field_reconstructed = SHBT(coefficients,RES)

%%% 最大次数の読み込み %%%
L_max = size(coefficients,1) - 1;

%%% SHBT(Spherical Harmonics Back Transform algorithm) %%%
f = zeros(RES(1), RES(2));
for iL = 0 : L_max
    for iM = -iL : iL    
        if coefficients(iL+1,iM+L_max+1) == 0
            continue
        end
        
        f = f + coefficients(iL+1,iM+L_max+1) * SH(iL,iM,RES);
    end
end
field_reconstructed = f;

end % function SHBT(s)





