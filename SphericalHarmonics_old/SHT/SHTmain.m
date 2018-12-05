% 球面調和関数を用いた回転不変な特徴量抽出
% 
% スカラー実数球面調和関数を扱う。
% 

% CONSTANTS
L_max = 3;
RES = [55,55];

tic

% prepare INPUT FIELD
field1 = SHTinput();

field2 = ones(RES(1),RES(2));
field3 = SH(0,0)/5 + SH(1,1)/7 + SH(2,1) + SH(3,1)/2 + SH(7,1)/10;
field4 = SH(0,0)/5 + SH(1,-1)/7 + SH(2,-1) + SH(3,-1)/2 + SH(7,-1)/10;
field5 = rotSphere(SH(2,-1,RES),60,0);

field = field5;

% 球面調和変換(SHT)
coefficients = SHT(field,RES,L_max);


% norm calculation

% spectrum = sum(abs(coef),2);

% 次数l ごとにノルム計算
spectrum = zeros(1,L_max+1);
for il = 1 : L_max+1
    COEF = zeros(size(coef));
    COEF(il,:) = coef(il,:);
    gg1 = SHBT(COEF);
    % plotSH(gg1,1)     % 各次数ごとに再構成したときの球面波
    spectrum(il) = SHnorm(gg1);
end

toc

%% プロット
plotSH(field,1)
plotCoef(coef)
spectrum'
% real(coef), imag(coef)

% figure
% bar3(spectrum)
% zlabel('L2-Norm'),ylabel('次数l')
% ax = gca;
% ax.YTickLabel = num2str([0:orderMax]');

