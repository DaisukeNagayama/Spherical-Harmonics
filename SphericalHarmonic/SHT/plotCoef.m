%%%%%%%%%%%%%%%%%%%
% �W����̃v���b�g
% 
% 
%  
%  [0 0 0 0.9 0 0 0;
%   0 0 0.3 0.2 0.4 0 0;
%   0 0.2 -0.1 0.0 0.2 0.0 0;
%   0.1 0.1 -0.1 -0.1 0.1 0.05 -0.05]
% 
%
function plotCoef(coef)

orderMax = size(coef,1) - 1;
CMax = 1; % CMax = max(max(abs(coef)));
outOfRange = zeros(orderMax+1,2*orderMax+1) + [ ... % �͈͊O���w�肷�邽�߂̍s��
    [flipud(tril(ones(orderMax))); zeros(1,orderMax)], ...
    zeros(orderMax+1,1), ...
    [triu(ones(orderMax)); zeros(1,orderMax)]   ];

CR = real(coef)/CMax;
CI = imag(coef)/CMax;

% �J���[�}�b�v�̗p��
sg = (0:1/127:1)/2;
sc = (1:-1/127:0)/2 + 0.5;
myColorMap = [[sg; sg; sc]'; fliplr([sc; sg; sg])'];

% �o��
figure
subplot(2,1,1)
imgCR = imshow(CR,[-1,1],'InitialMagnification',2000);
set(imgCR,'AlphaData',not(outOfRange));
whitebg('w')
title('���ʒ��a�X�y�N�g���i�����j')
xlabel('�x���i�ʐ��j degree'),ylabel('���� order')
cCR = colorbar;
cCR.TickLabels = {'-','0','+'};
axis on
ax = gca;
ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
ax.YTick = [1,fliplr(orderMax+1:-2:2)];
ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
ax.TickLength = [0 0];
ax.FontSize = 16;
colormap(ax,myColorMap)

subplot(2,1,2)
imgCI = imshow(CI,[-1,1],'InitialMagnification',2000);
set(imgCI,'AlphaData',not(outOfRange));
title('���ʒ��a�X�y�N�g���i�����j')
xlabel('�x���i�ʐ��j degree'),ylabel('���� order')
cCI = colorbar;
cCI.TickLabels = {'-','0','+'};
axis on
ax = gca;
ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
ax.YTick = [1,fliplr(orderMax+1:-2:2)];
ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
ax.TickLength = [0 0];
ax.FontSize = 16;
colormap(ax,myColorMap)

end
