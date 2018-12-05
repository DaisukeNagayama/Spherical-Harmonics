function plotfn(fn,base)
%% TSHT�̌��ʁA���Ԗڂ̊�ꂩ
%% �e���\�����ʒ��a�W�J�̌��ʂ�����
%   input : fn(l,m,base) �e��ꂲ�Ƃ̋��ʒ��a�֐��W�J�̌W����
%           base ���Ԗڂ̊��ɂ��Ă�
%
%%
N = (size(fn,2) - 1)/2;
n = N - size(fn,1) + 1;

for b = base
    
    C = fn(:,:,b);
    CMax = max(max(abs(C)));
    outOfRange = zeros(N+1,2*N+1) + [ ... % �͈͊O���w�肷�邽�߂̍s��
        [flipud(tril(ones(N))); zeros(1,N)], ...
        zeros(N+1,1), ...
        [triu(ones(N)); zeros(1,N)]   ];
    outOfRange = outOfRange(n+1:N+1,:);

    CR = real(C)/(CMax+eps);
    CI = imag(C)/(CMax+eps);

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
    colormap(myColorMap)
    cCR = colorbar;
    cCR.TickLabels = {'-','0','+'};
    axis on
    ax = gca;
    ax.XTick = [1:2:N,N+1,fliplr(2*N+1:-2:N+2)];
    ax.XTickLabel = num2str([-N:2:-1,0,fliplr(N:-2:1)]');
    ax.YTick = [1,fliplr(N+1:-2:2)];
    ax.YTickLabel = num2str([0,fliplr(N:-2:1)]');
    ax.TickLength = [0 0];
    ax.FontSize = 16;

    subplot(2,1,2)
    imgCI = imshow(CI,[-1,1],'InitialMagnification',2000);
    set(imgCI,'AlphaData',not(outOfRange));
    title('���ʒ��a�X�y�N�g���i�����j')
    xlabel('�x���i�ʐ��j degree'),ylabel('���� order')
    colormap(myColorMap)
    cCI = colorbar;
    cCI.TickLabels = {'-','0','+'};
    axis on
    ax = gca;
    ax.XTick = [1:2:N,N+1,fliplr(2*N+1:-2:N+2)];
    ax.XTickLabel = num2str([-N:2:-1,0,fliplr(N:-2:1)]');
    ax.YTick = [1,fliplr(N+1:-2:2)];
    ax.YTickLabel = num2str([0,fliplr(N:-2:1)]');
    ax.TickLength = [0 0];
    ax.FontSize = 16;

end

end