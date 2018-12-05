function plotfn(fn,base)
%% TSHTの結果、何番目の基底か
%% テンソル球面調和展開の結果を可視化
%   input : fn(l,m,base) 各基底ごとの球面調和関数展開の係数列
%           base 何番目の基底についてか
%
%%
N = (size(fn,2) - 1)/2;
n = N - size(fn,1) + 1;

for b = base
    
    C = fn(:,:,b);
    CMax = max(max(abs(C)));
    outOfRange = zeros(N+1,2*N+1) + [ ... % 範囲外を指定するための行列
        [flipud(tril(ones(N))); zeros(1,N)], ...
        zeros(N+1,1), ...
        [triu(ones(N)); zeros(1,N)]   ];
    outOfRange = outOfRange(n+1:N+1,:);

    CR = real(C)/(CMax+eps);
    CI = imag(C)/(CMax+eps);

    % カラーマップの用意
    sg = (0:1/127:1)/2;
    sc = (1:-1/127:0)/2 + 0.5;
    myColorMap = [[sg; sg; sc]'; fliplr([sc; sg; sg])'];

    % 出力
    figure
    subplot(2,1,1)
    imgCR = imshow(CR,[-1,1],'InitialMagnification',2000);
    set(imgCR,'AlphaData',not(outOfRange));
    whitebg('w')
    title('球面調和スペクトル（実部）')
    xlabel('度数（位数） degree'),ylabel('次数 order')
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
    title('球面調和スペクトル（虚部）')
    xlabel('度数（位数） degree'),ylabel('次数 order')
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