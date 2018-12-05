%%
% ベクトル球面調和変換(Vector Spherical Harmonics Transform)
% 
% 
% 通常の経度緯度の球面上の接線ベクトル場を扱う
% 
% 

%% 
orderMax = 5;

%%
% 元のベクトル場の用意
% [vt,vp] = VSHsample();

fB = zeros(orderMax+1, 2*orderMax+1);
fC = zeros(orderMax+1, 2*orderMax+1);
fB(4,orderMax+3) = 1;
fC(7,orderMax-4) = 1;
[vt,vp] = VSHBT(fB,fC);


%% スペクトルplot用のカラーマップを用意
    sg = (0:1/127:1)/2;
    sc = (1:-1/127:0)/2 + 0.5;
    myColorMap = sqrt([[sg; sg; sc]'; fliplr([sc; sg; sg])']);

%% スペクトルの範囲外を指定するための配列を用意
outOfRange = zeros(orderMax+1,2*orderMax+1) + [ ...
    [flipud(tril(ones(orderMax))); zeros(1,orderMax)], ...
    zeros(orderMax+1,1), ...
    [triu(ones(orderMax)); zeros(1,orderMax)]   ];
outOfRange(1,:) = 1;


%% スペクトル解析
    tic;
    [fB,fC] = VSHT(vt,vp,orderMax);
    toc

%% スペクトルから再構成
    tic;
    [wt,wp] = VSHBT(fB,fC);
    toc


%% 係数列を可視化
B = fB;             BMax = max(max(abs(B)));
BR = real(B)/BMax;  BI = imag(B)/BMax;
C = fC;             CMax = max(max(abs(C)));
CR = real(C)/CMax;  CI = imag(C)/CMax;

% figure(1)
%     % subplot(1,2,1)
%     subplot(2,2,1)
%     imgBR = imshow(BR,[-1,1],'InitialMagnification',2000);
%     set(imgBR,'AlphaData',not(outOfRange));
%     title('fB') %（実部）')
%     xlabel('度数（位数） degree'),ylabel('次数 order')
%     colormap(myColorMap)
%     colorbar
%     axis on
%     ax = gca;
%     ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
%     ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
%     ax.YTick = [1,fliplr(orderMax+1:-2:2)];
%     ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
%     ax.TickLength = [0 0];
%     ax.FontSize = 16;
% 
%     subplot(2,2,3)
%     imgBI = imshow(BI,[-1,1],'InitialMagnification',2000);
%     set(imgBI,'AlphaData',not(outOfRange));
%     title('fB（虚部）')
%     xlabel('度数（位数） degree'),ylabel('次数 order')
%     colormap(myColorMap)
%     colorbar
%     axis on
%     ax = gca;
%     ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
%     ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
%     ax.YTick = [1,fliplr(orderMax+1:-2:2)];
%     ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
%     ax.TickLength = [0 0];
%     ax.FontSize = 16;
% 
% 
%     % subplot(1,2,2)
%     subplot(2,2,2)
%     imgCR = imshow(CR,[-1,1],'InitialMagnification',2000);
%     set(imgCR,'AlphaData',not(outOfRange));
%     title('fC') %（実部）')
%     xlabel('度数（位数） degree'),ylabel('次数 order')
%     colormap(myColorMap)
%     colorbar
%     axis on
%     ax = gca;
%     ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
%     ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
%     ax.YTick = [1,fliplr(orderMax+1:-2:2)];
%     ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
%     ax.TickLength = [0 0];
%     ax.FontSize = 16;
% 
%     subplot(2,2,4)
%     imgCI = imshow(CI,[-1,1],'InitialMagnification',2000);
%     set(imgCI,'AlphaData',not(outOfRange));
%     title('fC（虚部）')
%     xlabel('度数（位数） degree'),ylabel('次数 order')
%     colormap(myColorMap)
%     colorbar
%     axis on
%     ax = gca;
%     ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
%     ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
%     ax.YTick = [1,fliplr(orderMax+1:-2:2)];
%     ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
%     ax.TickLength = [0 0];
%     ax.FontSize = 16;


%% 球面ベクトル場 の可視化
wt = real(wt); wp = real(wp);
vt = real(vt); vp = real(vp);
polarNum = size(wt,1);
azimuthNum = size(wt,2);

%%% 平面表示
% figure(2)
%     subplot(1,2,1)
%     quiver(1:azimuthNum,1:polarNum,wt,wp)
%     title('再構成ベクトル場'),xlabel('φ'),ylabel('θ')
%     axis tight
%     ax = gca;
%     ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
%     ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
%     ax.YTick = 1:(polarNum-1)/4:polarNum;
%     ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
%     ax.FontSize = 16;
% 
%     subplot(1,2,2)
%     quiver(1:azimuthNum,1:polarNum,vt,vp)
%     title('元のベクトル場'),xlabel('φ'),ylabel('θ')
%     axis tight
%     ax = gca;
%     ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
%     ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
%     ax.YTick = 1:(polarNum-1)/4:polarNum;
%     ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
%     ax.FontSize = 16;

figure
    subplot(1,2,1)
    quiver(1:azimuthNum,1:polarNum,vt,vp)
    title('ベクトル場１'),xlabel('φ'),ylabel('θ')
    axis tight
    ax = gca;
    ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
    ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
    ax.YTick = 1:(polarNum-1)/4:polarNum;
    ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
    ax.FontSize = 16;

    subplot(1,2,2)
    quiver(1:azimuthNum,1:polarNum,vt,vp)
    title('ベクトル場２'),xlabel('φ'),ylabel('θ')
    axis tight
    ax = gca;
    ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
    ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
    ax.YTick = 1:(polarNum-1)/4:polarNum;
    ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
    ax.FontSize = 16;
    
%% 球面表示
% rad = 101;      % ベクトル球の半径
% rad2 = 100;     % 地面の半径
% 
% % 風速
% tau_r = rad.*sin(theta);
% tau_x = tau_r.*cos(phi);
% tau_y = tau_r.*sin(phi);
% tau_z = rad.*cos(theta);
% tau_u = -wp.*sin(phi) - wt.*cos(theta).*cos(phi);
% tau_v = wp.*cos(phi) - wt.*cos(theta).*sin(phi);
% tau_w = wt.*sin(theta);
% 
% % 風速（原データ）
% rho_r = rad.*sin(theta);
% rho_x = rho_r.*cos(phi);
% rho_y = rho_r.*sin(phi);
% rho_z = rad.*cos(theta);
% rho_u = -vp.*sin(phi) - vt.*cos(theta).*cos(phi);
% rho_v = vp.*cos(phi) - vt.*cos(theta).*sin(phi);
% rho_w = vt.*sin(theta);
% 
% % 地面
% r2 = rad2.*sin(theta);
% x2 = r2.*cos(phi);
% y2 = r2.*sin(phi);
% z2 = rad2.*cos(theta);
% 
% % 描画
% figure
% subplot(1,2,1)
% quiver3(tau_x,tau_y,tau_z,tau_u,tau_v,tau_w, ...
%     'MaxHeadSize',1, ...
%     'LineWidth',2, ...
%     'AutoScaleFactor',2, ...
%     'Marker','.')
% axis equal off      % set axis equal and remove axis
% view(40,30)         % set viewpoint
% camzoom(1.5)        % zoom into scene %1.5
% hold on
% surf(x2,y2,z2, ...
%     'EdgeColor','none', ...
%     'FaceLighting','gouraud');
% colormap gray
% light
% hold off
% 
% subplot(1,2,2)
% quiver3(rho_x,rho_y,rho_z,rho_u,rho_v,rho_w, ...
%     'MaxHeadSize',1, ...
%     'LineWidth',2, ...
%     'AutoScaleFactor',2, ...
%     'Marker','.')
% axis equal off      % set axis equal and remove axis
% view(40,30)         % set viewpoint
% camzoom(1.5)        % zoom into scene %1.5
% hold on
% surf(x2,y2,z2, ...
%     'EdgeColor','none', ...
%     'FaceLighting','gouraud');
% colormap gray
% light
% hold off


