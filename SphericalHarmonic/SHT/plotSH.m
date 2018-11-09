%% 球面スカラー場を可視化する
% 
% 2次元行列 yy の値を原点からの距離と見なしてプロットする
%
% 実行テスト
%   plotSH(SH(3,2));
%   plotSH(SH(3,-2),1);
%   plotSH(SH(9,3),2);
% 
%%
function plotSH(yy,scale)
%% scale(拡大縮小) 指定が無ければ1
if nargin < 2; scale = 1; end

%% Option
radius = 0;     % オフセット半径
alpha = 0.5;    % 座標軸サイズ [-α,α; -α,α; -α,α]

%% 虚部の計算をするかどうか
isComplex = not(isreal(yy));

%% データに合わせて球面グリッドの用意
polarNum = size(yy,1);
azimuthNum = size(yy,2);  
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

%% カラーマップの作成
% 負          中間　　　　　　　正
% Blue        Gray             Red
% [0,0,1] ～ [0.5,0.5,0.5] ～ [1,0,0]

sg = (0 :  1/127 : 1)/2;        % 0～0.5
sc = (1 : -1/127 : 0)/2 + 0.5;  % 1～0.5
myColorMap = [[sg; sg; sc]'; fliplr([sc; sg; sg])'];


%% データの読み込み 実部と虚部に分ける
rho = radius + scale * real(yy);

if(isComplex)
    tau = radius + scale * imag(yy);
end

%% rho と tau を直交座標に変換
rho_r = rho.*sin(theta);
rho_x = abs(rho_r).*cos(phi);
rho_y = abs(rho_r).*sin(phi);
rho_z = abs(rho).*cos(theta);

if(isComplex)
    tau_r = tau.*sin(theta);
    tau_x = abs(tau_r).*cos(phi);
    tau_y = abs(tau_r).*sin(phi);
    tau_z = abs(tau).*cos(theta);
end

%% プロット

if(isComplex)
    subplot(1,2,1)
    sR = surf(rho_x,rho_y,rho_z,rho_r);
    % title('Real','fontsize',22)
    colormap(myColorMap)
    light
    lighting gouraud
    axis equal off
    axis([-1 1 -1 1 -1 1]*alpha)
    view(40,30)
    camzoom(1.5)

    subplot(1,2,2)
    sI = surf(tau_x,tau_y,tau_z,tau_r);
    % title('Imag','fontsize',22)
    colormap(myColorMap)
    light
    lighting gouraud
    axis equal off
    axis([-1 1 -1 1 -1 1]*alpha)
    view(40,30)
    camzoom(1.5)
else
    sR = surf(rho_x,rho_y,rho_z,rho_r); 
        % 'EdgeColor','none');
    % title('','fontsize',22)
    colormap(myColorMap)
    light
    lighting gouraud
    axis equal off
    box off
    axis([-1 1 -1 1 -1 1]*alpha)
    view(40,30)
    camzoom(1.5)

end

end



