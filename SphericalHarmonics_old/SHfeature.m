function SHfeature
%% 3次元ボクセルデータ
%% 回転平行移動不変な特徴量の計算
%   3dラドン変換
% → 動径方向fft から振幅成分を抽出
% → 各要素ごとに球面調和変換
% → 次数ごとに足し合わせてからL2-norm を計算
% 
% input : voxel(x,y,z) 3次元ボクセル画像
%     原点からの距離Z
%     y軸周りz軸→x軸方向の回転θ(0=<θ<π)
%     z軸周りx軸→y軸方向の回転φ(0=<φ<2π)
%     
% output: 特徴量()
% 
%
%% 原画像の用意
% Model 1
voxel = zeros(16,16,16);
voxel(1:16,7:10,7:10) = 1;
voxel(7:10,1:16,7:10) = 1;
voxel(7:10,7:10,1:16) = 1;

% Model 2
% voxel = zeros(16,16,16);
% voxel(1:16,7:10,7:10) = 1;
% voxel(7:10,1:16,7:10) = 1;
% voxel(7:10,7:10,1:16) = 1;
% voxel(7:10,1:4,1:16) = 1;

% Model 3
% voxel = zeros(16,16,16);
% voxel(1:16,7:10,7:10) = 1;
% voxel(7:10,1:16,7:10) = 10;
% voxel(7:10,7:10,1:16) = 1;
% voxel(7:10,1:4,1:16) = 1;

%% 投影角度の設定
theta = linspace(0,pi,8);
phi   = linspace(0,2*pi,8);

%% 回転してもはみ出さないように0埋め拡張と位置調整
oldSize = size(voxel);
newSize = ceil(sqrt(3)*size(voxel));
voxel(newSize(1),newSize(2),newSize(3)) = 0;
voxel = circshift(voxel, ceil(oldSize/2));

%% 回転, 平行移動
% 回転
% voxel = rot3d(voxel,pi/6,pi/6);
% 平行移動
% voxel = circshift(voxel,[2,3,-2]);

% モデルの表示
% plotVoxel(voxel)

%% 一次元投影の作成（３次元ラドン変換）
% Proj(Z,theta,phi)
disp('投影計算')
tic
Proj = projection3d(voxel, theta, phi);
toc

%% 動径方向にフーリエ変換 (平行移動不変にする)
% fProj(Z,theta,phi) -> fProj(theta,phi,Z)
disp('フーリエ変換')
tic

fProj = fft(Proj,[],1);
fProj = abs(fProj/size(fProj,1));           % パワー(振幅)スペクトルを抽出
fProj = fProj(1:ceil(size(fProj,1)/2),:,:); % 左右対称なので片側のみにする

fProj = permute(fProj,[2,3,1]);             % 次の処理のために並び替え

fProj = fProj(:,:,1:ceil(size(fProj,3)/8));

%%% パワースペクトル比較（中間層の取り出し）
figure
bar3(fProj( :,:,ceil(size(fProj,3)/2) ))
title('中間層')
xlabel('θ'),ylabel('φ'),zlabel('パワースペクトル')
ax = gca;
ax.XTick = [];
ax.YTick = [1,ceil(size(fProj,3)/2),size(fProj,3)]';
ax.YTickLabel = char({'0', 'π', '2π'});

toc

%% 球面調和変換(SHT)
disp('球面調和変換')
tic

orderMax = 6;
fn = zeros(orderMax + 1, 2 * orderMax + 1, size(fProj,3));
for iz = 1:size(fProj,3)
    fn(:,:,iz) = SHT(fProj(:,:,iz), orderMax);
end

%% fn の mについての総和を計算（回転不変にする）
spectrum = zeros(orderMax+1, size(fn,3));
for iz = 1:size(fn,3)
    for il = 1 : orderMax+1
        COEF = zeros(size(fn(:,:,1)));
        COEF(il,:) = fn(il,:,iz);
        gg1 = SHBT(COEF);
        % plotSH(gg1,1)     % 各次数ごとに再構成したときの球面波
        spectrum(il,iz) = SHnorm(gg1);
    end
end

toc

%% プロット

figure
bar3(spectrum)
xlabel('層'),ylabel('次数'),zlabel('ノルム')
ax = gca;
ax.YTick = [1:2:orderMax+1]';
ax.YTickLabel = num2str([0:2:orderMax]');

