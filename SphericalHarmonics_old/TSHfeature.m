function spectrum = TSHfeature()
%% 3次元ボクセルデータ
%% 回転平行移動不変な特徴量の計算
%   3dラドン変換
% → 動径方向fft から振幅成分を抽出
% → テンソル球面調和変換
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
% voxel = zeros(16,16,16);
% voxel(1:16,7:10,7:10) = 1;
% voxel(7:10,1:16,7:10) = 1;
% voxel(7:10,7:10,1:16) = 1;

% Model 2
% voxel = zeros(16,16,16);
% voxel(1:16,7:10,7:10) = 1;
% voxel(7:10,1:16,7:10) = 1;
% voxel(7:10,7:10,1:16) = 1;
% voxel(7:10,1:4,1:16) = 1;

% Model 3
voxel = zeros(16,16,16);
voxel(1:16,7:10,7:10) = 1;
voxel(7:10,1:16,7:10) = 10;
voxel(7:10,7:10,1:16) = 1;
voxel(7:10,1:4,1:16) = 1;

% voxel = zeros(32,32,32);
% voxel(3:29,3:29,16) = 1;

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
% voxel = rot3d(voxel,0,30);
% 平行移動
% voxel = circshift(voxel,[2,3,-2]);

% モデルの表示
% plotVoxel(voxel)

%% 一次元投影の作成（３次元ラドン変換）
% Proj(Z,theta,phi)
Proj = projection3d(voxel, theta, phi);

%% 動径方向にフーリエ変換 (平行移動不変にする)
fProj = fft(Proj,[],1);
fProj = abs(fProj/size(fProj,1));           % パワー(振幅)スペクトルを抽出
fProj = fProj(1:ceil(size(fProj,1)/2),:,:); % 左右対称なので片側のみにする

%%% パワースペクトル比較（θ=π/2 部分の取り出し）
% figure
% bar3(permute(fProj( :,ceil(size(fProj,2)/2),: ),[1,3,2]))
% title('θ = π/2')
% xlabel('φ'),ylabel('周波数'),zlabel('パワースペクトル')
% ax = gca;
% ax.YTick = [];
% ax.XTick = [1,ceil(size(fProj,3)/2),size(fProj,3)]';
% ax.XTickLabel = char({'0', 'π', '2π'});

% figure
% for i = 1:49
% %     bar3(permute( Proj(i,:,:),[2,3,1] ));           % θ-φ表示
% %     bar3(permute( Proj(:,mod(i,9)+1,:),[1,3,2] ));  % Z-φ表示
%     bar3(permute( Proj(:,:,mod(i,16)+1),[1,2,3] )); % Z-θ表示
%     axis([0 inf 0 inf 0 200])
%     drawnow
%     pause(0.2)
% end



%% テンソル球面調和変換(TSHT)
n = nextpow2(size(fProj,1));
N = 2*n;
fn = TSHT(permute(fProj,[2,3,1]),N);

%% fn の mについての総和を計算（回転不変にする）
spectrum = permute(sum(abs(fn),2),[1,3,2]);

%%% 修正
% spectrum = spectrum(:,1:16);
% spectrum(2:2:6,:) = 100*spectrum(2:2:6,:);

%% プロット

figure
bar3(spectrum)
xlabel('基底'),ylabel('次数'),zlabel('振幅')
ax = gca;
ax.YTickLabel = num2str([abs(n):N]');










