function norm = SHnorm(yy)
%% 球面スカラー場のL2ノルムを求める
%       input : SH  球面スカラー場
% 
%       2乗してから球面で積分
%       テスト用: SHnorm(SH(2,0)), SHnorm(SH(2,1))
% 
%% サンプリング点の数を読み込み
polarNum = size(yy,1);    % ルジャンドル多項式の次数（ガウス緯度を何点とるのか）
azimuthNum = size(yy,2);  % 経度方向に何点とるか(始点(0)と終点(2π)を両方カウント)

%% 緯度経度
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

% ガウス緯度μと経度λの計算
syms X;
LegPolynomial = legendreP(polarNum, X);     % polarNum次のルジャンドル多項式
mu_ = vpasolve(LegPolynomial == 0);          % ルジャンドル多項式の根がガウス緯度になる
mu_ = fliplr(double(mu_)');                     % シンボリックから倍精度数へ変換          
lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;    % 経度の刻み幅
[lambda,mu] = meshgrid(lambda_,mu_);

%% Gauss-Legendre 積分に使うガウス重みの計算
d_LegPolynomial = diff(LegPolynomial, X);
weight = double((((1-mu_.^2) .* (subs(d_LegPolynomial, mu_)).^2)/2).^-1);

%% 入力球面波を一般緯度からガウス緯度に変換する
yy = spline(theta(:,1), yy', acos(mu(:,1)))';

%% ノルム計算
% 2乗してから球面で積分
G = abs(yy).^2;

% 台形求積（φ方向に積分）
G = trapz(G,2)' /2/pi;

% ガウス・ルジャンドル積分（θ方向に積分）
norm = sum(weight .* G);

%% 回転不変性のテスト
% 
% de = 30;
% da = 30;
% testMat = zeros(180/de + 1, 360/da + 1);
% for iEle = 0:de:180
%     for iAzi = 0:da:360
%         testMat(iEle/de+1,iAzi/da+1) = SHnorm(rotSphere(SH(1,1), iEle, iAzi));
%     end
% end
% 

end

