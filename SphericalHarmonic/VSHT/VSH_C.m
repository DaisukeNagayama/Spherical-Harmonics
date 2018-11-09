%% ベクトル球面調和関数 Clm を 計算する
% VSH_B.m と VSH_C.m はセット
% 
% "Computational Harmonic Analysis for Tensor Fields on The Two-Sphere"
% 
%%
function Clm = VSH_C(l,m)
%% 球面グリッドの交点の数
polarNum = 41;      % 北極から南極 (0～π)
azimuthNum = 41;    % 一周 (0～2π)

%% 球面グリッドの作成
theta_ = 0 : pi/(polarNum - 1) : pi;     % polar angle
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi; % azimuth angle
[phi,theta] = meshgrid(phi_,theta_);      % define the grid

%% スカラー球面調和関数
Ylm = zeros(polarNum, azimuthNum,3);
for k = 1:3
    kl = l + k - 2;     % kl = l-1, l, l+1
    
    % ルジャンドル陪関数の用意
    Plm = legendre(kl, cos(theta(:,1))); % 各行に Pl0 から Pll までが並んだ行列
    Plm = Plm(abs(m) + 1,:)';           % Plm の行を取り出して縦にする
    Plm = repmat(Plm,1,size(phi,2));    % 球面グリッドのサイズに合わせる

    % 符号と正規化
    sign = (-1)^((m+abs(m))/2);
    normSH = sqrt( 4*pi / (2*kl+1) * factorial(kl+abs(m)) / factorial(kl-abs(m)) );

    % 複素 球面調和関数
    Ylm(:,:,k) = sign/normSH * Plm .* exp(1i * m * phi);
end

%% ベクトル球面調和関数
term0 = (sqrt(l*(l+1)) * sin(theta) + eps).^-1;
term1 = (l+1)*(l+m) / (2*l+1);
term2 = 1i*m;
term3 = -l*(l-m+1) / (2*l+1);

Clm = zeros(polarNum, azimuthNum,2);
Clm(:,:,1) = term0 .* (term2 * Ylm(:,:,2));                         % θ成分
Clm(:,:,2) = term0 .* (term1 * Ylm(:,:,1) + term3 * Ylm(:,:,3));    % φ成分


end
