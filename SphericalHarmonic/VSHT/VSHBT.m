%% ベクトル球面調和 逆変換 Vector Spherical Harmonics Back Transform
% 2つの係数列に従ってベクトル球面調和関数を足し合わせて球面ベクトル場を得る
% 
% ベクトル場 {vt,vp} を求める
% 原データ {wt,wp}
% 
% 
% 1 =< l, |m| =< l =< orderMax について求める
%  Blm
%  θ成分 → Ylm のθ方向に微分したものを√(l(l+1))で割る
%  φ成分 → Ylm に im/sinθ を掛けたものを√(l(l+1))で割る
%  Clm
%   中心向き単位ベクトル -er とBlm の外積
%   θ成分 → Blm のφ成分
%   φ成分 → Blm のθ成分のマイナス
% 
%
% 実行時
%  [vt,vp] = VSHsample(); [fB,fC] = VSHT(vt,vp,10); VSHBT(fB,fC);
%
%  % fB(l=2,m=1) を生成してsampleと比較実験 → 結果: 一致、問題なし
%  fC = zeros(11,21);
%  fB = zeros(11,21); fB(3,12) = 1;
%  VSHBT(fB,fC);
% 
%%
function [wt,wp] = VSHBT(fB,fC)
%% 最大次数の読み込み
orderMax = size(fB,1) - 1;

%% 刻み数を指定
polarNum = 41;      % theta  ルジャンドル多項式の次数（ガウス緯度を何点とるのか）
azimuthNum = 41;    % phi  経度方向に何点とるか(始点(0)と終点(2π)を両方カウント)

% 通常の緯度経度
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

% %%% ガウス緯度μと経度λの計算 %%%
% syms X;
% LegPolynomial = legendreP(polarNum, X); % polarNum次のルジャンドル多項式
% mu_ = vpasolve(LegPolynomial == 0);      % ルジャンドル多項式の根がガウス緯度になる
% mu_ = fliplr(double(mu_)');                % シンボリックから倍精度数へ変換
% lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;
% [lambda,mu] = meshgrid(lambda_,mu_);

%% VSHBT の本体
wt_prov = zeros(polarNum,azimuthNum);
wp_prov = zeros(polarNum,azimuthNum);
for order = 1 : orderMax
for degree = -order : order
    
    N = (-1)^((degree+abs(degree))/2) * sqrt( (2*order+1) / (4*pi) * ...
            factorial(order-abs(degree)) / factorial(order+abs(degree)) );
    
    % Ylm
    Plm = legendre(order,cos(theta(:,1)));
    Plm = Plm(abs(degree) + 1,:)';
    Plm = repmat(Plm,1,size(phi,2));
    Ylm = N * Plm .* exp(1i * degree * phi);
    
    % Blm
    AA = 1/sqrt(order * (order+1));
    dtYlm = diff(Ylm,1,1) * 2*pi;
    dtYlmW = [dtYlm(1,:); dtYlm];
    Blm_t = AA * dtYlmW;
    Blm_p = AA * 1i * degree * Ylm./(sin(theta)+eps);
    
    % Clm
    Clm_t = Blm_p;
    Clm_p = -Blm_t;
    
    % 足していく
    wt_prov = wt_prov + ( ...
        fB(order + 1,degree + orderMax + 1) * Blm_t ...
        + fC(order + 1,degree + orderMax + 1) * Clm_t );
    wp_prov = wp_prov + ( ...
        fB(order + 1,degree + orderMax + 1) * Blm_p ...
        + fC(order + 1,degree + orderMax + 1) * Clm_p );
end
end
wt = wt_prov;
wp = wp_prov;


% % 再構成画像をガウス緯度から一般緯度に変換する
% wt = spline(pi - acos(mu(:,1)), wt', theta(:,1))';
% wp = spline(pi - acos(mu(:,1)), wp', theta(:,1))';


end % function vectorSHBT(s)


