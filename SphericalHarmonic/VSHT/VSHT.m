%% ベクトル球面調和変換 Vector Spherical Harmonics Transform
% 球面ベクトル場をベクトル球面調和関数に級数展開して2つの係数列を得る
% 
% 手順： SHT使ってあとは足し算
% 
% input: v_θ,v_φ ベクトル場のθ方向成分とφ方向成分  
%   次数 l = N まで求めるとする
%   
%   g_θ = <v_θ/sin(θ), Ylm>   % θ方向についての係数列
%   g_φ = <v_φ/sin(θ), Ylm>   % φ方向についての係数列
%   
%   各 l,m について(l = 1~N)
%     fB, fC <-- g_θ, g_φ
% output: fB, fC
% 
%% テスト実行時
%  [vt,vp] = VSHsample(); VSHT(vt,vp,10);
%  
%%  
% TODO: VSHBT は問題ないみたいなのでこっちが間違えている可能性高い
%       SHT に入れる 1/sin(θ) を飛ばしたらそれっぽくなった → なぜ？
%       位数がマイナスだとうまくいかないのでそこらへん？
%       次数が低いほど大きな値になってしまう
%
%% 
function [fB,fC] = VSHT(vt,vp,orderMax)
%% サンプリング点の数を読み込み
polarNum = size(vt,1);    % ルジャンドル多項式の次数（ガウス緯度を何点とるのか）
azimuthNum = size(vt,2);  % 経度方向に何点とるか(始点(0)と終点(2π)を両方カウント)

%% 通常の緯度経度
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

%% ガウス緯度μと経度λの計算
syms X;
LegPolynomial = legendreP(polarNum, X);     % polarNum次のルジャンドル多項式
mu_ = vpasolve(LegPolynomial == 0);          % ルジャンドル多項式の根がガウス緯度になる
mu_ = fliplr(double(mu_)');                    % シンボリックから倍精度数へ変換          
lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;    % 経度の刻み幅
[lambda,mu] = meshgrid(lambda_,mu_);

%% Gauss-Legendre 積分に使うガウス重みの計算
d_LegPolynomial = diff(LegPolynomial, X);
weight = double((((1-mu_.^2) .* (subs(d_LegPolynomial, mu_)).^2)/2).^-1);

%% 入力球面波を一般緯度からガウス緯度に変換する
vt = spline(theta(:,1), vt', acos(mu(:,1)))';
vp = spline(theta(:,1), vp', acos(mu(:,1)))';


%% SHT(Spherical Harmonics Transform algorithm)
orderMaxPlus = orderMax + 1;
gt = zeros(orderMaxPlus + 1, 2 * orderMaxPlus + 1);
gp = zeros(orderMaxPlus + 1, 2 * orderMaxPlus + 1);

for order = 0 : orderMaxPlus
    % Plm(μj) ルジャンドル陪関数
    Plms = legendre(order, mu(:,1)); % orderごとにまとめて生成
    for degree = -order : order
        % degree に対応する Plm を取り出す
        Plm = Plms(abs(degree) + 1,:);
        
        % 台形積分 横方向に積分
        GT = vt./(sin(acos(mu)) + eps) .* exp(-1i * degree * lambda);
        GT = trapz(lambda_, GT, 2)';
        GP = vp./(sin(acos(mu)) + eps) .* exp(-1i * degree * lambda);
        GP = trapz(lambda_, GP, 2)';

        % 符号と正規化
        sign = (-1)^((degree+abs(degree))/2);
        normSH = sqrt( (4*pi)/(2*order+1) * factorial(order+abs(degree)) / factorial(order-abs(degree)) );
        
        
        % ガウス・ルジャンドル積分の実行 式(B.44)
        gt(order+1, degree+orderMaxPlus+1) = sign/normSH * sum(weight .* GT .* Plm);
        gp(order+1, degree+orderMaxPlus+1) = sign/normSH * sum(weight .* GP .* Plm);
    end
end

%% fB, fC の計算
[DEGREE,ORDER] = meshgrid(-orderMaxPlus:orderMaxPlus, 0:orderMaxPlus);

% AA = ones(size(DEGREE));
% BB = ones(size(DEGREE));
% CC = ones(size(DEGREE));
AA = sqrt(ORDER.*(ORDER+1)).^-1;
BB = ORDER.*(ORDER-DEGREE+1) ./ (2*ORDER+1);
CC = (ORDER+1).*(ORDER+DEGREE) ./ (2*ORDER+1);

fB = AA.*( ...
    BB .* circshift(gt,[-1,0]) ...
    - CC .* circshift(gt,[1,0]) ...
    - 1i * ORDER .* gp);
fC = AA.*( ...
    -BB .* circshift(gp,[-1,0]) ...
    + CC .* circshift(gp,[1,0]) ...
    - 1i .* ORDER .* gt);

%% fB,fC を整形する
fB = fB(1:orderMaxPlus, 2:2*orderMaxPlus);
fC = fC(1:orderMaxPlus, 2:2*orderMaxPlus);

outOfRange = zeros(orderMax+1,2*orderMax+1) + [ ... % 範囲外を指定するための行列
    [flipud(tril(ones(orderMax))); zeros(1,orderMax)], ...
    zeros(orderMax+1,1), ...
    [triu(ones(orderMax)); zeros(1,orderMax)]   ];
outOfRange(1,:) = 1;

fB(1,:) = 0;
fB(outOfRange == 1) = 0;
fC(1,:) = 0;
fC(outOfRange == 1) = 0;

end % function SHT()





