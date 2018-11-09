% Spherical Harmonic Transform
% 
% INPUTS:
% 
% field - Spherical field, [RES(1) x RES(2)]
% RES   - Vector of # of points to use [#Theta x #Phi points],[1x2] or [2x1]
% L_max - Max of SHT degree L, [1x1]
% 
% OUTPUTS:
% 
% coefficients  - Spherical harmonic coefficients, [L_max+1 x 2*L_max+1]
% 
% 
% 入力波と球面調和関数の内積をとってるだけ <f,Ylm>
% 
% 内積計算をするため縦横に積分する必要がある
% 緯度方向はガウス緯度でとってガウス・ルジャンドル積分公式
% 経度方向は台形積分。周期的なので離散フーリエ変換に書き換え可能
% 
% <対応関係>
%  ガウス緯度μ   <---> 極角 theta のコサイン
%  経度 λ　　    <---> 方位角 phi
% 
% MATLABのルジャンドル陪関数、完全正規化版と通常版との関係式
%  Plm('norm') = Plm() * N;
%  N = sqrt( (2*l+1) / (2*pi) * factorial(l-abs(m)) / factorial(l+abs(m)) );
%  
% 
% 参考Webページ
%  英語版 Wikipedia
%   https://en.wikipedia.org/wiki/Spherical_harmonics
%  B.1.3節「球面のスペクトル法」特に式(B.41)から式(B.44)
%   https://www.gfd-dennou.org/member/hiroki/file/BPmodel/Note_BPM.pdf
% 
%
% 実行テスト
%  SHT(SH(2,-1)+SH(1,1),2)
% 
% Nagayama Daisuke, 2016

%
function coefficients = SHT(field,RES,L_max)

if ((RES(1) ~= size(field,1)) || (RES(2) ~= size(field,2)))
    disp('error(SHT): The size of FIELD and RES must be equal.');
    return
end

THETA = linspace(0,  pi,RES(1));  % Azimuthal/Longitude/Circumferential
PHI   = linspace(0,2*pi,RES(2));  % Altitude /Latitude /Elevation
[PHI,THETA] = meshgrid(PHI,THETA);

% ガウス緯度μと経度λの計算
syms X;
LegPolynomial = legendreP(RES(1), X);     % polarNum次のルジャンドル多項式
mu_ = vpasolve(LegPolynomial == 0);          % ルジャンドル多項式の根がガウス緯度になる
mu_ = fliplr(double(mu_)');                     % シンボリックから倍精度数へ変換          
lambda_ = 0 : 2*pi / (RES(2) - 1) : 2*pi;    % 経度の刻み幅
[lambda,mu] = meshgrid(lambda_,mu_);

% Gauss-Legendre 積分に使うガウス重みの計算
d_LegPolynomial = diff(LegPolynomial, X);
weight = double((((1-mu_.^2) .* (subs(d_LegPolynomial, mu_)).^2)/2).^-1);

% 入力球面波を一般緯度からガウス緯度に変換する
field = spline(theta(:,1), field', acos(mu(:,1)))';


% SHT(Spherical Harmonic Transform algorithm)
% 係数格納用の変数、球面調和関数のピラミッド図みたいな直角二等辺三角形状に格納
coefficients = zeros(L_max + 1, 2 * L_max + 1);

for order = 0 : L_max    % 球面調和関数を低次から順に
    % Plm(μj) ルジャンドル陪関数
    Plms = legendre(order, mu(:,1)); % orderごとにまとめて生成
    
    for degree = -order : order % ループで回して係数計算する        
        
        % degree に対応する Plm を取り出す
        Plm = Plms(abs(degree) + 1,:);
        
        % 台形積分 横方向に積分
        G = field .* exp(-1i * degree * lambda);
        G = trapz(lambda_, G, 2)';
        
        % 符号と正規化
        sign = (-1)^((degree+abs(degree))/2);
        normSH = sqrt( (4*pi)/(2*order+1) * factorial(order+abs(degree)) / factorial(order-abs(degree)) );
        
        % ガウス・ルジャンドル積分 縦方向に積分 
        coefficients(order+1, degree+L_max+1) = ...
            sign/normSH * sum(weight .* Plm .* G);
        
    end
end

end % function SHT()



