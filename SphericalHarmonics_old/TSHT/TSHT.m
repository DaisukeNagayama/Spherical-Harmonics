function fn = TSHT(f,N)
%% テンソル場, 級数展開の最大次数
%% テンソル場を球面調和関数で級数展開してスペクトルに変換する
%   input : f(θ,φ,base) テンソル場
%           N 級数展開の最大次数(l=0~N の球面調和関数を用いて展開する)
%   
%   output: fn(l,m,base) 各基底ごとの球面調和関数展開の係数列
%   
%   
%   
%% 実行テスト
% AAA = rotSphere(SH(1,0),0,0);
% BBB = rotSphere(SH(2,-1),0,0);
% test_data = cat(3,AAA,BBB); 
% test_fn = TSHT(test_data,4);
% plotSH(test_data(:,:,1)); plotSH(test_data(:,:,2));
% 
%% 定数の準備

n = nextpow2(size(f,3));  % テンソル場の階数

%% サンプリング点の数
elevationNum = size(f,1);    % ルジャンドル多項式の次数（ガウス緯度を何点とるのか）
azimuthNum = size(f,2);  % 経度方向に何点とるか(始点(0)と終点(2π)を両方カウント)

%% 通常の緯度経度
theta_ = 0 : pi/(elevationNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

%% ガウス緯度μと経度λの計算
syms X;
LegPolynomial = legendreP(elevationNum, X);     % elevationNum次のルジャンドル多項式
mu_ = vpasolve(LegPolynomial == 0);          % ルジャンドル多項式の根がガウス緯度になる
mu_ = fliplr(double(mu_)');                    % シンボリックから倍精度数へ変換          
lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;    % 経度の刻み幅
[lambda,mu] = meshgrid(lambda_,mu_);

%% Gauss-Legendre 積分に使うガウス重みの計算
d_LegPolynomial = diff(LegPolynomial, X);
weight = double((((1-mu_.^2) .* (subs(d_LegPolynomial, mu_)).^2)/2).^-1);


tic;
%% スケーリング定数の計算
% C(l,p,m,n)
%     access to C(l,p,m,n)
%     C(l -abs(n)+1, p +1, m +N+1, n +1);
C = makeC(n,N);

%% 緯度変換
%%% 入力球面波を一般緯度からガウス緯度に変換する %%%
for base = 1:size(f,3)
    f(:,:,base) = spline(theta(:,1), f(:,:,base)', acos(mu(:,1)))';
end

%% 補助関数gnの計算（内積計算）
% gn(l,m,base)
%     access to gn(l,m,base)
%     gn(l +1, m +N+abs(n)+1, base +1);
gn = zeros(N+abs(n)+1, 2*(N-abs(n))+1, 2^abs(n));
for l = 0:N+abs(n)
    % Plm(μj) ルジャンドル陪関数
    Plms = legendre(l, mu(:,1)); % lごとにまとめて生成
    
    for m = -l:l
        % degree に対応する Plm を取り出す
        Plm = Plms(abs(m) + 1,:);
        
        norm = 4*pi/(2*l+1) * factorial(l+abs(m)) / factorial(l-abs(m));
        
        meson = 1/norm * exp(-1i*m*lambda) ./(sin(acos(mu)).^abs(n) + eps);
        
        for base = 1:size(f,3)
            G = meson .* f(:,:,base);
            G = sum(G,2)' / sqrt(azimuthNum);
            
            gn(l +1, m +N+abs(n)+1,base) = sum(weight .* G .* Plm);
        end
    end
end

%% スペクトル計算 
% fn(l,m,base)
%     access to fn(l,m,base)
%     fn(l -abs(n)+1, m +N+1, base +1);
fn = zeros(N-abs(n)+1, 2*N+1, 2^abs(n));
for l = abs(n):N
    for m = -l:l
        fn_ = 0;
        for p = l-n:l+n
            fn_ = fn_ + conj( C(l -abs(n)+1, p +1, m +N+1, n +1) ) * gn(p +1,m +N+abs(n)+1,:);
%             fn_ = fn_ + gn(p +1,m +N+abs(n)+1,:);
        end
        fn(l -abs(n)+1, m +N+1, :) = fn_;
    end
end


end % function tenssor

%% function スケーリング定数Cの計算
function C = makeC(n,N)

l = abs(n):N;
p = 0:abs(n)+N;
m = -N:N;

% C(l,p,m,n)
C = zeros(length(l), length(p), length(m), abs(n)+1);

% n = 0 のとき
for jl = l
    C(jl -abs(n)+1, jl +1, :, 0 +1) = 1;
end

% n が 1以上のときを順に計算
for jn = 0:n-1
    for jl = l
    for jm = m
        for jp = p
            term1 = -jm * C(jl -abs(n)+1,jp +1,jm +N+1,jn +1);
            if jp-1 < p(1)
                term2 = 0;
            else
                term2 = (jp-jm)*(jp-2*jn-1)/(2*jp-1) * C(jl -abs(n)+1,jp-1 +1,jm +N+1,jn +1);
            end
            if jp+1 > p(length(p))
                term3 = 0;
            else
                term3 = -(jp+jm+1)*(jp+2+2*jn)/(2*jp+3) * C(jl -abs(n)+1,jp+1 +1,jm +N+1,jn +1);
            end
            
            C(jl -abs(n)+1,jp +1,jm +N+1,jn+1 +1) = ...
                1/(eps + sqrt(jl*(jl+1) - jn*(jn+1))) * ( term1 + term2 + term3 );
        end
    end
    end
end

end % function makeC



