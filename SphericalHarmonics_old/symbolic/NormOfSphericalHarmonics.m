%% 球面調和関数のノルム計算（内積計算）
% l,m の値が一致するもの同士のときは 1 になり、それ以外のときは 0 になる (正規直交性）
% 

syms theta phi

%% 球面調和関数の用意
Y00 = 1/2*sqrt(1/pi);

Y1m1 = 1/2*sqrt(3/2/pi)*sin(theta)*exp(-1i*phi);
Y10  = 1/2*sqrt(3/pi)*cos(theta);
Y11  = -1/2*sqrt(3/2/pi)*sin(theta)*exp(1i*phi);

Y2m2 = 1/4*sqrt(15/2/pi)*sin(theta)^2*exp(-2i*phi);
Y2m1 = 1/2*sqrt(15/2/pi)*sin(theta)*cos(theta)*exp(-1i*phi);
Y20 = 1/4*sqrt(5/pi)*(3*cos(theta)^2 - 1);
Y21 = -1/2*sqrt(15/2/pi)*sin(theta)*cos(theta)*exp(1i*phi);
Y22 = 1/4*sqrt(15/2/pi)*sin(theta)^2*exp(2i*phi);

%% 球面積分で内積計算
term1 = Y21;
term2 = Y21;

YY = int(int( ...
        term1 * conj(term2) * sin(theta), ...
    phi, 0, 2*pi), ... 
    theta, 0, pi);

double(YY)
