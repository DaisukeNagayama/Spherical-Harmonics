% 球面調和関数の積（合成）
% 
% ２つの球面調和関数の積を公式をもとに計算する
% 球面調和関数の線形和で表せるらしい
% 係数は、Clebsch-Gordan 係数や Wigner 3-j 記号 を用いて計算する
% 
%   近い形状にはなるけれど一致しない
%   何か解釈間違い？
% 
%%
% [Clebsch Gordan 係数]
% |l1 - l2| <= L <= l1 + l2
% m1 + m2 + M = 0
% 
% SH(l1,m1) * SH(l2,m2)
% = sum[L,M]
%   sqrt( (2*l1+1)*(2*l2+1) / 4*pi*(2*L+1) )
%   * ClebschGordan(l1,l2,L,0,0,0) * ClebschGordan(l1,l2,L,m1,m2,M)
%   * SH(L,M)
% 
% [Wigner 3-j symbol]
% 
% SH(l1,m1) * SH(l2,m2)
% = sum[l,m]
%   sqrt( (2*l1+1)*(2*l2+1)*(2*L+1) / 4*pi )
%   * Wigner3j(l1,l2,L,0,0,0)
%   * Wigner3j(l1,l2,L,m1,m2,-M)
%   * conj(SH(L,M))
% 
%%
l1 = 1; m1 = 1;
l2 = 2; m2 = -1;

L = abs(l1 - l2) : l1 + l2;
M = m1 + m2;

%% Clebsch Gordan
% %{
right = zeros(size(SH(0,0)));
for il = L
    right = right ...
        + sqrt( (2*l1+1)*(2*l2+1) / 4*pi*(2*il+1) ) ...
        * ClebschGordan(l1,l2,il,0,0,0) * ClebschGordan(l1,l2,il,m1,m2,M) ...
        * SH(il,M);
end

left = SH(l1,m1).*SH(l2,m2);
%}

%% Wigner 3-j
%{
right = zeros(size(SH(0,0)));
for il = L
    for im = -il:il
        right = right ...
            + sqrt( (2*l1+1)*(2*l2+1)*(2*il+1) / 4*pi ) ...
            * Wigner3j(l1,l2,il,0,0,0) * Wigner3j(l1,l2,il,m1,m2,-im) ...
            * conj(SH(il,im));
    end
end

left = SH(l1,m1).*SH(l2,m2);
%}

%% Result

figure, plotSH(right)
figure, plotSH(left)
