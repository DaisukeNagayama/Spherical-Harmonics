%% シンボリックで球面調和関数
% 

syms x
assume(x >= -1 & x <= 1)
syms theta phi
assume(theta >= 0 & theta <= pi)
assume(phi >= 0 & phi <= 2*pi)
% syms l m
% assume(in(l,'integer') & l >= 0)
% assume(in(m,'integer') )% & m >= -l & m <= l)

%% l,m の指定
l = 2;
m = 2;


%% 球面調和関数の計算
% ルジャンドル陪関数 Plm = legendre(l, cos(theta));
syms Pl(x) Plm(x)
Pl(x) = 1/(2^l*factorial(l)) * diff(((x^2 - 1)^l), x, l);
Plm(x) = (-1)^abs(m) * (1-x^2)^(abs(m)/2) * diff(Pl(x),x, abs(m));

% 符号と係数
sign = (-1)^((m+abs(m))/2);
normSH = sqrt( 4*pi / (2*l+1) * factorial(l+abs(m)) / factorial(l-abs(m)) );

% 複素数 球面調和関数
Ylm = sign/normSH * Plm(cos(theta)) .* exp(1i * m * phi);

% 
Blm = 1/sqrt(l(l+1)) * (diff(

%% 球面積分で内積計算
yy = int(int( ...
        Ylm * conj(Ylm) * sin(theta), ...
    phi, 0, 2*pi), ... 
    theta, 0, pi);

double(yy)

%%
function C = calcC(n,l,p,m)

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
            
            C = 1/(eps + sqrt(jl*(jl+1) - jn*(jn+1))) ...
                * ( term1 * calcC(n,l,p,m) ...
                    + term2 * calcC(n,l + term3 );
        end
    end
    end
end

end % function makeC
