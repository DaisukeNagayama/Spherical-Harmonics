%% ƒVƒ“ƒ{ƒŠƒbƒN‚Å‹…–Ê’²˜a•ÏŠ·
% 
% 
%%
syms x
assume(x >= -1 & x <= 1)
syms theta phi
assume(theta >= 0 & theta <= pi)
assume(phi >= 0 & phi <= 2*pi)
% syms l m
% assume(in(l,'integer') & l >= 0)
% assume(in(m,'integer') )% & m >= -l & m <= l)
syms Pl(x) Plm(x)

%% input
yy = sin(2*phi);
N = 2;

%% ‰ñ“]s—ñ
% Rx = [1 0 0;0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
% Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
% Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

% yy = Rz * Ry * yy;

%% l,m ‚É‚Â‚¢‚Äƒ‹[ƒv
f = zeros(1,(N+1)^2);
for l = 0:N
    for m = -l:l
        %% ‹…–Ê’²˜aŠÖ”‚ÌŒvŽZ
        % ƒ‹ƒWƒƒƒ“ƒhƒ‹”†ŠÖ” Plm = legendre(l, cos(theta));       
        Pl(x) = 1/(2^l*factorial(l)) * diff(((x^2 - 1)^l), x, l);
        Plm(x) = (-1)^abs(m) * (1-x^2)^(abs(m)/2) * diff(Pl(x),x, abs(m));

        % •„†‚ÆŒW”
        sign = (-1)^((m+abs(m))/2);
        normSH = sqrt( 4*pi / (2*l+1) * factorial(l+abs(m)) / factorial(l-abs(m)) );

        % •¡‘f” ‹…–Ê’²˜aŠÖ”
        Ylm = sign/normSH * Plm(cos(theta)) .* exp(1i * m * phi);

        %% ‹…–ÊÏ•ª‚Å“àÏŒvŽZ        
        f(l^2+m+1) = double(int(int( ...
                yy * conj(Ylm) * sin(theta), ...
            phi, 0, 2*pi), ... 
            theta, 0, pi));
    end
end

f




