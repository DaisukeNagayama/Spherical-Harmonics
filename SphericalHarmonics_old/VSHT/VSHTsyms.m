% ƒVƒ“ƒ{ƒŠƒbƒN‚ÅŒŸŽZ—ûK

syms Y P
syms x
syms Bt Bp
syms lap

syms l m
assume(in(l,'integer') & l>=0)
assume(in(m,'integer') & m>=-l & m<=l)
syms theta phi
assume(in(theta,'real') & theta>0 & theta<pi)
assume(in(phi,'real') & phi>0 & phi<2*pi)


P(x,l,m) = 1/(2^l*factorial(l)) * (1-x^2)^(abs(m)/2) * diff((x^2-1)^l,l+abs(m));
Y(l,m) = P(cos(theta),l,m)*exp(1i *m*phi);

Bt(l,m) = 1/sqrt(l*(l+1)) * diff(Y,theta);
Bp(l,m) = 1/sqrt(l*(l+1)) * 1/sin(theta)*diff(Y,phi);

Btlap(l,m) = -( diff(Bt,theta,2) + cot(theta)*diff(Bt,theta) ...
                + 1/(sin(theta))^2 * diff(Bt,phi,2) ...
                - Bt/(sin(theta))^2 ...
                ) + 2*cot(theta)/sin(theta) * diff(Bp,phi);
Bplap(l,m) = -( diff(Bp,theta,2) + cot(theta)*diff(Bp,theta) ...
                + 1/(sin(theta))^2 * diff(Bp,phi,2) ...
                - Bp/(sin(theta))^2 ...
                ) - 2*cot(theta)/sin(theta) * diff(Bt,phi);

