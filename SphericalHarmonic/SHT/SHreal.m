%%%%%%%%%%%%
% ³‹K‰» ŽÀ” ‹…–Ê’²˜aŠÖ” Ylm ‚ð ŒvŽZ‚·‚é
% 
% ŽÀ”‹…–Ê’²˜aŠÖ”‚Í —§•û’²˜aŠÖ” ‚â Tesseral spherical harmonics ‚Æ‚àŒÄ‚Ô‚ç‚µ‚¢
% 
% 

function Ylm = SHreal(l,m)

% ‹…–ÊƒOƒŠƒbƒh‚ÌŒð“_‚Ì”
polarNum = 41;      % –k‹É‚©‚ç“ì‹É (0`ƒÎ)
azimuthNum = 41;    % ˆêŽü (0`2ƒÎ)

% ‹…–ÊƒOƒŠƒbƒh‚Ìì¬
theta_ = 0 : pi/(polarNum - 1) : pi;     % polar angle
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi; % azimuth angle
[phi,theta] = meshgrid(phi_,theta_);      % define the grid

% ƒ‹ƒWƒƒƒ“ƒhƒ‹”†ŠÖ”‚Ì—pˆÓ
Plm = legendre(l, cos(theta(:,1))); % Šes‚É Pl0 ‚©‚ç Pll ‚Ü‚Å‚ª•À‚ñ‚¾s—ñ
Plm = Plm(abs(m) + 1,:)';           % Plm ‚Ìs‚ðŽæ‚èo‚µ‚Äc‚É‚·‚é
Plm = repmat(Plm,1,size(phi,2));    % ‹…–ÊƒOƒŠƒbƒh‚ÌƒTƒCƒY‚É‡‚í‚¹‚é

% •„†‚Æ³‹K‰»
sign = (-1)^((m+abs(m))/2);
normSH = sqrt( 4*pi / (2*l+1) * factorial(l+abs(m)) / factorial(l-abs(m)) );

% ŽÀ” ‹…–Ê’²˜aŠÖ”
if m > 0
    Ylm(:,:) = sqrt(2) * sign/normSH * Plm .* cos(abs(m) * phi);
elseif m < 0 
    Ylm(:,:) = sqrt(2) * sign/normSH * Plm .* sin(abs(m) * phi);
else
    Ylm(:,:) = sign/normSH * Plm;
end

end
