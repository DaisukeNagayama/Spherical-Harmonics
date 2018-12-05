% Calculate Spherical Harmonics
% 
% 
% INPUTS: 
% 
% L     - Spherical harmonic degree, [1x1]
% M     - Spherical harmonic order,  [1x1]
% RES   - Vector of # of points to use [#Theta x #Phi points],[1x2] or [2x1] 
% 
% 
% OUTPUTS:
% 
% Ylm   - Spherical harmonics field 
% 
% �������ʒ��a�֐��� �������a�֐� �� Tesseral spherical harmonics �Ƃ��ĂԂ炵��
% 
% Nagayama 2016

function Ylm = SH(L,M,RES)

if nargin==0
  L=3;
  M=2;
  RES=[55 55];
end

if nargin<3
  RES=[25 25];
end

% prepare spherical coodinates
theta_ = linspace(0,  pi,RES(1));
phi_   = linspace(0,2*pi,RES(2));
[PHI,THETA] = meshgrid(phi_,theta_);

% calculation of Associated Legendre functoins 
Plm = legendre(L, cos(THETA(:,1))); % �e�s�� Pl0 ���� Pll �܂ł����񂾍s��
Plm = Plm(abs(M) + 1,:)';           % Plm �̍s�����o���ďc�ɂ���
Plm = repmat(Plm,1,size(PHI,2));    % ���ʃO���b�h�̃T�C�Y�ɍ��킹��

% (complex) Spherical harmonics
sign = (-1)^((M+abs(M))/2);
normSH = sqrt( 4*pi / (2*L+1) * factorial(L+abs(M)) / factorial(L-abs(M)) );
Ylm = sign/normSH * Plm .* exp(1i * M * PHI);

% (real) spherical harmonics
% if m > 0
%     Ylm = sqrt(2) * sign/normSH * Plm .* cos(abs(M) * PHI);
% elseif m < 0 
%     Ylm = sqrt(2) * sign/normSH * Plm .* sin(abs(M) * PHI);
% else
%     Ylm(:,:) = sign/normSH * Plm;
% end

end
