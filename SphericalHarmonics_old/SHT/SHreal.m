%%%%%%%%%%%%
% ���K�� ���� ���ʒ��a�֐� Ylm �� �v�Z����
% 
% �������ʒ��a�֐��� �������a�֐� �� Tesseral spherical harmonics �Ƃ��ĂԂ炵��
% 
% 

function Ylm = SHreal(l,m)

% ���ʃO���b�h�̌�_�̐�
polarNum = 41;      % �k�ɂ����� (0�`��)
azimuthNum = 41;    % ��� (0�`2��)

% ���ʃO���b�h�̍쐬
theta_ = 0 : pi/(polarNum - 1) : pi;     % polar angle
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi; % azimuth angle
[phi,theta] = meshgrid(phi_,theta_);      % define the grid

% ���W�����h�����֐��̗p��
Plm = legendre(l, cos(theta(:,1))); % �e�s�� Pl0 ���� Pll �܂ł����񂾍s��
Plm = Plm(abs(m) + 1,:)';           % Plm �̍s�����o���ďc�ɂ���
Plm = repmat(Plm,1,size(phi,2));    % ���ʃO���b�h�̃T�C�Y�ɍ��킹��

% �����Ɛ��K��
sign = (-1)^((m+abs(m))/2);
normSH = sqrt( 4*pi / (2*l+1) * factorial(l+abs(m)) / factorial(l-abs(m)) );

% ���� ���ʒ��a�֐�
if m > 0
    Ylm(:,:) = sqrt(2) * sign/normSH * Plm .* cos(abs(m) * phi);
elseif m < 0 
    Ylm(:,:) = sqrt(2) * sign/normSH * Plm .* sin(abs(m) * phi);
else
    Ylm(:,:) = sign/normSH * Plm;
end

end
