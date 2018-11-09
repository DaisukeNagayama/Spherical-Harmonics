%% �x�N�g�����ʒ��a�֐� Clm �� �v�Z����
% VSH_B.m �� VSH_C.m �̓Z�b�g
% 
% "Computational Harmonic Analysis for Tensor Fields on The Two-Sphere"
% 
%%
function Clm = VSH_C(l,m)
%% ���ʃO���b�h�̌�_�̐�
polarNum = 41;      % �k�ɂ����� (0�`��)
azimuthNum = 41;    % ��� (0�`2��)

%% ���ʃO���b�h�̍쐬
theta_ = 0 : pi/(polarNum - 1) : pi;     % polar angle
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi; % azimuth angle
[phi,theta] = meshgrid(phi_,theta_);      % define the grid

%% �X�J���[���ʒ��a�֐�
Ylm = zeros(polarNum, azimuthNum,3);
for k = 1:3
    kl = l + k - 2;     % kl = l-1, l, l+1
    
    % ���W�����h�����֐��̗p��
    Plm = legendre(kl, cos(theta(:,1))); % �e�s�� Pl0 ���� Pll �܂ł����񂾍s��
    Plm = Plm(abs(m) + 1,:)';           % Plm �̍s�����o���ďc�ɂ���
    Plm = repmat(Plm,1,size(phi,2));    % ���ʃO���b�h�̃T�C�Y�ɍ��킹��

    % �����Ɛ��K��
    sign = (-1)^((m+abs(m))/2);
    normSH = sqrt( 4*pi / (2*kl+1) * factorial(kl+abs(m)) / factorial(kl-abs(m)) );

    % ���f ���ʒ��a�֐�
    Ylm(:,:,k) = sign/normSH * Plm .* exp(1i * m * phi);
end

%% �x�N�g�����ʒ��a�֐�
term0 = (sqrt(l*(l+1)) * sin(theta) + eps).^-1;
term1 = (l+1)*(l+m) / (2*l+1);
term2 = 1i*m;
term3 = -l*(l-m+1) / (2*l+1);

Clm = zeros(polarNum, azimuthNum,2);
Clm(:,:,1) = term0 .* (term2 * Ylm(:,:,2));                         % �Ɛ���
Clm(:,:,2) = term0 .* (term1 * Ylm(:,:,1) + term3 * Ylm(:,:,3));    % �Ӑ���


end
