function norm = SHnorm(yy)
%% ���ʃX�J���[���L2�m���������߂�
%       input : SH  ���ʃX�J���[��
% 
%       2�悵�Ă��狅�ʂŐϕ�
%       �e�X�g�p: SHnorm(SH(2,0)), SHnorm(SH(2,1))
% 
%% �T���v�����O�_�̐���ǂݍ���
polarNum = size(yy,1);    % ���W�����h���������̎����i�K�E�X�ܓx�����_�Ƃ�̂��j
azimuthNum = size(yy,2);  % �o�x�����ɉ��_�Ƃ邩(�n�_(0)�ƏI�_(2��)�𗼕��J�E���g)

%% �ܓx�o�x
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

% �K�E�X�ܓx�ʂƌo�x�ɂ̌v�Z
syms X;
LegPolynomial = legendreP(polarNum, X);     % polarNum���̃��W�����h��������
mu_ = vpasolve(LegPolynomial == 0);          % ���W�����h���������̍����K�E�X�ܓx�ɂȂ�
mu_ = fliplr(double(mu_)');                     % �V���{���b�N����{���x���֕ϊ�          
lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;    % �o�x�̍��ݕ�
[lambda,mu] = meshgrid(lambda_,mu_);

%% Gauss-Legendre �ϕ��Ɏg���K�E�X�d�݂̌v�Z
d_LegPolynomial = diff(LegPolynomial, X);
weight = double((((1-mu_.^2) .* (subs(d_LegPolynomial, mu_)).^2)/2).^-1);

%% ���͋��ʔg����ʈܓx����K�E�X�ܓx�ɕϊ�����
yy = spline(theta(:,1), yy', acos(mu(:,1)))';

%% �m�����v�Z
% 2�悵�Ă��狅�ʂŐϕ�
G = abs(yy).^2;

% ��`���ρi�ӕ����ɐϕ��j
G = trapz(G,2)' /2/pi;

% �K�E�X�E���W�����h���ϕ��i�ƕ����ɐϕ��j
norm = sum(weight .* G);

%% ��]�s�ϐ��̃e�X�g
% 
% de = 30;
% da = 30;
% testMat = zeros(180/de + 1, 360/da + 1);
% for iEle = 0:de:180
%     for iAzi = 0:da:360
%         testMat(iEle/de+1,iAzi/da+1) = SHnorm(rotSphere(SH(1,1), iEle, iAzi));
%     end
% end
% 

end

