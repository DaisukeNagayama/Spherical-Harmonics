%% �x�N�g�����ʒ��a �t�ϊ� Vector Spherical Harmonics Back Transform
% 2�̌W����ɏ]���ăx�N�g�����ʒ��a�֐��𑫂����킹�ċ��ʃx�N�g����𓾂�
% 
% �x�N�g���� {vt,vp} �����߂�
% ���f�[�^ {wt,wp}
% 
% 
% 1 =< l, |m| =< l =< orderMax �ɂ��ċ��߂�
%  Blm
%  �Ɛ��� �� Ylm �̃ƕ����ɔ����������̂���(l(l+1))�Ŋ���
%  �Ӑ��� �� Ylm �� im/sin�� ���|�������̂���(l(l+1))�Ŋ���
%  Clm
%   ���S�����P�ʃx�N�g�� -er ��Blm �̊O��
%   �Ɛ��� �� Blm �̃Ӑ���
%   �Ӑ��� �� Blm �̃Ɛ����̃}�C�i�X
% 
%
% ���s��
%  [vt,vp] = VSHsample(); [fB,fC] = VSHT(vt,vp,10); VSHBT(fB,fC);
%
%  % fB(l=2,m=1) �𐶐�����sample�Ɣ�r���� �� ����: ��v�A���Ȃ�
%  fC = zeros(11,21);
%  fB = zeros(11,21); fB(3,12) = 1;
%  VSHBT(fB,fC);
% 
%%
function [wt,wp] = VSHBT(fB,fC)
%% �ő原���̓ǂݍ���
orderMax = size(fB,1) - 1;

%% ���ݐ����w��
polarNum = 41;      % theta  ���W�����h���������̎����i�K�E�X�ܓx�����_�Ƃ�̂��j
azimuthNum = 41;    % phi  �o�x�����ɉ��_�Ƃ邩(�n�_(0)�ƏI�_(2��)�𗼕��J�E���g)

% �ʏ�̈ܓx�o�x
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

% %%% �K�E�X�ܓx�ʂƌo�x�ɂ̌v�Z %%%
% syms X;
% LegPolynomial = legendreP(polarNum, X); % polarNum���̃��W�����h��������
% mu_ = vpasolve(LegPolynomial == 0);      % ���W�����h���������̍����K�E�X�ܓx�ɂȂ�
% mu_ = fliplr(double(mu_)');                % �V���{���b�N����{���x���֕ϊ�
% lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;
% [lambda,mu] = meshgrid(lambda_,mu_);

%% VSHBT �̖{��
wt_prov = zeros(polarNum,azimuthNum);
wp_prov = zeros(polarNum,azimuthNum);
for order = 1 : orderMax
for degree = -order : order
    
    N = (-1)^((degree+abs(degree))/2) * sqrt( (2*order+1) / (4*pi) * ...
            factorial(order-abs(degree)) / factorial(order+abs(degree)) );
    
    % Ylm
    Plm = legendre(order,cos(theta(:,1)));
    Plm = Plm(abs(degree) + 1,:)';
    Plm = repmat(Plm,1,size(phi,2));
    Ylm = N * Plm .* exp(1i * degree * phi);
    
    % Blm
    AA = 1/sqrt(order * (order+1));
    dtYlm = diff(Ylm,1,1) * 2*pi;
    dtYlmW = [dtYlm(1,:); dtYlm];
    Blm_t = AA * dtYlmW;
    Blm_p = AA * 1i * degree * Ylm./(sin(theta)+eps);
    
    % Clm
    Clm_t = Blm_p;
    Clm_p = -Blm_t;
    
    % �����Ă���
    wt_prov = wt_prov + ( ...
        fB(order + 1,degree + orderMax + 1) * Blm_t ...
        + fC(order + 1,degree + orderMax + 1) * Clm_t );
    wp_prov = wp_prov + ( ...
        fB(order + 1,degree + orderMax + 1) * Blm_p ...
        + fC(order + 1,degree + orderMax + 1) * Clm_p );
end
end
wt = wt_prov;
wp = wp_prov;


% % �č\���摜���K�E�X�ܓx�����ʈܓx�ɕϊ�����
% wt = spline(pi - acos(mu(:,1)), wt', theta(:,1))';
% wp = spline(pi - acos(mu(:,1)), wp', theta(:,1))';


end % function vectorSHBT(s)


