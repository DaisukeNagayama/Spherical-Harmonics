%% �x�N�g�����ʒ��a�ϊ� Vector Spherical Harmonics Transform
% ���ʃx�N�g������x�N�g�����ʒ��a�֐��ɋ����W�J����2�̌W����𓾂�
% 
% �菇�F SHT�g���Ă��Ƃ͑����Z
% 
% input: v_��,v_�� �x�N�g����̃ƕ��������ƃӕ�������  
%   ���� l = N �܂ŋ��߂�Ƃ���
%   
%   g_�� = <v_��/sin(��), Ylm>   % �ƕ����ɂ��Ă̌W����
%   g_�� = <v_��/sin(��), Ylm>   % �ӕ����ɂ��Ă̌W����
%   
%   �e l,m �ɂ���(l = 1~N)
%     fB, fC <-- g_��, g_��
% output: fB, fC
% 
%% �e�X�g���s��
%  [vt,vp] = VSHsample(); VSHT(vt,vp,10);
%  
%%  
% TODO: VSHBT �͖��Ȃ��݂����Ȃ̂ł��������ԈႦ�Ă���\������
%       SHT �ɓ���� 1/sin(��) ���΂����炻����ۂ��Ȃ��� �� �Ȃ��H
%       �ʐ����}�C�i�X���Ƃ��܂������Ȃ��̂ł�����ւ�H
%       �������Ⴂ�قǑ傫�Ȓl�ɂȂ��Ă��܂�
%
%% 
function [fB,fC] = VSHT(vt,vp,orderMax)
%% �T���v�����O�_�̐���ǂݍ���
polarNum = size(vt,1);    % ���W�����h���������̎����i�K�E�X�ܓx�����_�Ƃ�̂��j
azimuthNum = size(vt,2);  % �o�x�����ɉ��_�Ƃ邩(�n�_(0)�ƏI�_(2��)�𗼕��J�E���g)

%% �ʏ�̈ܓx�o�x
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

%% �K�E�X�ܓx�ʂƌo�x�ɂ̌v�Z
syms X;
LegPolynomial = legendreP(polarNum, X);     % polarNum���̃��W�����h��������
mu_ = vpasolve(LegPolynomial == 0);          % ���W�����h���������̍����K�E�X�ܓx�ɂȂ�
mu_ = fliplr(double(mu_)');                    % �V���{���b�N����{���x���֕ϊ�          
lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;    % �o�x�̍��ݕ�
[lambda,mu] = meshgrid(lambda_,mu_);

%% Gauss-Legendre �ϕ��Ɏg���K�E�X�d�݂̌v�Z
d_LegPolynomial = diff(LegPolynomial, X);
weight = double((((1-mu_.^2) .* (subs(d_LegPolynomial, mu_)).^2)/2).^-1);

%% ���͋��ʔg����ʈܓx����K�E�X�ܓx�ɕϊ�����
vt = spline(theta(:,1), vt', acos(mu(:,1)))';
vp = spline(theta(:,1), vp', acos(mu(:,1)))';


%% SHT(Spherical Harmonics Transform algorithm)
orderMaxPlus = orderMax + 1;
gt = zeros(orderMaxPlus + 1, 2 * orderMaxPlus + 1);
gp = zeros(orderMaxPlus + 1, 2 * orderMaxPlus + 1);

for order = 0 : orderMaxPlus
    % Plm(��j) ���W�����h�����֐�
    Plms = legendre(order, mu(:,1)); % order���Ƃɂ܂Ƃ߂Đ���
    for degree = -order : order
        % degree �ɑΉ����� Plm �����o��
        Plm = Plms(abs(degree) + 1,:);
        
        % ��`�ϕ� �������ɐϕ�
        GT = vt./(sin(acos(mu)) + eps) .* exp(-1i * degree * lambda);
        GT = trapz(lambda_, GT, 2)';
        GP = vp./(sin(acos(mu)) + eps) .* exp(-1i * degree * lambda);
        GP = trapz(lambda_, GP, 2)';

        % �����Ɛ��K��
        sign = (-1)^((degree+abs(degree))/2);
        normSH = sqrt( (4*pi)/(2*order+1) * factorial(order+abs(degree)) / factorial(order-abs(degree)) );
        
        
        % �K�E�X�E���W�����h���ϕ��̎��s ��(B.44)
        gt(order+1, degree+orderMaxPlus+1) = sign/normSH * sum(weight .* GT .* Plm);
        gp(order+1, degree+orderMaxPlus+1) = sign/normSH * sum(weight .* GP .* Plm);
    end
end

%% fB, fC �̌v�Z
[DEGREE,ORDER] = meshgrid(-orderMaxPlus:orderMaxPlus, 0:orderMaxPlus);

% AA = ones(size(DEGREE));
% BB = ones(size(DEGREE));
% CC = ones(size(DEGREE));
AA = sqrt(ORDER.*(ORDER+1)).^-1;
BB = ORDER.*(ORDER-DEGREE+1) ./ (2*ORDER+1);
CC = (ORDER+1).*(ORDER+DEGREE) ./ (2*ORDER+1);

fB = AA.*( ...
    BB .* circshift(gt,[-1,0]) ...
    - CC .* circshift(gt,[1,0]) ...
    - 1i * ORDER .* gp);
fC = AA.*( ...
    -BB .* circshift(gp,[-1,0]) ...
    + CC .* circshift(gp,[1,0]) ...
    - 1i .* ORDER .* gt);

%% fB,fC �𐮌`����
fB = fB(1:orderMaxPlus, 2:2*orderMaxPlus);
fC = fC(1:orderMaxPlus, 2:2*orderMaxPlus);

outOfRange = zeros(orderMax+1,2*orderMax+1) + [ ... % �͈͊O���w�肷�邽�߂̍s��
    [flipud(tril(ones(orderMax))); zeros(1,orderMax)], ...
    zeros(orderMax+1,1), ...
    [triu(ones(orderMax)); zeros(1,orderMax)]   ];
outOfRange(1,:) = 1;

fB(1,:) = 0;
fB(outOfRange == 1) = 0;
fC(1,:) = 0;
fC(outOfRange == 1) = 0;

end % function SHT()





