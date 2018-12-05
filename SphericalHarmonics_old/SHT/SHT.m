% Spherical Harmonic Transform
% 
% INPUTS:
% 
% field - Spherical field, [RES(1) x RES(2)]
% RES   - Vector of # of points to use [#Theta x #Phi points],[1x2] or [2x1]
% L_max - Max of SHT degree L, [1x1]
% 
% OUTPUTS:
% 
% coefficients  - Spherical harmonic coefficients, [L_max+1 x 2*L_max+1]
% 
% 
% ���͔g�Ƌ��ʒ��a�֐��̓��ς��Ƃ��Ă邾�� <f,Ylm>
% 
% ���όv�Z�����邽�ߏc���ɐϕ�����K�v������
% �ܓx�����̓K�E�X�ܓx�łƂ��ăK�E�X�E���W�����h���ϕ�����
% �o�x�����͑�`�ϕ��B�����I�Ȃ̂ŗ��U�t�[���G�ϊ��ɏ��������\
% 
% <�Ή��֌W>
%  �K�E�X�ܓx��   <---> �Ɋp theta �̃R�T�C��
%  �o�x �Ɂ@�@    <---> ���ʊp phi
% 
% MATLAB�̃��W�����h�����֐��A���S���K���łƒʏ�łƂ̊֌W��
%  Plm('norm') = Plm() * N;
%  N = sqrt( (2*l+1) / (2*pi) * factorial(l-abs(m)) / factorial(l+abs(m)) );
%  
% 
% �Q�lWeb�y�[�W
%  �p��� Wikipedia
%   https://en.wikipedia.org/wiki/Spherical_harmonics
%  B.1.3�߁u���ʂ̃X�y�N�g���@�v���Ɏ�(B.41)���玮(B.44)
%   https://www.gfd-dennou.org/member/hiroki/file/BPmodel/Note_BPM.pdf
% 
%
% ���s�e�X�g
%  SHT(SH(2,-1)+SH(1,1),2)
% 
% Nagayama Daisuke, 2016

%
function coefficients = SHT(field,RES,L_max)

if ((RES(1) ~= size(field,1)) || (RES(2) ~= size(field,2)))
    disp('error(SHT): The size of FIELD and RES must be equal.');
    return
end

THETA = linspace(0,  pi,RES(1));  % Azimuthal/Longitude/Circumferential
PHI   = linspace(0,2*pi,RES(2));  % Altitude /Latitude /Elevation
[PHI,THETA] = meshgrid(PHI,THETA);

% �K�E�X�ܓx�ʂƌo�x�ɂ̌v�Z
syms X;
LegPolynomial = legendreP(RES(1), X);     % polarNum���̃��W�����h��������
mu_ = vpasolve(LegPolynomial == 0);          % ���W�����h���������̍����K�E�X�ܓx�ɂȂ�
mu_ = fliplr(double(mu_)');                     % �V���{���b�N����{���x���֕ϊ�          
lambda_ = 0 : 2*pi / (RES(2) - 1) : 2*pi;    % �o�x�̍��ݕ�
[lambda,mu] = meshgrid(lambda_,mu_);

% Gauss-Legendre �ϕ��Ɏg���K�E�X�d�݂̌v�Z
d_LegPolynomial = diff(LegPolynomial, X);
weight = double((((1-mu_.^2) .* (subs(d_LegPolynomial, mu_)).^2)/2).^-1);

% ���͋��ʔg����ʈܓx����K�E�X�ܓx�ɕϊ�����
field = spline(theta(:,1), field', acos(mu(:,1)))';


% SHT(Spherical Harmonic Transform algorithm)
% �W���i�[�p�̕ϐ��A���ʒ��a�֐��̃s���~�b�h�}�݂����Ȓ��p�񓙕ӎO�p�`��Ɋi�[
coefficients = zeros(L_max + 1, 2 * L_max + 1);

for order = 0 : L_max    % ���ʒ��a�֐���᎟���珇��
    % Plm(��j) ���W�����h�����֐�
    Plms = legendre(order, mu(:,1)); % order���Ƃɂ܂Ƃ߂Đ���
    
    for degree = -order : order % ���[�v�ŉ񂵂ČW���v�Z����        
        
        % degree �ɑΉ����� Plm �����o��
        Plm = Plms(abs(degree) + 1,:);
        
        % ��`�ϕ� �������ɐϕ�
        G = field .* exp(-1i * degree * lambda);
        G = trapz(lambda_, G, 2)';
        
        % �����Ɛ��K��
        sign = (-1)^((degree+abs(degree))/2);
        normSH = sqrt( (4*pi)/(2*order+1) * factorial(order+abs(degree)) / factorial(order-abs(degree)) );
        
        % �K�E�X�E���W�����h���ϕ� �c�����ɐϕ� 
        coefficients(order+1, degree+L_max+1) = ...
            sign/normSH * sum(weight .* Plm .* G);
        
    end
end

end % function SHT()



