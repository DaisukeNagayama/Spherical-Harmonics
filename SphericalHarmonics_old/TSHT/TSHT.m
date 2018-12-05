function fn = TSHT(f,N)
%% �e���\����, �����W�J�̍ő原��
%% �e���\��������ʒ��a�֐��ŋ����W�J���ăX�y�N�g���ɕϊ�����
%   input : f(��,��,base) �e���\����
%           N �����W�J�̍ő原��(l=0~N �̋��ʒ��a�֐���p���ēW�J����)
%   
%   output: fn(l,m,base) �e��ꂲ�Ƃ̋��ʒ��a�֐��W�J�̌W����
%   
%   
%   
%% ���s�e�X�g
% AAA = rotSphere(SH(1,0),0,0);
% BBB = rotSphere(SH(2,-1),0,0);
% test_data = cat(3,AAA,BBB); 
% test_fn = TSHT(test_data,4);
% plotSH(test_data(:,:,1)); plotSH(test_data(:,:,2));
% 
%% �萔�̏���

n = nextpow2(size(f,3));  % �e���\����̊K��

%% �T���v�����O�_�̐�
elevationNum = size(f,1);    % ���W�����h���������̎����i�K�E�X�ܓx�����_�Ƃ�̂��j
azimuthNum = size(f,2);  % �o�x�����ɉ��_�Ƃ邩(�n�_(0)�ƏI�_(2��)�𗼕��J�E���g)

%% �ʏ�̈ܓx�o�x
theta_ = 0 : pi/(elevationNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

%% �K�E�X�ܓx�ʂƌo�x�ɂ̌v�Z
syms X;
LegPolynomial = legendreP(elevationNum, X);     % elevationNum���̃��W�����h��������
mu_ = vpasolve(LegPolynomial == 0);          % ���W�����h���������̍����K�E�X�ܓx�ɂȂ�
mu_ = fliplr(double(mu_)');                    % �V���{���b�N����{���x���֕ϊ�          
lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;    % �o�x�̍��ݕ�
[lambda,mu] = meshgrid(lambda_,mu_);

%% Gauss-Legendre �ϕ��Ɏg���K�E�X�d�݂̌v�Z
d_LegPolynomial = diff(LegPolynomial, X);
weight = double((((1-mu_.^2) .* (subs(d_LegPolynomial, mu_)).^2)/2).^-1);


tic;
%% �X�P�[�����O�萔�̌v�Z
% C(l,p,m,n)
%     access to C(l,p,m,n)
%     C(l -abs(n)+1, p +1, m +N+1, n +1);
C = makeC(n,N);

%% �ܓx�ϊ�
%%% ���͋��ʔg����ʈܓx����K�E�X�ܓx�ɕϊ����� %%%
for base = 1:size(f,3)
    f(:,:,base) = spline(theta(:,1), f(:,:,base)', acos(mu(:,1)))';
end

%% �⏕�֐�gn�̌v�Z�i���όv�Z�j
% gn(l,m,base)
%     access to gn(l,m,base)
%     gn(l +1, m +N+abs(n)+1, base +1);
gn = zeros(N+abs(n)+1, 2*(N-abs(n))+1, 2^abs(n));
for l = 0:N+abs(n)
    % Plm(��j) ���W�����h�����֐�
    Plms = legendre(l, mu(:,1)); % l���Ƃɂ܂Ƃ߂Đ���
    
    for m = -l:l
        % degree �ɑΉ����� Plm �����o��
        Plm = Plms(abs(m) + 1,:);
        
        norm = 4*pi/(2*l+1) * factorial(l+abs(m)) / factorial(l-abs(m));
        
        meson = 1/norm * exp(-1i*m*lambda) ./(sin(acos(mu)).^abs(n) + eps);
        
        for base = 1:size(f,3)
            G = meson .* f(:,:,base);
            G = sum(G,2)' / sqrt(azimuthNum);
            
            gn(l +1, m +N+abs(n)+1,base) = sum(weight .* G .* Plm);
        end
    end
end

%% �X�y�N�g���v�Z 
% fn(l,m,base)
%     access to fn(l,m,base)
%     fn(l -abs(n)+1, m +N+1, base +1);
fn = zeros(N-abs(n)+1, 2*N+1, 2^abs(n));
for l = abs(n):N
    for m = -l:l
        fn_ = 0;
        for p = l-n:l+n
            fn_ = fn_ + conj( C(l -abs(n)+1, p +1, m +N+1, n +1) ) * gn(p +1,m +N+abs(n)+1,:);
%             fn_ = fn_ + gn(p +1,m +N+abs(n)+1,:);
        end
        fn(l -abs(n)+1, m +N+1, :) = fn_;
    end
end


end % function tenssor

%% function �X�P�[�����O�萔C�̌v�Z
function C = makeC(n,N)

l = abs(n):N;
p = 0:abs(n)+N;
m = -N:N;

% C(l,p,m,n)
C = zeros(length(l), length(p), length(m), abs(n)+1);

% n = 0 �̂Ƃ�
for jl = l
    C(jl -abs(n)+1, jl +1, :, 0 +1) = 1;
end

% n �� 1�ȏ�̂Ƃ������Ɍv�Z
for jn = 0:n-1
    for jl = l
    for jm = m
        for jp = p
            term1 = -jm * C(jl -abs(n)+1,jp +1,jm +N+1,jn +1);
            if jp-1 < p(1)
                term2 = 0;
            else
                term2 = (jp-jm)*(jp-2*jn-1)/(2*jp-1) * C(jl -abs(n)+1,jp-1 +1,jm +N+1,jn +1);
            end
            if jp+1 > p(length(p))
                term3 = 0;
            else
                term3 = -(jp+jm+1)*(jp+2+2*jn)/(2*jp+3) * C(jl -abs(n)+1,jp+1 +1,jm +N+1,jn +1);
            end
            
            C(jl -abs(n)+1,jp +1,jm +N+1,jn+1 +1) = ...
                1/(eps + sqrt(jl*(jl+1) - jn*(jn+1))) * ( term1 + term2 + term3 );
        end
    end
    end
end

end % function makeC



