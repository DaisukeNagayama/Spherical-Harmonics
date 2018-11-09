%%
%
%

%%
function [wt,wp] = VSHsample()
%% ����
% order,degree ���w�� (1 =< order, |degree| =< order)
B_OD = [2,-1; 8,3];   % [order1, degree1; order2, degree2; ...]
B_COEF = [1,1];   % [coef1, coef2, ...]
C_OD = [3,-3];   % [order1, degree1; order2, degree2; ...]
C_COEF = [0];   % [coef1, coef2, ...]


%% ���ʃO���b�h�̗p��
polarNum = 41;    % theta, mu   �ܓx�̍��ݐ�
azimuthNum = 41;  % phi, lambda �o�x�̍��ݐ�

% �ʏ�̈ܓx�o�x
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

% �K�E�X�ܓx�ʂƌo�x�ɂ̌v�Z
syms X;
LegPolynomial = legendreP(polarNum, X); % polarNum���̃��W�����h��������
mu_ = vpasolve(LegPolynomial == 0);      % ���W�����h���������̍����K�E�X�ܓx�ɂȂ�
mu_ = fliplr(double(mu_)');                        % �V���{���b�N����{���x���֕ϊ�
lambda_ = 0 : 2*pi / (azimuthNum - 1) : 2*pi;
[lambda,mu] = meshgrid(lambda_,mu_);


%% VSHBP �̖{��
wt_prov = zeros(polarNum,azimuthNum);
wp_prov = zeros(polarNum,azimuthNum);

for j = 1:size(B_OD,1)      % Blm
    order = B_OD(j,1);
    degree = B_OD(j,2);
    
    % N
    N = (-1)^((degree+abs(degree))/2) * sqrt( (2*order+1) / (4*pi) * ...
            factorial(order-abs(degree)) / factorial(order+abs(degree)) );
    
    % Ylm
    Plm = legendre(order,mu(:,1));
    Plm = Plm(abs(degree) + 1,:)';
    Plm = repmat(Plm,1,size(lambda,2));
    Ylm = N * Plm .* exp(1i * degree * lambda);
    
    % Blm
    A = 1/sqrt(order * (order+1));
    dtYlm = diff(Ylm,1,1) * 2*pi;
    dtYlmW = [dtYlm(1,:); dtYlm];
    Blm_t = A * dtYlmW;
    Blm_p = A * 1i * degree * Ylm./(sin(acos(mu)) + eps);
    
    % �������킹�Ă���
    wt_prov = wt_prov + B_COEF(j) * Blm_t;
    wp_prov = wp_prov + B_COEF(j) * Blm_p;
end

for j = 1:size(C_OD,1)      % Clm
    order = C_OD(j,1);
    degree = C_OD(j,2);
    
    % N
    N = (-1)^((degree+abs(degree))/2) * sqrt( (2*order+1) / (4*pi) * ...
            factorial(order-abs(degree)) / factorial(order+abs(degree)) );
    
    % Ylm
    Plm = legendre(order,mu(:,1));
    Plm = Plm(abs(degree) + 1,:)';
    Plm = repmat(Plm,1,size(lambda,2));
    Ylm = N * Plm .* exp(1i * degree * lambda);
    
    % Blm
    AA = 1/sqrt(order * (order+1));
    dtYlm = diff(Ylm,1,1) * 2*pi;
    dtYlmW = [dtYlm(1,:); dtYlm];
    Blm_t = AA * dtYlmW;
    Blm_p = AA * 1i * degree * Ylm./(sin(acos(mu)) + eps);
    
    % Clm
    Clm_t = Blm_p;
    Clm_p = -Blm_t;
    
    % �������킹�Ă���
    wt_prov = wt_prov + C_COEF(j) * Clm_t;
    wp_prov = wp_prov + C_COEF(j) * Clm_p;
end

wt = wt_prov;
wp = wp_prov;

% �č\���摜���K�E�X�ܓx�����ʈܓx�ɕϊ�����
wt = spline(acos(mu(:,1)), wt', theta(:,1))';
wp = spline(acos(mu(:,1)), wp', theta(:,1))';


return  % �ȉ��Aplot��
%% �v���b�g
wtMax = max(max(abs(wt))); wpMax = max(max(abs(wp)));
wtR = real(wt)/wtMax; wpR = real(wp)/wpMax;
wtI = imag(wt)/wtMax; wpI = imag(wp)/wpMax;

%% ���ʃx�N�g���� �̉���

figure
quiver(1:azimuthNum,1:polarNum,wtR,wpR)
xlabel('��'),ylabel('��')
axis tight
ax = gca;
ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
ax.YTick = 1:(polarNum-1)/4:polarNum;
ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
ax.FontSize = 16;


% %%% 2�����摜�ŕ\��
% 
% figure
% subplot(2,2,1)
% imshow(wtR,[-1,1],'InitialMagnification',2000);
% title('wtR')
% 
% subplot(2,2,2)
% imshow(wpR,[-1,1],'InitialMagnification',2000);
% title('wpR')
% 
% subplot(2,2,3)
% imshow(wtI,[-1,1],'InitialMagnification',2000);
% title('wtI')
% 
% subplot(2,2,4)
% imshow(wpI,[-1,1],'InitialMagnification',2000);
% title('wpI')



% %%% ���ʂŕ\��
% rad = 101;      % �x�N�g�����̔��a
% rad2 = 100;     % �n�ʂ̔��a
% 
% % ����
% wt = real(wt);      wp = real(wp);
% tau_r = rad.*sin(theta);
% tau_x = tau_r.*cos(phi);
% tau_y = tau_r.*sin(phi);
% tau_z = rad.*cos(theta);
% tau_u = -wp.*sin(phi) - wt.*cos(theta).*cos(phi);
% tau_v = wp.*cos(phi) - wt.*cos(theta).*sin(phi);
% tau_w = wt.*sin(theta);
% 
% % �n��
% r2 = rad2.*sin(theta);
% x2 = r2.*cos(phi);
% y2 = r2.*sin(phi);
% z2 = rad2.*cos(theta);
% 
% % �`��
% figure
% quiver3(tau_x,tau_y,tau_z,tau_u,tau_v,tau_w, ...
%     'MaxHeadSize',1, ...
%     'LineWidth',3, ...
%     'AutoScaleFactor',1, ...
%     'Marker','.')
% axis equal off      % set axis equal and remove axis
% view(40,30)         % set viewpoint
% camzoom(1.5)        % zoom into scene %1.5
% hold on
% surf(x2,y2,z2, ...
%     'EdgeColor','none', ...
%     'FaceLighting','gouraud');
% colormap gray
% light
% hold off




end % function VSHsample

