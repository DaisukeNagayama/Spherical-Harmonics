%% ���ʃx�N�g�������������
% 
% ���ʃx�N�g���� vv(�ƍ��W, �Ӎ��W, �ǂ���̐�����)
% Ex. vv(45,30,1) -> ��=45, ��= 30, e_�Ɛ���
%
% ���s�e�X�g
%   plotVSH(cat(3,VSH_B(3,-2),VSH_C(3,-2)));
%   plotVSH(VSH(3,-2),1);
%   plotVSH(VSH(9,3),2);
% 
%%
function plotVSH(vv,scale)
%% scale(�g��k��) �w�肪�������1
if nargin < 2; scale = 1; end

%% Option
% radius = 0;     % �I�t�Z�b�g���a
% alpha = 0.5;    % ���W���T�C�Y [-��,��; -��,��; -��,��]

%% �����̌v�Z�����邩�ǂ���
% isComplex = not(isreal(vv(:,:,1)));


%% ���ʃx�N�g���� �̉���
rvt = real(vv(:,:,1)); rvp = real(vv(:,:,2));
ivt = imag(vv(:,:,1)); ivp = imag(vv(:,:,2));
polarNum = size(rvt,1);
azimuthNum = size(rvt,2);
[X,Y] = meshgrid(1:azimuthNum,1:polarNum);
theta_ = 0 : pi/(polarNum - 1) : pi;
phi_ = 0 : 2*pi/(azimuthNum - 1) : 2*pi;
[phi,theta] = meshgrid(phi_,theta_);

%% ���ʕ\��
figure
    subplot(1,2,1)
    quiver(X,Y,rvt,rvp)
    title('�����x�N�g����'),xlabel('��'),ylabel('��')
    axis tight
%     ax = gca;
%     ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
%     ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
%     ax.YTick = 1:(polarNum-1)/4:polarNum;
%     ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
%     ax.FontSize = 16;

    subplot(1,2,2)
    quiver(X,Y,ivt,ivp)
    title('�����x�N�g����'),xlabel('��'),ylabel('��')
    axis tight
%     ax = gca;
%     ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
%     ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
%     ax.YTick = 1:(polarNum-1)/4:polarNum;
%     ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
%     ax.FontSize = 16;
    
%% ���ʕ\��
rad = 101;      % �x�N�g�����̔��a
rad2 = 100;     % �n�ʂ̔��a

% �����i�����j
tau_r = rad.*sin(theta);
tau_x = tau_r.*cos(phi);
tau_y = tau_r.*sin(phi);
tau_z = rad.*cos(theta);
tau_u = -rvp.*sin(phi) - rvt.*cos(theta).*cos(phi);
tau_v = rvp.*cos(phi) - rvt.*cos(theta).*sin(phi);
tau_w = rvt.*sin(theta);

% �����i�����j
rho_r = rad.*sin(theta);
rho_x = rho_r.*cos(phi);
rho_y = rho_r.*sin(phi);
rho_z = rad.*cos(theta);
rho_u = -ivp.*sin(phi) - ivt.*cos(theta).*cos(phi);
rho_v = ivp.*cos(phi) - ivt.*cos(theta).*sin(phi);
rho_w = ivt.*sin(theta);

% �n��
r2 = rad2.*sin(theta);
x2 = r2.*cos(phi);
y2 = r2.*sin(phi);
z2 = rad2.*cos(theta);

% �`��
figure
subplot(1,2,1)
quiver3(tau_x,tau_y,tau_z,tau_u,tau_v,tau_w, ...
    'MaxHeadSize',1, ...
    'LineWidth',2, ...
    'AutoScaleFactor',2, ...
    'Marker','.')
axis equal off      % set axis equal and remove axis
view(40,30)         % set viewpoint
camzoom(1.5)        % zoom into scene %1.5
hold on
surf(x2,y2,z2, ...
    'EdgeColor','none', ...
    'FaceLighting','gouraud');
colormap gray
light
hold off

subplot(1,2,2)
quiver3(rho_x,rho_y,rho_z,rho_u,rho_v,rho_w, ...
    'MaxHeadSize',1, ...
    'LineWidth',2, ...
    'AutoScaleFactor',2, ...
    'Marker','.')
axis equal off      % set axis equal and remove axis
view(40,30)         % set viewpoint
camzoom(1.5)        % zoom into scene %1.5
hold on
surf(x2,y2,z2, ...
    'EdgeColor','none', ...
    'FaceLighting','gouraud');
colormap gray
light
hold off


