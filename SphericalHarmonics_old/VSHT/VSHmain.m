%%
% �x�N�g�����ʒ��a�ϊ�(Vector Spherical Harmonics Transform)
% 
% 
% �ʏ�̌o�x�ܓx�̋��ʏ�̐ڐ��x�N�g���������
% 
% 

%% 
orderMax = 5;

%%
% ���̃x�N�g����̗p��
% [vt,vp] = VSHsample();

fB = zeros(orderMax+1, 2*orderMax+1);
fC = zeros(orderMax+1, 2*orderMax+1);
fB(4,orderMax+3) = 1;
fC(7,orderMax-4) = 1;
[vt,vp] = VSHBT(fB,fC);


%% �X�y�N�g��plot�p�̃J���[�}�b�v��p��
    sg = (0:1/127:1)/2;
    sc = (1:-1/127:0)/2 + 0.5;
    myColorMap = sqrt([[sg; sg; sc]'; fliplr([sc; sg; sg])']);

%% �X�y�N�g���͈̔͊O���w�肷�邽�߂̔z���p��
outOfRange = zeros(orderMax+1,2*orderMax+1) + [ ...
    [flipud(tril(ones(orderMax))); zeros(1,orderMax)], ...
    zeros(orderMax+1,1), ...
    [triu(ones(orderMax)); zeros(1,orderMax)]   ];
outOfRange(1,:) = 1;


%% �X�y�N�g�����
    tic;
    [fB,fC] = VSHT(vt,vp,orderMax);
    toc

%% �X�y�N�g������č\��
    tic;
    [wt,wp] = VSHBT(fB,fC);
    toc


%% �W���������
B = fB;             BMax = max(max(abs(B)));
BR = real(B)/BMax;  BI = imag(B)/BMax;
C = fC;             CMax = max(max(abs(C)));
CR = real(C)/CMax;  CI = imag(C)/CMax;

% figure(1)
%     % subplot(1,2,1)
%     subplot(2,2,1)
%     imgBR = imshow(BR,[-1,1],'InitialMagnification',2000);
%     set(imgBR,'AlphaData',not(outOfRange));
%     title('fB') %�i�����j')
%     xlabel('�x���i�ʐ��j degree'),ylabel('���� order')
%     colormap(myColorMap)
%     colorbar
%     axis on
%     ax = gca;
%     ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
%     ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
%     ax.YTick = [1,fliplr(orderMax+1:-2:2)];
%     ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
%     ax.TickLength = [0 0];
%     ax.FontSize = 16;
% 
%     subplot(2,2,3)
%     imgBI = imshow(BI,[-1,1],'InitialMagnification',2000);
%     set(imgBI,'AlphaData',not(outOfRange));
%     title('fB�i�����j')
%     xlabel('�x���i�ʐ��j degree'),ylabel('���� order')
%     colormap(myColorMap)
%     colorbar
%     axis on
%     ax = gca;
%     ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
%     ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
%     ax.YTick = [1,fliplr(orderMax+1:-2:2)];
%     ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
%     ax.TickLength = [0 0];
%     ax.FontSize = 16;
% 
% 
%     % subplot(1,2,2)
%     subplot(2,2,2)
%     imgCR = imshow(CR,[-1,1],'InitialMagnification',2000);
%     set(imgCR,'AlphaData',not(outOfRange));
%     title('fC') %�i�����j')
%     xlabel('�x���i�ʐ��j degree'),ylabel('���� order')
%     colormap(myColorMap)
%     colorbar
%     axis on
%     ax = gca;
%     ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
%     ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
%     ax.YTick = [1,fliplr(orderMax+1:-2:2)];
%     ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
%     ax.TickLength = [0 0];
%     ax.FontSize = 16;
% 
%     subplot(2,2,4)
%     imgCI = imshow(CI,[-1,1],'InitialMagnification',2000);
%     set(imgCI,'AlphaData',not(outOfRange));
%     title('fC�i�����j')
%     xlabel('�x���i�ʐ��j degree'),ylabel('���� order')
%     colormap(myColorMap)
%     colorbar
%     axis on
%     ax = gca;
%     ax.XTick = [1:2:orderMax,orderMax+1,fliplr(2*orderMax+1:-2:orderMax+2)];
%     ax.XTickLabel = num2str([-orderMax:2:-1,0,fliplr(orderMax:-2:1)]');
%     ax.YTick = [1,fliplr(orderMax+1:-2:2)];
%     ax.YTickLabel = num2str([0,fliplr(orderMax:-2:1)]');
%     ax.TickLength = [0 0];
%     ax.FontSize = 16;


%% ���ʃx�N�g���� �̉���
wt = real(wt); wp = real(wp);
vt = real(vt); vp = real(vp);
polarNum = size(wt,1);
azimuthNum = size(wt,2);

%%% ���ʕ\��
% figure(2)
%     subplot(1,2,1)
%     quiver(1:azimuthNum,1:polarNum,wt,wp)
%     title('�č\���x�N�g����'),xlabel('��'),ylabel('��')
%     axis tight
%     ax = gca;
%     ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
%     ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
%     ax.YTick = 1:(polarNum-1)/4:polarNum;
%     ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
%     ax.FontSize = 16;
% 
%     subplot(1,2,2)
%     quiver(1:azimuthNum,1:polarNum,vt,vp)
%     title('���̃x�N�g����'),xlabel('��'),ylabel('��')
%     axis tight
%     ax = gca;
%     ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
%     ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
%     ax.YTick = 1:(polarNum-1)/4:polarNum;
%     ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
%     ax.FontSize = 16;

figure
    subplot(1,2,1)
    quiver(1:azimuthNum,1:polarNum,vt,vp)
    title('�x�N�g����P'),xlabel('��'),ylabel('��')
    axis tight
    ax = gca;
    ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
    ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
    ax.YTick = 1:(polarNum-1)/4:polarNum;
    ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
    ax.FontSize = 16;

    subplot(1,2,2)
    quiver(1:azimuthNum,1:polarNum,vt,vp)
    title('�x�N�g����Q'),xlabel('��'),ylabel('��')
    axis tight
    ax = gca;
    ax.XTick = 1:(azimuthNum-1)/8:azimuthNum;
    ax.XTickLabel = strcat(num2str([0:1/4:2]'),'\pi');
    ax.YTick = 1:(polarNum-1)/4:polarNum;
    ax.YTickLabel = strcat(num2str([0:1/4:1]'),'\pi');
    ax.FontSize = 16;
    
%% ���ʕ\��
% rad = 101;      % �x�N�g�����̔��a
% rad2 = 100;     % �n�ʂ̔��a
% 
% % ����
% tau_r = rad.*sin(theta);
% tau_x = tau_r.*cos(phi);
% tau_y = tau_r.*sin(phi);
% tau_z = rad.*cos(theta);
% tau_u = -wp.*sin(phi) - wt.*cos(theta).*cos(phi);
% tau_v = wp.*cos(phi) - wt.*cos(theta).*sin(phi);
% tau_w = wt.*sin(theta);
% 
% % �����i���f�[�^�j
% rho_r = rad.*sin(theta);
% rho_x = rho_r.*cos(phi);
% rho_y = rho_r.*sin(phi);
% rho_z = rad.*cos(theta);
% rho_u = -vp.*sin(phi) - vt.*cos(theta).*cos(phi);
% rho_v = vp.*cos(phi) - vt.*cos(theta).*sin(phi);
% rho_w = vt.*sin(theta);
% 
% % �n��
% r2 = rad2.*sin(theta);
% x2 = r2.*cos(phi);
% y2 = r2.*sin(phi);
% z2 = rad2.*cos(theta);
% 
% % �`��
% figure
% subplot(1,2,1)
% quiver3(tau_x,tau_y,tau_z,tau_u,tau_v,tau_w, ...
%     'MaxHeadSize',1, ...
%     'LineWidth',2, ...
%     'AutoScaleFactor',2, ...
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
% 
% subplot(1,2,2)
% quiver3(rho_x,rho_y,rho_z,rho_u,rho_v,rho_w, ...
%     'MaxHeadSize',1, ...
%     'LineWidth',2, ...
%     'AutoScaleFactor',2, ...
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


