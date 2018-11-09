function spectrum = TSHfeature()
%% 3�����{�N�Z���f�[�^
%% ��]���s�ړ��s�ςȓ����ʂ̌v�Z
%   3d���h���ϊ�
% �� ���a����fft ����U�������𒊏o
% �� �e���\�����ʒ��a�ϊ�
% �� �������Ƃɑ������킹�Ă���L2-norm ���v�Z
% 
% input : voxel(x,y,z) 3�����{�N�Z���摜
%     ���_����̋���Z
%     y������z����x�������̉�]��(0=<��<��)
%     z������x����y�������̉�]��(0=<��<2��)
%     
% output: ������()
% 
%
%% ���摜�̗p��
% Model 1
% voxel = zeros(16,16,16);
% voxel(1:16,7:10,7:10) = 1;
% voxel(7:10,1:16,7:10) = 1;
% voxel(7:10,7:10,1:16) = 1;

% Model 2
% voxel = zeros(16,16,16);
% voxel(1:16,7:10,7:10) = 1;
% voxel(7:10,1:16,7:10) = 1;
% voxel(7:10,7:10,1:16) = 1;
% voxel(7:10,1:4,1:16) = 1;

% Model 3
voxel = zeros(16,16,16);
voxel(1:16,7:10,7:10) = 1;
voxel(7:10,1:16,7:10) = 10;
voxel(7:10,7:10,1:16) = 1;
voxel(7:10,1:4,1:16) = 1;

% voxel = zeros(32,32,32);
% voxel(3:29,3:29,16) = 1;

%% ���e�p�x�̐ݒ�
theta = linspace(0,pi,8);
phi   = linspace(0,2*pi,8);

%% ��]���Ă��͂ݏo���Ȃ��悤��0���ߊg���ƈʒu����
oldSize = size(voxel);
newSize = ceil(sqrt(3)*size(voxel));
voxel(newSize(1),newSize(2),newSize(3)) = 0;
voxel = circshift(voxel, ceil(oldSize/2));

%% ��], ���s�ړ�
% ��]
% voxel = rot3d(voxel,0,30);
% ���s�ړ�
% voxel = circshift(voxel,[2,3,-2]);

% ���f���̕\��
% plotVoxel(voxel)

%% �ꎟ�����e�̍쐬�i�R�������h���ϊ��j
% Proj(Z,theta,phi)
Proj = projection3d(voxel, theta, phi);

%% ���a�����Ƀt�[���G�ϊ� (���s�ړ��s�ςɂ���)
fProj = fft(Proj,[],1);
fProj = abs(fProj/size(fProj,1));           % �p���[(�U��)�X�y�N�g���𒊏o
fProj = fProj(1:ceil(size(fProj,1)/2),:,:); % ���E�Ώ̂Ȃ̂ŕБ��݂̂ɂ���

%%% �p���[�X�y�N�g����r�i��=��/2 �����̎��o���j
% figure
% bar3(permute(fProj( :,ceil(size(fProj,2)/2),: ),[1,3,2]))
% title('�� = ��/2')
% xlabel('��'),ylabel('���g��'),zlabel('�p���[�X�y�N�g��')
% ax = gca;
% ax.YTick = [];
% ax.XTick = [1,ceil(size(fProj,3)/2),size(fProj,3)]';
% ax.XTickLabel = char({'0', '��', '2��'});

% figure
% for i = 1:49
% %     bar3(permute( Proj(i,:,:),[2,3,1] ));           % ��-�ӕ\��
% %     bar3(permute( Proj(:,mod(i,9)+1,:),[1,3,2] ));  % Z-�ӕ\��
%     bar3(permute( Proj(:,:,mod(i,16)+1),[1,2,3] )); % Z-�ƕ\��
%     axis([0 inf 0 inf 0 200])
%     drawnow
%     pause(0.2)
% end



%% �e���\�����ʒ��a�ϊ�(TSHT)
n = nextpow2(size(fProj,1));
N = 2*n;
fn = TSHT(permute(fProj,[2,3,1]),N);

%% fn �� m�ɂ��Ă̑��a���v�Z�i��]�s�ςɂ���j
spectrum = permute(sum(abs(fn),2),[1,3,2]);

%%% �C��
% spectrum = spectrum(:,1:16);
% spectrum(2:2:6,:) = 100*spectrum(2:2:6,:);

%% �v���b�g

figure
bar3(spectrum)
xlabel('���'),ylabel('����'),zlabel('�U��')
ax = gca;
ax.YTickLabel = num2str([abs(n):N]');










