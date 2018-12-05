function SHfeature
%% 3�����{�N�Z���f�[�^
%% ��]���s�ړ��s�ςȓ����ʂ̌v�Z
%   3d���h���ϊ�
% �� ���a����fft ����U�������𒊏o
% �� �e�v�f���Ƃɋ��ʒ��a�ϊ�
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
voxel = zeros(16,16,16);
voxel(1:16,7:10,7:10) = 1;
voxel(7:10,1:16,7:10) = 1;
voxel(7:10,7:10,1:16) = 1;

% Model 2
% voxel = zeros(16,16,16);
% voxel(1:16,7:10,7:10) = 1;
% voxel(7:10,1:16,7:10) = 1;
% voxel(7:10,7:10,1:16) = 1;
% voxel(7:10,1:4,1:16) = 1;

% Model 3
% voxel = zeros(16,16,16);
% voxel(1:16,7:10,7:10) = 1;
% voxel(7:10,1:16,7:10) = 10;
% voxel(7:10,7:10,1:16) = 1;
% voxel(7:10,1:4,1:16) = 1;

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
% voxel = rot3d(voxel,pi/6,pi/6);
% ���s�ړ�
% voxel = circshift(voxel,[2,3,-2]);

% ���f���̕\��
% plotVoxel(voxel)

%% �ꎟ�����e�̍쐬�i�R�������h���ϊ��j
% Proj(Z,theta,phi)
disp('���e�v�Z')
tic
Proj = projection3d(voxel, theta, phi);
toc

%% ���a�����Ƀt�[���G�ϊ� (���s�ړ��s�ςɂ���)
% fProj(Z,theta,phi) -> fProj(theta,phi,Z)
disp('�t�[���G�ϊ�')
tic

fProj = fft(Proj,[],1);
fProj = abs(fProj/size(fProj,1));           % �p���[(�U��)�X�y�N�g���𒊏o
fProj = fProj(1:ceil(size(fProj,1)/2),:,:); % ���E�Ώ̂Ȃ̂ŕБ��݂̂ɂ���

fProj = permute(fProj,[2,3,1]);             % ���̏����̂��߂ɕ��ёւ�

fProj = fProj(:,:,1:ceil(size(fProj,3)/8));

%%% �p���[�X�y�N�g����r�i���ԑw�̎��o���j
figure
bar3(fProj( :,:,ceil(size(fProj,3)/2) ))
title('���ԑw')
xlabel('��'),ylabel('��'),zlabel('�p���[�X�y�N�g��')
ax = gca;
ax.XTick = [];
ax.YTick = [1,ceil(size(fProj,3)/2),size(fProj,3)]';
ax.YTickLabel = char({'0', '��', '2��'});

toc

%% ���ʒ��a�ϊ�(SHT)
disp('���ʒ��a�ϊ�')
tic

orderMax = 6;
fn = zeros(orderMax + 1, 2 * orderMax + 1, size(fProj,3));
for iz = 1:size(fProj,3)
    fn(:,:,iz) = SHT(fProj(:,:,iz), orderMax);
end

%% fn �� m�ɂ��Ă̑��a���v�Z�i��]�s�ςɂ���j
spectrum = zeros(orderMax+1, size(fn,3));
for iz = 1:size(fn,3)
    for il = 1 : orderMax+1
        COEF = zeros(size(fn(:,:,1)));
        COEF(il,:) = fn(il,:,iz);
        gg1 = SHBT(COEF);
        % plotSH(gg1,1)     % �e�������Ƃɍč\�������Ƃ��̋��ʔg
        spectrum(il,iz) = SHnorm(gg1);
    end
end

toc

%% �v���b�g

figure
bar3(spectrum)
xlabel('�w'),ylabel('����'),zlabel('�m����')
ax = gca;
ax.YTick = [1:2:orderMax+1]';
ax.YTickLabel = num2str([0:2:orderMax]');

