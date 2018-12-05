% ���ʒ��a�֐����m�̊p�^���ʂ̍���
% �ϕ\���̋K�񕪉��iClebsch-Gordan �W���e�[�u���̍쐬�j
% 
%
% �p�^���ʂ��K�񕪉�
%   V_j1 (X) V_j2 = V_|j1-j2| (+) ... (+) V_(j1+j2)
% 
% ���ʒ��a�֐����m�̐ρA�Ή�����v�f���Ƃ̐ςł���
% SH1.*SH2
% 
% �E�ӂ��v�Z���āA���ӂƓ����ɂȂ邩�m�F
%
% 
% �Ώ̐�
%   CG(j1,j2,J,-m1,-m2,-M) = (-1)^(j1+j2-J) * CG(j1,j2,J,m1,m2,M)
% 
% 
% ������
% �N���l�b�J�[�ρi�ЎR�����M�p����ȁj
% ���\��r�A�uX1,X2���ꂼ��L2�m�����v�uX1X2�A�_�}�[����L2�m�����v
% �u������w�摜��͂̂��߂̋��ʃe���\���㐔�v�̎����A�ȒP�ȗ�
% 
% 
% �Ȃɂ����߂�����Ă��Ȃ�
% ���a�ƒ��ρA�x�N�g����Ԃ̒u�������A�Ή��t��
% CG�W���̋�̓I�ȈӖ����m�F
% 
%% ����
% SH(l,m) = SH(l^2+l+1+m)
SH1 = [1, 1,1,1];  
SH2 = [1, 1,1,1];

%% �萔�Ƃ��̗p��
% �W�J����
B1 = sqrt(length(SH1)) - 1;
B2 = sqrt(length(SH2)) - 1;

%% �E�ӂ̌v�Z
% �e
% CG �W���̃��[��
% j,m�͐�����������
% m1+m2 = M
% |j1-j2| =< J =< j1+j2
% M =< |J|

xSH = SH1'*SH2; % �ϕ\��

%% ��̗�A��������
% V0 (X) V0 �̓W�J
% SH( 00; 00) -> CG(1,1,2, 1,1,2)
% 
% 
% V1 (X) V1 �̓W�J     SH(lm;lm) = SH(lm) (X) SH(lm)
% SH( 11; 11) -> CG(1,1,2, 1,1,2)
% 
% SH( 11; 10) -> CG(1,1,2, 1,0,1) + CG(1,1,1, 1,0,1)
% SH( 10; 11) -> CG(1,1,2, 0,1,1) + CG(1,1,1, 0,1,1)
% 
% SH( 11;1-1) -> CG(1,1,2, 1,-1,0) + CG(1,1,1, 1,-1,0) + CG(1,1,0, 1,-1,0)
% SH( 10; 10) -> CG(1,1,2, 0,0,0) + CG(1,1,1, 0,0,0) + CG(1,1,0, 0,0,0)
% SH(1-1; 11) -> CG(1,1,2, -1,1,0) + CG(1,1,2, -1,1,0) + CG(1,1,2, -1,1,0)
% 
% SH(1-1; 10) -> CG(1,1,2, -1,0,-1) + CG(1,1,1, -1,0,-1)
% SH( 10;1-1) -> CG(1,1,2, 0,-1,-1) + CG(1,1,1, 0,-1,-1)
% 
% SH(1-1;1-1) -> CG(1,1,2, -1,-1,-2)
% 

%% Clebsch-Gordan table �̌v�Z
CGcoef = zeros(size(xSH));  % Clebsch-Gordan �W���̕\

CGcoef(1,:) = 1;    % V0 �̐ϕ\���ɂ��Ă�1
CGcoef(:,1) = 1;    % 

for ij1 = 1:B1      % Vi (X) Vj �̃u���b�N�w��
    for ij2 = 1:B2  % 
        for im1 = -ij1:ij1      % �u���b�N���K�񕪉��������Ƃ̃}�X�ڎw��
            for im2 = -ij2:ij2  % 
                idx1 = ij1^2 + ij1 + im1 + 1;
                idx2 = ij2^2 + ij2 + im2 + 1;
                M = im1+im2;
                
                % SH(ij1,im1) (X) SH(ij2,im2) �ɂ��Ă̌W�����v�Z����
                S = 0;
                for iL = ij1+ij2 : -1 : max(abs(M),abs(ij1-ij2))
                    CG = ClebschGordan(ij1,ij2,iL,im1,im2,M);
                    S = S + CG;
                    disp([ij1 ij2 iL, im1 im2 M])
                    disp(CG)
                end
                CGcoef(idx1,idx2) = S;
                disp(['SH(', num2str(ij1), num2str(im1), ') (X) SH(', num2str(ij2) ,num2str(im2), ')'])
                disp([idx1,idx2, S])
                disp('---')
            end
        end
    end
end

C = xSH.*CGcoef;

%% �č\��
% �W���s������ƂɁA���ʏ���č\������


R = zeros(size(SH(0,0)));
for il1 = 1:B1      % Vi (X) Vj �̃u���b�N�w��
    for il2 = 1:B2  % 
        for im1 = -ij1:ij1      % �u���b�N���K�񕪉��������Ƃ̃}�X�ڎw��
            for im2 = -il2:il2  % 
                idx1 = il1^2 + il1 + im1 + 1;
                idx2 = il2^2 + il2 + im2 + 1;
                M = im1+im2;
                for iL = il1+il2 : -1 : max(abs(M),abs(il1-il2))    % L �̎w��
                    term = sqrt((2*il1 + 1)*(2*il2 + 1) / (4*pi*(2*iL+1)));

                    R = R ...
                        + term ...
                        * ClebschGordan(il1,il2,iL,0,0,0) ...
                        * C(idx1,idx2) ...
                        * SH(iL,M); 
                end
            end
        end
    end
end
figure, plotSH(R)

% for ij1 = 1:B1      % Vi (X) Vj �̃u���b�N�w��
%     for ij2 = 1:B2  % 
%         for im1 = -ij1:ij1      % �u���b�N���K�񕪉��������Ƃ̃}�X�ڎw��
%             for im2 = -ij2:ij2  % 
%                 R = R + C(ij1^2+ij1+1+im1,ij2^2+ij2+1+im2) ...
%                     * SH(ij1,im1).*SH(ij2,im2); 
%             end
%         end
%     end
% end

