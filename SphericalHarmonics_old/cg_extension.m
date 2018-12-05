function [C] = cg_extension(coef1,coef2)
% CG_EXTENSION  �N���u�V���E�S���_���W�����g�p���Ċ���W�J
%   [C] = CG_EXTENSION(coef1,coef2)
%       coef1: ���ʒ��a�W���P
%       coef2: ���ʒ��a�W���Q
%       C    : ����W�J�������ʂ̋��ʒ��a�W��
%   
%   coef1 �̋��ʏ�� coef2 �̋��ʏ���|�������ʂ����߂�
%   
%   Yl1m1 * Yl2m2 = sum[L,M]{
%                       ��((2l1+1)(2l2+1)/(4��(2L+1)))
%                       �ql1 0 l2 0 | L 0�r
%                       �ql1 m1 l2 m2 | L M�r
%                       YLM
%                   }

%   1/��4�� �{�ɂȂ��Ă�H
% 
% coef1 = [0 1 0; 0 0 0];
% coef2 = [0 0 1 0 0; 0 0 1 0 0; 0 1 1 0 0];
% C = cg_extension(coef1,coef2)
% SHnorm(SHBT(C))

%%
L1 = size(coef1,1)-1;  % coef1 �̍ő原�� l1
L2 = size(coef2,1)-1;  % coef2 �̍ő原�� l2

L = abs(L1 - L2):L1 + L2;  % C �̎�肤�鎟�� L
C = zeros(L(end)+1,2*L(end)+1);

for il1 = 0:L1
    for il2 = 0:L2
        for im1 = -il1:il1
            for im2 = -il2:il2
                iL = il1 + il2;
                term1 = sqrt((2*il1+1)*(2*il2+1)/(4*pi*(2*iL+1)));
                term2 = Wigner3j(il1,il2,iL,0,0,0);
                for iM = -iL:iL
                    C(iL+1,L(end)+iM+1) = C(iL+1,L(end)+iM+1)...
                        + term1 * term2 ...
                            * Wigner3j(il1,il2,iL,im1,im2,iM) ...
                            * coef1(il1+1,L1+im1+1) ...
                            * coef2(il2+1,L2+im2+1);
                end
            end
        end
    end
end



