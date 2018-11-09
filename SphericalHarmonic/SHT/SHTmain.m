% ���ʒ��a�֐���p������]�s�ςȓ����ʒ��o
% 
% �X�J���[�������ʒ��a�֐��������B
% 

% CONSTANTS
L_max = 3;
RES = [55,55];

tic

% prepare INPUT FIELD
field1 = SHTinput();

field2 = ones(RES(1),RES(2));
field3 = SH(0,0)/5 + SH(1,1)/7 + SH(2,1) + SH(3,1)/2 + SH(7,1)/10;
field4 = SH(0,0)/5 + SH(1,-1)/7 + SH(2,-1) + SH(3,-1)/2 + SH(7,-1)/10;
field5 = rotSphere(SH(2,-1,RES),60,0);

field = field5;

% ���ʒ��a�ϊ�(SHT)
coefficients = SHT(field,RES,L_max);


% norm calculation

% spectrum = sum(abs(coef),2);

% ����l ���ƂɃm�����v�Z
spectrum = zeros(1,L_max+1);
for il = 1 : L_max+1
    COEF = zeros(size(coef));
    COEF(il,:) = coef(il,:);
    gg1 = SHBT(COEF);
    % plotSH(gg1,1)     % �e�������Ƃɍč\�������Ƃ��̋��ʔg
    spectrum(il) = SHnorm(gg1);
end

toc

%% �v���b�g
plotSH(field,1)
plotCoef(coef)
spectrum'
% real(coef), imag(coef)

% figure
% bar3(spectrum)
% zlabel('L2-Norm'),ylabel('����l')
% ax = gca;
% ax.YTickLabel = num2str([0:orderMax]');

