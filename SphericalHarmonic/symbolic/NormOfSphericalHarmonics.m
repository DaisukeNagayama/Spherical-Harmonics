%% ���ʒ��a�֐��̃m�����v�Z�i���όv�Z�j
% l,m �̒l����v������̓��m�̂Ƃ��� 1 �ɂȂ�A����ȊO�̂Ƃ��� 0 �ɂȂ� (���K���𐫁j
% 

syms theta phi

%% ���ʒ��a�֐��̗p��
Y00 = 1/2*sqrt(1/pi);

Y1m1 = 1/2*sqrt(3/2/pi)*sin(theta)*exp(-1i*phi);
Y10  = 1/2*sqrt(3/pi)*cos(theta);
Y11  = -1/2*sqrt(3/2/pi)*sin(theta)*exp(1i*phi);

Y2m2 = 1/4*sqrt(15/2/pi)*sin(theta)^2*exp(-2i*phi);
Y2m1 = 1/2*sqrt(15/2/pi)*sin(theta)*cos(theta)*exp(-1i*phi);
Y20 = 1/4*sqrt(5/pi)*(3*cos(theta)^2 - 1);
Y21 = -1/2*sqrt(15/2/pi)*sin(theta)*cos(theta)*exp(1i*phi);
Y22 = 1/4*sqrt(15/2/pi)*sin(theta)^2*exp(2i*phi);

%% ���ʐϕ��œ��όv�Z
YY = int(int( ...
        Y10 * conj(Y10) * sin(theta), ...
    phi, 0, 2*pi), ... 
    theta, 0, pi);

double(YY)