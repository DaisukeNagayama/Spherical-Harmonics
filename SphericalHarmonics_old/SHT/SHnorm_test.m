


%% l,m 総当たりテスト
% orderMax = 5;
% testResult = zeros(orderMax+1,2*orderMax+1);
% tic
% for il = 0:orderMax
%     for im = -il:il
%         % testResult(il+1,im+orderMax+1) = SHnorm(SH(il,im));
%         testResult(il+1,im+orderMax+1) = SHnorm(rotSphere(SH(il,im),0,30));
%     end
% end
% toc

%% 回転テスト
yy = SH(2,1);
rotEve = 0:30:180;
rotAzi = 0:30:180;

NE = length(rotEve);
NA = length(rotAzi);
testResult = zeros(NE,NA);
tic
for ie = 1:NE
    for ia = 1:NA
        testResult(ie,ia) = SHnorm(rotSphere(yy,rotEve(ie),rotAzi(ia)));
    end
end
toc

%% testResult 出力
% figure 
% subplot(1,3,1), surf(real(testResult))
% subplot(1,3,2), surf(imag(testResult))
% subplot(1,3,3), surf(abs(testResult))

figure, surf(abs(testResult))
