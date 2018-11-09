%% 球面調和関数の回転不変性を可視化する
%
%
%   テスト用
%       yy = rotSphere(SHreal(1,1),0,45); SHRotInv(yy,3);
%       yy = rotSphere(SHreal(1,1)+SHreal(2,-1),20,-30); SHRotInv(yy,3);
%       SHRotInv(2,-1,10,30,3);
%       SHRotInv(3,-2,10,30,3);
%
%       SHRotInv(2,-1,10,30,2);
%       SHRotInv(1,0, 20,0,2);
%
%%
function SHRotInv(yy,orderMax)

coef = SHT(yy,orderMax);

figure
    plotSH(yy,0.75)

figure
for il = 0:orderMax
    for im = -il:il
        x = orderMax + im + 1;
        y = il + 1;
        xy = x + (y - 1) * (2*orderMax+1);
        % xy = il + im + (2*il + 1) * orderMax + 1;   % 上の３行をまとめた式
        subplot(orderMax+1,2*orderMax+1,xy)
        plotSH(SHreal(il,im),abs(coef(y,x)))
    end
end

end