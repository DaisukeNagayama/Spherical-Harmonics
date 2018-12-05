% 球面調和関数同士の角運動量の合成
% 積表現の規約分解（Clebsch-Gordan 係数テーブルの作成）
% 
%
% 角運動量を規約分解
%   V_j1 (X) V_j2 = V_|j1-j2| (+) ... (+) V_(j1+j2)
% 
% 球面調和関数同士の積、対応する要素ごとの積でいい
% SH1.*SH2
% 
% 右辺を計算して、左辺と同じになるか確認
%
% 
% 対称性
%   CG(j1,j2,J,-m1,-m2,-M) = (-1)^(j1+j2-J) * CG(j1,j2,J,m1,m2,M)
% 
% 
% ＊＊＊
% クロネッカー積（片山さんを信用するな）
% 性能比較、「X1,X2それぞれL2ノルム」「X1X2アダマール積L2ノルム」
% 「生物医学画像解析のための球面テンソル代数」の実装、簡単な例
% 
% 
% なにか解釈が足りていない
% 直和と直積、ベクトル空間の置き換え、対応付け
% CG係数の具体的な意味を確認
% 
%% 入力
% SH(l,m) = SH(l^2+l+1+m)
SH1 = [1, 1,1,1];  
SH2 = [1, 1,1,1];

%% 定数とかの用意
% 展開次数
B1 = sqrt(length(SH1)) - 1;
B2 = sqrt(length(SH2)) - 1;

%% 右辺の計算
% 各
% CG 係数のルール
% j,mは整数か半整数
% m1+m2 = M
% |j1-j2| =< J =< j1+j2
% M =< |J|

xSH = SH1'*SH2; % 積表現

%% 具体例、メモ書き
% V0 (X) V0 の展開
% SH( 00; 00) -> CG(1,1,2, 1,1,2)
% 
% 
% V1 (X) V1 の展開     SH(lm;lm) = SH(lm) (X) SH(lm)
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

%% Clebsch-Gordan table の計算
CGcoef = zeros(size(xSH));  % Clebsch-Gordan 係数の表

CGcoef(1,:) = 1;    % V0 の積表現については1
CGcoef(:,1) = 1;    % 

for ij1 = 1:B1      % Vi (X) Vj のブロック指定
    for ij2 = 1:B2  % 
        for im1 = -ij1:ij1      % ブロックを規約分解したあとのマス目指定
            for im2 = -ij2:ij2  % 
                idx1 = ij1^2 + ij1 + im1 + 1;
                idx2 = ij2^2 + ij2 + im2 + 1;
                M = im1+im2;
                
                % SH(ij1,im1) (X) SH(ij2,im2) についての係数を計算する
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

%% 再構成
% 係数行列をもとに、球面場を再構成する


R = zeros(size(SH(0,0)));
for il1 = 1:B1      % Vi (X) Vj のブロック指定
    for il2 = 1:B2  % 
        for im1 = -ij1:ij1      % ブロックを規約分解したあとのマス目指定
            for im2 = -il2:il2  % 
                idx1 = il1^2 + il1 + im1 + 1;
                idx2 = il2^2 + il2 + im2 + 1;
                M = im1+im2;
                for iL = il1+il2 : -1 : max(abs(M),abs(il1-il2))    % L の指定
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

% for ij1 = 1:B1      % Vi (X) Vj のブロック指定
%     for ij2 = 1:B2  % 
%         for im1 = -ij1:ij1      % ブロックを規約分解したあとのマス目指定
%             for im2 = -ij2:ij2  % 
%                 R = R + C(ij1^2+ij1+1+im1,ij2^2+ij2+1+im2) ...
%                     * SH(ij1,im1).*SH(ij2,im2); 
%             end
%         end
%     end
% end

