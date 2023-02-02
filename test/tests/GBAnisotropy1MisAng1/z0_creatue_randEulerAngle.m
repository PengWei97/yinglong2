% 用于生成随机取向的欧拉角
% F:\PhD\a_CopperThinFilm_AGG_2022\脚本\z0_creatue_randEulerAngle.m
close all 
clc;clear;
count = 10;
V=zeros(count,4);

for k=1:count
    V(k,1)= roundn(unifrnd(0,45),-2); % phi 
    V(k,2)= 0; %roundn(unifrnd(0,90),-2); % theta
    V(k,3)= 0; %roundn(unifrnd(0,45),-2); % Phi
    V(k,4) = 1;
end

% selectGrainID = [62 90 41 191 122 182 4 152 149 2];
% 
% for k = 1:length(selectGrainID)
%     V(selectGrainID(k)+1,1) =  90.0; %roundn(unifrnd(89,91),-2);
% end



%-----------------------------输出---------------------------------------
filename = strcat('grn_',num2str(count),'_2D_3.txt')
f=fullfile('.\',filename)
[r,c]=size(V);
fid=fopen(filename,'w');
fprintf(fid,'Texture File\n');
fprintf(fid,'File generated from MATLAB\n');
fprintf(fid,'B ');
fprintf(fid,'%d\n',r);

for i=1:r
    fprintf(fid,'   ');
    for j=1:c
        if j==c
            fprintf(fid,'%3.2f\n',V(i,j));%如果是最后一个，就换行
        else
            if V(i,j)>0 && V(i,j)<10
               fprintf(fid,' '); 
            end
            fprintf(fid,'%4.2f',V(i,j));%如果不是最后一个，就tab
            fprintf(fid,'   ');
        end
    end

end
fclose(fid);