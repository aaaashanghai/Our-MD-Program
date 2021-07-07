clear all
% Copyright Team Chi, Liu and Tang @ Tsinghua Univ., proposed for project
% for Computational Physics Course in 2021 Spring.

%几行/列/层原胞
rows = 11;
colums = 12;
layers = 1;
N = rows*colums*layers ; %原胞数目
%在x、y、z方向上采用PBC，这是PBC的长度
xlength = 0.142*colums*1.5;
ylength = (0.142*2*rows)*sqrt(3)/2 ;
zlength = layers*10 ;
length = [xlength ylength zlength];

% 考虑到石墨烯的原胞含有两个原子的特性，我觉得用2*N来代表原胞数目
% 是一个合适的表示手法，故更改下列定义
% 同时我打算以不同的字母代表不同原子，具体来说：
% mass = zeros(2*N+2*N1,1),2*N1就用来表示多的杂质原子，这样的话感觉还行？

mass = zeros(2*N,1); %用来存粒子的质量
coord = zeros(2*N,3); %用来存粒子的位置
veloc = zeros(2*N,3); %用来存粒子的速度
accel = zeros(2*N,3); %用来存粒子的加速度
tmax = 5*10^3 ; %最大时间步数
timestep = 0.0000001 ; %时间步
parttype = zeros(2*N,1); %用来存储粒子种类
rcut = 0.2; %粒子截断半径，半径之外不计算作用力;0611更新，这里改了下。
POTType = 1; %势能种类，1为LJ
ALGOType = 1; %算法种类，1为Verlet

%以下初始化
if ALGOType == 1
    coordprev = zeros (2*N,1); %Verlet算法需要存储上一时刻的位置，需要定义新的数组
end

for i=1:2*N
    mass(i)=1.993*10^(-26); %赋予质量，可以根据粒子种类parttype来赋值
end

%%
%这里是给出初始排列的地方
imperf = 0.001;
for i=1:N
    if rem(i,2)==0
        coord(2*i-1,1)=imperf*rand() + 0.426/2*rem(floor((i-1)),colums)-0.071;
        coord(2*i-1,2)=imperf*rand() + 0.142*sqrt(3)*(rem(floor((i-1)/colums),colums));
        coord(2*i-1,3)=imperf*rand();
        coord(2*i,1)=imperf*rand() + 0.426/2*rem(floor((i-1)),colums)+0.071;
        coord(2*i,2)=imperf*rand() + 0.142*sqrt(3)*(rem(floor((i-1)/colums),colums));
        coord(2*i,3)=imperf*rand();
    else
        coord(2*i-1,1)=imperf*rand() + 0.426/2*(rem(floor((i-1)),colums))-0.071;
        coord(2*i-1,2)=imperf*rand() + 0.142*sqrt(3)*(rem(floor((i-1)/colums),colums)+0.5);
        coord(2*i-1,3)=imperf*rand();
        coord(2*i,1)=imperf*rand() + 0.426/2*(rem(floor((i-1)),colums))+0.071;
        coord(2*i,2)=imperf*rand() + 0.142*sqrt(3)*(rem(floor((i-1)/colums),colums)+0.5);
        coord(2*i,3)=imperf*rand();
    end
        %赋予石墨烯六角晶格结构，给出的是平面结构，单位是nm
        %我将2*N改为了原胞数量，因为石墨烯一个原胞内部有两个原子
        %围绕使用简单的rand()来给出，包括平面内和z方向
        %可以用别的方法以及别的围绕
        %这里显然还可以继续优化，比如给定初始原子间距为参数啥的
end
scatter3(coord(:,1),coord(:,2),coord(:,3),10,'filled')
view(0,90)
%%

for i=1:2*N
    for j=1:3
        veloc(i,j)=0; %赋予初速度，可以用别的方法
    end    
end    

if ALGOType == 1
    for i=1:2*N
        for j=1:3
            coordprev(i,j) = coord(i,j) - timestep*veloc(i,j); %Verlet算法中用上一步位置表示初速度
        end
    end
end
%%
%2021/06/19更新，我们至少要考虑次近邻的原子，所以此处给出次近邻
%2021/06/20更新，次近邻就是shit
%这里给出LJ参数，我们给定势阱深度V0和势阱位置r0给出对应的AB，
%V = Ar^(-12)-Br^(-6),距离长度单位为nm，势阱能量单位为J
%势能的给定参考摩尔键能
r0 = 0.142;
V0 = 2000*10^3/(6.02*10^23);

% B = V0*54/29*r0^6;
% A = B*27*29/1462*r0^6;

B = V0*2*r0^6;
A = B*0.5*r0^6;

%%
%边界条件拉力给出，现在只管x方向拉力,lambda为力的线密度，单位N/m
lambda = 10^(-9);
fx_per_atom =lambda*0.142*sqrt(3)*10^(-9); 

%%
for t = 1 : tmax
    %以下计算力与势能
    accel = zeros(2*N,3); %清空加速度

    if POTType == 1 %LJ势能
        [sxi,xmin] = sort(coord(:,1),'ascend');
        [sxa,xmax] = sort(coord(:,1),'descend');
        [syi,ymin] = sort(coord(:,2),'ascend');
        [sya,ymax] = sort(coord(:,2),'descend');
        for i = 1:2*N
            for j=i+1:2*N %利用牛顿第三定律节约一半计算量;
                %首先判定距离是否小于rcut
                if (abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut)
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %计算PBC下i到j的位矢
                            distvec=zeros(3,1);
                            for k=1:3
                                if abs(coord(i,k)-coord(j,k))<rcut
                                    distvec(k)=coord(j,k)-coord(i,k);
                                else
                                    if coord(i,k)<coord(j,k)
                                        distvec(k)=-(length(k)+coord(i,k)-coord(j,k));
                                    else
                                        distvec(k)=length(k)+coord(j,k)-coord(i,k);
                                    end
                                end
                            end                                                   
                            % 下面这一行是PBC下粒子i与j的距离
                            dist = norm(distvec);
                            if dist < rcut %还是要判断一下是不是在截断半径内
                                %计算力的时候还可以判断一下原子类型再决定参数，这里暂时略去这个过程 
                                absforce = -12*A*(1/dist)^12+6*B*(1/dist)^6; %力的大小，负为排斥，正为吸引
                                force = distvec * absforce / dist; %i受到的来自j的力矢量
                                accel (i,:) = accel(i,:)+transpose(force/mass(i));
                                accel (j,:) = accel(j,:)-transpose(force/mass(j));
                            end
                        end
                    end
                end        
            end
%             %给一个z方向的约束，不然粒子貌似很容易跑？
%             kz = 10^(-12);
%             accel (i,:) = accel(i,:)+transpose(-kz*coord(i,3)*[0;0;1]/mass(i));
            %这里给一个边界条件
            xminjud = ismember(i,xmin(1:colums));
            xmaxjud = ismember(i,xmax(1:colums));
            yminjud = ismember(i,ymin(1:rows));
            ymaxjud = ismember(i,ymax(1:rows)); 
            if xminjud
                accel(i,:)= accel(i,:)-[1,0,0]*fx_per_atom/mass(i);
            elseif xmaxjud
                accel(i,:)= accel(i,:)+[1,0,0]*fx_per_atom/mass(i);
            end
        end
    end

    %以下更新相空间状态

    if ALGOType == 1 %Verlet
        coordnext = zeros(2*N,3); %用来暂存下一次的坐标
        coordnext = 2*coord - coordprev + accel * timestep^2;
        coordprev = coord;
        coord = coordnext;
%         for i = 1:2*N
%             for j=1:3
%                 if coord(i,j)>(length(j))
%                     coord(i,j)=coord(i,j)-length(j);
%                 else
%                     if coord(i,j)<-(length(j)/2)
%                         coord(i,j)=coord(i,j)+length(j);
%                     end
%                 end
%             end
%         end
    end
    if rem(t,10)==0
        disp(t)
        scatter3(coord(:,1),coord(:,2),coord(:,3),10,'filled')
%         axis([-xlength/2 xlength/2 -ylength/2 ylength/2 -0.1 0.1]);
        view(0,90)
        mov(t)=getframe;
    end
end

%%
%输出gif文件
outputname = 'xf_10m9.gif';
for i = 1:1:tmax/10
    I = frame2im(mov(10*i));
    [I,map] = rgb2ind(I,256);
    if i == 1
        imwrite(I,map,outputname,'gif','Loopcount',inf,'DelayTime',0.005);
    else
        imwrite(I,map,outputname,'gif','WriteMode','append','DelayTime',0.005);
    end
end
