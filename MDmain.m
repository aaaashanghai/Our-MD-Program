clear all
% Copyright Team Chi, Liu and Tang @ Tsinghua Univ., proposed for project
% for Computational Physics Course in 2021 Spring.

%同一个目录下有变量表，方便二次开发
%如果要修改LJ势能参数，在主循环与函数accelcorr都要修改

%以下单位均为SI，可以考虑自行进行单位变换
N = 125 ; %粒子数目
%在x、y、z方向上采用PBC，这是PBC的长度；
%如果希望边界为OBC，将对应方向改为较大的数值即可（可能影响视觉化输出模块）
xlength = 10 ;
ylength = 10 ;
zlength = 10 ;
length = [xlength ylength zlength];
mass = zeros(N,1); %用来存粒子的质量
coord = zeros(N,3); %用来存粒子的位置
veloc = zeros(N,3); %用来存粒子的速度
accel = zeros(N,3); %用来存粒子的加速度
tmax = 1000000 ; %最大时间步数
timestep = 0.0001 ; %时间步
parttype = zeros(N,1); %用来存储粒子种类
rcut = 10 ; %粒子截断半径，半径之外不计算作用力
POTType = 1; %势能种类，1为LJ,2为Morse
ALGOType = 2; %算法种类，1为Verlet，2为PC4，3为PC6

%以下初始化
if ALGOType == 1
    coordprev = zeros (N,3); %Verlet算法需要存储上一时刻的位置，需要定义新的数组
end

if ALGOType == 2
    phasespace = zeros (N,3,4); %PC4算法存储位置、速度、加速度与加速度导数
end

if ALGOType == 3
    phasespace = zeros (N,3,6); %PC6算法存储位置至位置的五阶导数
end

for i=1:N
    mass(i)=1; %赋予质量，可以根据粒子种类parttype来赋值
end

%这组参数下，如果希望看到客观的运动，可以将以下三个坐标第一项系数改大/改小
for i=1:N
    coord(i,1)=rand()*0.1 + 2*floor((i-1)/25);
    coord(i,2)=rand()*0.1 + 2*rem(floor((i-1)/5),5);
    coord(i,3)=rand()*0.1 + 2*rem(i,5);%赋予位置，可以用别的方法
end

if ALGOType == 2 || ALGOType == 3
   phasespace(:,:,1)=coord; 
end

for i=1:N
    for j=1:3
        veloc(i,j)=0; %赋予初速度，可以用别的方法
    end    
end

if ALGOType == 2 || ALGOType == 3
   phasespace(:,:,2)=veloc; 
end

if ALGOType == 1
    for i=1:N
        for j=1:3
            coordprev(i,j) = coord(i,j) - timestep*veloc(i,j); %Verlet算法中用上一步位置表示初速度
        end
    end
end

%以下主循环
for t = 1 : tmax
    %以下计算力与势能
    %计算力和势能第一步是PBC下判断是否近邻，多层循环是一种减小计算欧氏距离次数的策略
    accel = zeros(N,3); %清空加速度
    
    if POTType == 1 %LJ势能
        for i = 1:N          
            for j=i+1:N %利用牛顿第三定律节约一半计算量
                %首先判定距离是否小于rcut
                if abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %计算PBC下i到j的位矢
                            distvec=zeros(3,1);
                            for k=1:3
                                if abs(coord(i,k)-coord(j,k))<rcut
                                    distvec(k)=coord(j,k)-coord(i,k);
                                else
                                    if coord(i,k)<coord(j,k)
                                        distvec(k)=-length(k)-coord(i,k)+coord(j,k);
                                    else
                                        distvec(k)=length(k)+coord(j,k)-coord(i,k);
                                    end
                                end
                            end                                                   
                            % 下面这一行是PBC下粒子i与j的距离
                            dist = norm(distvec);
                            if dist < rcut %还是要判断一下是不是在截断半径内
                                %计算力的时候还可以判断一下原子类型再决定参数，这里暂时略去这个过程
                                sigma = 2^(5/6) ; %随便取的参数
                                epsilon = 1;
                                absforce = -24*epsilon*(2/dist*(sigma/dist)^12-1/dist*(sigma/dist)^6); %力的大小，负为排斥，正为吸引
                                force = distvec * absforce / dist; %i受到的来自j的力矢量
                                accel (i,:) = accel(i,:)+transpose(force/mass(i));
                                accel (j,:) = accel(j,:)-transpose(force/mass(j));
                            end
                        end
                    end
                end        
            end        
        end
    end
    
    if POTType == 2 %Morse势能
        for i = 1:N          
            for j=i+1:N %利用牛顿第三定律节约一半计算量
                %首先判定距离是否小于rcut
                if abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %计算PBC下i到j的位矢
                            distvec=zeros(3,1);
                            for k=1:3
                                if abs(coord(i,k)-coord(j,k))<rcut
                                    distvec(k)=coord(j,k)-coord(i,k);
                                else
                                    if coord(i,k)<coord(j,k)
                                        distvec(k)=-length(k)-coord(i,k)+coord(j,k);
                                    else
                                        distvec(k)=length(k)+coord(j,k)-coord(i,k);
                                    end
                                end
                            end                                                   
                            % 下面这一行是PBC下粒子i与j的距离
                            dist = norm(distvec);
                            if dist < rcut %还是要判断一下是不是在截断半径内
                                %计算力的时候还可以判断一下原子类型再决定参数，这里暂时略去这个过程
                                a = 1 ; %随便取的参数
                                D = 1 ; %随便取的参数
                                r0 = 1 ;
                                absforce = 2*a*d*(exp(-2*a*(dist-r0))-exp(-a*(dist-r0))); %力的大小，负为排斥，正为吸引
                                force = distvec * absforce / dist; %i受到的来自j的力矢量
                                accel (i,:) = accel(i,:)+transpose(force/mass(i));
                                accel (j,:) = accel(j,:)-transpose(force/mass(j));
                            end
                        end
                    end
                end        
            end        
        end
    end
    
    
    if ALGOType == 2 || 3
        phasespace(:,:,3)=accel*timestep^2/2;
    end
    

    %以下更新相空间状态

    if ALGOType == 1 %Verlet
        coordnext = zeros(N,3); %用来暂存下一次的坐标
        coordnext = 2*coord - coordprev + accel * timestep^2;
        coordprev = coord;
        coord = coordnext;
        %以下是对周期性边界条件的处理：
        %最多使得粒子module一倍PBC，如果更高倍则说明此时速度*时间步已经太大
        for i = 1:N
            for j=1:3
                if coord(i,j)>length(j)
                    coord(i,j)=coord(i,j)-length(j);
                else
                    if coord(i,j)<0
                        coord(i,j)=coord(i,j)+length(j);
                    end
                end
            end
        end  
    end
    
    if ALGOType == 2 % PC4 Normal Form    
        B=[1 1 1 1;0 1 2 3; 0 0 1 3; 0 0 0 1];
        phasenext = zeros(N,3,4); %用来暂存下一次的相空间坐标
        for temp=1:N
            for dimtemp=1:3               
                tempphase(:)=phasespace(temp,dimtemp,:);
                phasenext(temp,dimtemp,:)= B*transpose(tempphase);
            end
        end
        coordtemp=phasenext(:,:,1);
        accelcorr=Calcacceltemp(coordtemp,mass,N,POTType,rcut,xlength,ylength,zlength);
        diffaccel=(accelcorr-phasenext(:,:,3))*timestep^2/2;
        c=[1/6;5/6;1;1/3];
        for temp = 1:4
            phasenext(:,:,temp)=phasenext(:,:,temp)+c(temp)*diffaccel(:,:);
        end
        phasespace=phasenext;
        coord=phasespace(:,:,1);
        veloc=phasespace(:,:,2)/timestep;
        %以下是对周期性边界条件的处理：
        %最多使得粒子module一倍PBC，如果更高倍则说明此时速度*时间步已经太大
        for i = 1:N
            for j=1:3
                if coord(i,j)>length(j)
                    coord(i,j)=coord(i,j)-length(j);
                else
                    if coord(i,j)<0
                        coord(i,j)=coord(i,j)+length(j);
                    end
                end
            end
        end
        phasespace(:,:,1)=coord;
    end
    
    if ALGOType == 3 % PC6 Normal Form
        B=[1 1 1 1 1 1;0 1 2 3 4 5; 0 0 1 3 6 10; 0 0 0 1 4 10; 0 0 0 0 1 5; 0 0 0 0 0 1];
        phasenext = zeros(N,3,6); %用来暂存下一次的相空间坐标
        for temp=1:N
            for dimtemp=1:3               
                tempphase(:)=phasespace(temp,dimtemp,:);
                phasenext(temp,dimtemp,:)= B*transpose(tempphase);
            end
        end
        coordtemp=phasenext(:,:,1);
        accelcorr=Calcacceltemp(coordtemp,mass,N,POTType,rcut,xlength,ylength,zlength);
        diffaccel=(accelcorr-phasenext(:,:,3))*timestep^2/2;
        c=[3/20;251/360;1;11/18;1/6;1/60];
        for temp = 1:6
            phasenext(:,:,temp)=phasenext(:,:,temp)+c(temp)*diffaccel(:,:);
        end
        phasespace=phasenext;
        coord=phasespace(:,:,1);
        veloc=phasespace(:,:,2)/timestep;
        %以下是对周期性边界条件的处理：
        %最多使得粒子module一倍PBC，如果更高倍则说明此时速度*时间步已经太大
        for i = 1:N
            for j=1:3
                if coord(i,j)>length(j)
                    coord(i,j)=coord(i,j)-length(j);
                else
                    if coord(i,j)<0
                        coord(i,j)=coord(i,j)+length(j);
                    end
                end
            end
        end
        phasespace(:,:,1)=coord;
    end
    
    %以下视觉化输出模块，当前设定为每100个时间步输出一张图像，
    %为了提高效率建议删除这段
    t
    if rem(t,1)==0
        scatter3(coord(:,1),coord(:,2),coord(:,3),1)
        view(45,45)
        mov(t)=getframe;
    end
end

%以下是因为多步法（PC4与PC6）均需要额外计算一次受力，而这次计算中我们不希望改变中间变量
function accelcorr=Calcacceltemp(coordtemp,mass,N,POTType,rcut,xlength,ylength,zlength)
    coord=coordtemp;
    accel = zeros(N,3); %清空加速度

    if POTType == 1 %LJ势能
        for i = 1:N          
            for j=i+1:N %利用牛顿第三定律节约一半计算量
                %首先判定距离是否小于rcut
                if abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %计算PBC下i到j的位矢
                            distvec=zeros(3,1);
                            for k=1:3
                                if abs(coord(i,k)-coord(j,k))<rcut
                                    distvec(k)=coord(j,k)-coord(i,k);
                                else
                                    if coord(i,k)<coord(j,k)
                                        distvec(k)=-length(k)-coord(i,k)+coord(j,k);
                                    else
                                        distvec(k)=length(k)+coord(j,k)-coord(i,k);
                                    end
                                end
                            end                                                   
                            % 下面这一行是PBC下粒子i与j的距离
                            dist = norm(distvec);
                            if dist < rcut %还是要判断一下是不是在截断半径内
                                %计算力的时候还可以判断一下原子类型再决定参数，这里暂时略去这个过程
                                sigma = 2^(5/6) ; %随便取的参数
                                epsilon = 1;
                                absforce = -24*epsilon*(2/dist*(sigma/dist)^12-1/dist*(sigma/dist)^6); %力的大小，负为排斥，正为吸引
                                force = distvec * absforce / dist; %i受到的来自j的力矢量
                                accel (i,:) = accel(i,:)+transpose(force/mass(i));
                                accel (j,:) = accel(j,:)-transpose(force/mass(j));
                            end
                        end
                    end
                end        
            end        
        end
    end
    if POTType == 2 %Morse势能
        for i = 1:N          
            for j=i+1:N %利用牛顿第三定律节约一半计算量
                %首先判定距离是否小于rcut
                if abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %计算PBC下i到j的位矢
                            distvec=zeros(3,1);
                            for k=1:3
                                if abs(coord(i,k)-coord(j,k))<rcut
                                    distvec(k)=coord(j,k)-coord(i,k);
                                else
                                    if coord(i,k)<coord(j,k)
                                        distvec(k)=-length(k)-coord(i,k)+coord(j,k);
                                    else
                                        distvec(k)=length(k)+coord(j,k)-coord(i,k);
                                    end
                                end
                            end                                                   
                            % 下面这一行是PBC下粒子i与j的距离
                            dist = norm(distvec);
                            if dist < rcut %还是要判断一下是不是在截断半径内
                                %计算力的时候还可以判断一下原子类型再决定参数，这里暂时略去这个过程
                                a = 1 ; %随便取的参数
                                D = 1 ; %随便取的参数
                                r0 = 1 ;
                                absforce = 2*a*d*(exp(-2*a*(dist-r0))-exp(-a*(dist-r0))); %力的大小，负为排斥，正为吸引
                                force = distvec * absforce / dist; %i受到的来自j的力矢量
                                accel (i,:) = accel(i,:)+transpose(force/mass(i));
                                accel (j,:) = accel(j,:)-transpose(force/mass(j));
                            end
                        end
                    end
                end        
            end        
        end
    end
    
    accelcorr=accel;
    
end