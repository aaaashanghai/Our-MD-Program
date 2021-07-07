clear all
% Copyright Team Chi, Liu and Tang @ Tsinghua Univ., proposed for project
% for Computational Physics Course in 2021 Spring.

%����/��/��ԭ��
rows = 11;
colums = 12;
layers = 1;
N = rows*colums*layers ; %ԭ����Ŀ
%��x��y��z�����ϲ���PBC������PBC�ĳ���
xlength = 0.142*colums*1.5;
ylength = (0.142*2*rows)*sqrt(3)/2 ;
zlength = layers*10 ;
length = [xlength ylength zlength];

% ���ǵ�ʯīϩ��ԭ����������ԭ�ӵ����ԣ��Ҿ�����2*N������ԭ����Ŀ
% ��һ�����ʵı�ʾ�ַ����ʸ������ж���
% ͬʱ�Ҵ����Բ�ͬ����ĸ����ͬԭ�ӣ�������˵��
% mass = zeros(2*N+2*N1,1),2*N1��������ʾ�������ԭ�ӣ������Ļ��о����У�

mass = zeros(2*N,1); %���������ӵ�����
coord = zeros(2*N,3); %���������ӵ�λ��
veloc = zeros(2*N,3); %���������ӵ��ٶ�
accel = zeros(2*N,3); %���������ӵļ��ٶ�
tmax = 5*10^3 ; %���ʱ�䲽��
timestep = 0.0000001 ; %ʱ�䲽
parttype = zeros(2*N,1); %�����洢��������
rcut = 0.2; %���ӽضϰ뾶���뾶֮�ⲻ����������;0611���£���������¡�
POTType = 1; %�������࣬1ΪLJ
ALGOType = 1; %�㷨���࣬1ΪVerlet

%���³�ʼ��
if ALGOType == 1
    coordprev = zeros (2*N,1); %Verlet�㷨��Ҫ�洢��һʱ�̵�λ�ã���Ҫ�����µ�����
end

for i=1:2*N
    mass(i)=1.993*10^(-26); %�������������Ը�����������parttype����ֵ
end

%%
%�����Ǹ�����ʼ���еĵط�
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
        %����ʯīϩ���Ǿ���ṹ����������ƽ��ṹ����λ��nm
        %�ҽ�2*N��Ϊ��ԭ����������Ϊʯīϩһ��ԭ���ڲ�������ԭ��
        %Χ��ʹ�ü򵥵�rand()������������ƽ���ں�z����
        %�����ñ�ķ����Լ����Χ��
        %������Ȼ�����Լ����Ż������������ʼԭ�Ӽ��Ϊ����ɶ��
end
scatter3(coord(:,1),coord(:,2),coord(:,3),10,'filled')
view(0,90)
%%

for i=1:2*N
    for j=1:3
        veloc(i,j)=0; %������ٶȣ������ñ�ķ���
    end    
end    

if ALGOType == 1
    for i=1:2*N
        for j=1:3
            coordprev(i,j) = coord(i,j) - timestep*veloc(i,j); %Verlet�㷨������һ��λ�ñ�ʾ���ٶ�
        end
    end
end
%%
%2021/06/19���£���������Ҫ���Ǵν��ڵ�ԭ�ӣ����Դ˴������ν���
%2021/06/20���£��ν��ھ���shit
%�������LJ���������Ǹ����������V0������λ��r0������Ӧ��AB��
%V = Ar^(-12)-Br^(-6),���볤�ȵ�λΪnm������������λΪJ
%���ܵĸ����ο�Ħ������
r0 = 0.142;
V0 = 2000*10^3/(6.02*10^23);

% B = V0*54/29*r0^6;
% A = B*27*29/1462*r0^6;

B = V0*2*r0^6;
A = B*0.5*r0^6;

%%
%�߽�������������������ֻ��x��������,lambdaΪ�������ܶȣ���λN/m
lambda = 10^(-9);
fx_per_atom =lambda*0.142*sqrt(3)*10^(-9); 

%%
for t = 1 : tmax
    %���¼�����������
    accel = zeros(2*N,3); %��ռ��ٶ�

    if POTType == 1 %LJ����
        [sxi,xmin] = sort(coord(:,1),'ascend');
        [sxa,xmax] = sort(coord(:,1),'descend');
        [syi,ymin] = sort(coord(:,2),'ascend');
        [sya,ymax] = sort(coord(:,2),'descend');
        for i = 1:2*N
            for j=i+1:2*N %����ţ�ٵ������ɽ�Լһ�������;
                %�����ж������Ƿ�С��rcut
                if (abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut)
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %����PBC��i��j��λʸ
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
                            % ������һ����PBC������i��j�ľ���
                            dist = norm(distvec);
                            if dist < rcut %����Ҫ�ж�һ���ǲ����ڽضϰ뾶��
                                %��������ʱ�򻹿����ж�һ��ԭ�������پ���������������ʱ��ȥ������� 
                                absforce = -12*A*(1/dist)^12+6*B*(1/dist)^6; %���Ĵ�С����Ϊ�ų⣬��Ϊ����
                                force = distvec * absforce / dist; %i�ܵ�������j����ʸ��
                                accel (i,:) = accel(i,:)+transpose(force/mass(i));
                                accel (j,:) = accel(j,:)-transpose(force/mass(j));
                            end
                        end
                    end
                end        
            end
%             %��һ��z�����Լ������Ȼ����ò�ƺ������ܣ�
%             kz = 10^(-12);
%             accel (i,:) = accel(i,:)+transpose(-kz*coord(i,3)*[0;0;1]/mass(i));
            %�����һ���߽�����
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

    %���¸�����ռ�״̬

    if ALGOType == 1 %Verlet
        coordnext = zeros(2*N,3); %�����ݴ���һ�ε�����
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
%���gif�ļ�
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
