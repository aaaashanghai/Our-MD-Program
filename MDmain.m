clear all
% Copyright Team Chi, Liu and Tang @ Tsinghua Univ., proposed for project
% for Computational Physics Course in 2021 Spring.

%ͬһ��Ŀ¼���б�����������ο���
%���Ҫ�޸�LJ���ܲ���������ѭ���뺯��accelcorr��Ҫ�޸�

%���µ�λ��ΪSI�����Կ������н��е�λ�任
N = 125 ; %������Ŀ
%��x��y��z�����ϲ���PBC������PBC�ĳ��ȣ�
%���ϣ���߽�ΪOBC������Ӧ�����Ϊ�ϴ����ֵ���ɣ�����Ӱ���Ӿ������ģ�飩
xlength = 10 ;
ylength = 10 ;
zlength = 10 ;
length = [xlength ylength zlength];
mass = zeros(N,1); %���������ӵ�����
coord = zeros(N,3); %���������ӵ�λ��
veloc = zeros(N,3); %���������ӵ��ٶ�
accel = zeros(N,3); %���������ӵļ��ٶ�
tmax = 1000000 ; %���ʱ�䲽��
timestep = 0.0001 ; %ʱ�䲽
parttype = zeros(N,1); %�����洢��������
rcut = 10 ; %���ӽضϰ뾶���뾶֮�ⲻ����������
POTType = 1; %�������࣬1ΪLJ,2ΪMorse
ALGOType = 2; %�㷨���࣬1ΪVerlet��2ΪPC4��3ΪPC6

%���³�ʼ��
if ALGOType == 1
    coordprev = zeros (N,3); %Verlet�㷨��Ҫ�洢��һʱ�̵�λ�ã���Ҫ�����µ�����
end

if ALGOType == 2
    phasespace = zeros (N,3,4); %PC4�㷨�洢λ�á��ٶȡ����ٶ�����ٶȵ���
end

if ALGOType == 3
    phasespace = zeros (N,3,6); %PC6�㷨�洢λ����λ�õ���׵���
end

for i=1:N
    mass(i)=1; %�������������Ը�����������parttype����ֵ
end

%��������£����ϣ�������͹۵��˶������Խ��������������һ��ϵ���Ĵ�/��С
for i=1:N
    coord(i,1)=rand()*0.1 + 2*floor((i-1)/25);
    coord(i,2)=rand()*0.1 + 2*rem(floor((i-1)/5),5);
    coord(i,3)=rand()*0.1 + 2*rem(i,5);%����λ�ã������ñ�ķ���
end

if ALGOType == 2 || ALGOType == 3
   phasespace(:,:,1)=coord; 
end

for i=1:N
    for j=1:3
        veloc(i,j)=0; %������ٶȣ������ñ�ķ���
    end    
end

if ALGOType == 2 || ALGOType == 3
   phasespace(:,:,2)=veloc; 
end

if ALGOType == 1
    for i=1:N
        for j=1:3
            coordprev(i,j) = coord(i,j) - timestep*veloc(i,j); %Verlet�㷨������һ��λ�ñ�ʾ���ٶ�
        end
    end
end

%������ѭ��
for t = 1 : tmax
    %���¼�����������
    %�����������ܵ�һ����PBC���ж��Ƿ���ڣ����ѭ����һ�ּ�С����ŷ�Ͼ�������Ĳ���
    accel = zeros(N,3); %��ռ��ٶ�
    
    if POTType == 1 %LJ����
        for i = 1:N          
            for j=i+1:N %����ţ�ٵ������ɽ�Լһ�������
                %�����ж������Ƿ�С��rcut
                if abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %����PBC��i��j��λʸ
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
                            % ������һ����PBC������i��j�ľ���
                            dist = norm(distvec);
                            if dist < rcut %����Ҫ�ж�һ���ǲ����ڽضϰ뾶��
                                %��������ʱ�򻹿����ж�һ��ԭ�������پ���������������ʱ��ȥ�������
                                sigma = 2^(5/6) ; %���ȡ�Ĳ���
                                epsilon = 1;
                                absforce = -24*epsilon*(2/dist*(sigma/dist)^12-1/dist*(sigma/dist)^6); %���Ĵ�С����Ϊ�ų⣬��Ϊ����
                                force = distvec * absforce / dist; %i�ܵ�������j����ʸ��
                                accel (i,:) = accel(i,:)+transpose(force/mass(i));
                                accel (j,:) = accel(j,:)-transpose(force/mass(j));
                            end
                        end
                    end
                end        
            end        
        end
    end
    
    if POTType == 2 %Morse����
        for i = 1:N          
            for j=i+1:N %����ţ�ٵ������ɽ�Լһ�������
                %�����ж������Ƿ�С��rcut
                if abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %����PBC��i��j��λʸ
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
                            % ������һ����PBC������i��j�ľ���
                            dist = norm(distvec);
                            if dist < rcut %����Ҫ�ж�һ���ǲ����ڽضϰ뾶��
                                %��������ʱ�򻹿����ж�һ��ԭ�������پ���������������ʱ��ȥ�������
                                a = 1 ; %���ȡ�Ĳ���
                                D = 1 ; %���ȡ�Ĳ���
                                r0 = 1 ;
                                absforce = 2*a*d*(exp(-2*a*(dist-r0))-exp(-a*(dist-r0))); %���Ĵ�С����Ϊ�ų⣬��Ϊ����
                                force = distvec * absforce / dist; %i�ܵ�������j����ʸ��
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
    

    %���¸�����ռ�״̬

    if ALGOType == 1 %Verlet
        coordnext = zeros(N,3); %�����ݴ���һ�ε�����
        coordnext = 2*coord - coordprev + accel * timestep^2;
        coordprev = coord;
        coord = coordnext;
        %�����Ƕ������Ա߽������Ĵ���
        %���ʹ������moduleһ��PBC��������߱���˵����ʱ�ٶ�*ʱ�䲽�Ѿ�̫��
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
        phasenext = zeros(N,3,4); %�����ݴ���һ�ε���ռ�����
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
        %�����Ƕ������Ա߽������Ĵ���
        %���ʹ������moduleһ��PBC��������߱���˵����ʱ�ٶ�*ʱ�䲽�Ѿ�̫��
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
        phasenext = zeros(N,3,6); %�����ݴ���һ�ε���ռ�����
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
        %�����Ƕ������Ա߽������Ĵ���
        %���ʹ������moduleһ��PBC��������߱���˵����ʱ�ٶ�*ʱ�䲽�Ѿ�̫��
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
    
    %�����Ӿ������ģ�飬��ǰ�趨Ϊÿ100��ʱ�䲽���һ��ͼ��
    %Ϊ�����Ч�ʽ���ɾ�����
    t
    if rem(t,1)==0
        scatter3(coord(:,1),coord(:,2),coord(:,3),1)
        view(45,45)
        mov(t)=getframe;
    end
end

%��������Ϊ�ಽ����PC4��PC6������Ҫ�������һ������������μ��������ǲ�ϣ���ı��м����
function accelcorr=Calcacceltemp(coordtemp,mass,N,POTType,rcut,xlength,ylength,zlength)
    coord=coordtemp;
    accel = zeros(N,3); %��ռ��ٶ�

    if POTType == 1 %LJ����
        for i = 1:N          
            for j=i+1:N %����ţ�ٵ������ɽ�Լһ�������
                %�����ж������Ƿ�С��rcut
                if abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %����PBC��i��j��λʸ
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
                            % ������һ����PBC������i��j�ľ���
                            dist = norm(distvec);
                            if dist < rcut %����Ҫ�ж�һ���ǲ����ڽضϰ뾶��
                                %��������ʱ�򻹿����ж�һ��ԭ�������پ���������������ʱ��ȥ�������
                                sigma = 2^(5/6) ; %���ȡ�Ĳ���
                                epsilon = 1;
                                absforce = -24*epsilon*(2/dist*(sigma/dist)^12-1/dist*(sigma/dist)^6); %���Ĵ�С����Ϊ�ų⣬��Ϊ����
                                force = distvec * absforce / dist; %i�ܵ�������j����ʸ��
                                accel (i,:) = accel(i,:)+transpose(force/mass(i));
                                accel (j,:) = accel(j,:)-transpose(force/mass(j));
                            end
                        end
                    end
                end        
            end        
        end
    end
    if POTType == 2 %Morse����
        for i = 1:N          
            for j=i+1:N %����ţ�ٵ������ɽ�Լһ�������
                %�����ж������Ƿ�С��rcut
                if abs(coord(i,1)-coord(j,1))<rcut || abs(coord(i,1)-coord(j,1)) > xlength - rcut
                    if abs(coord(i,2)-coord(j,2))<rcut || abs(coord(i,2)-coord(j,2)) > ylength - rcut
                        if abs(coord(i,3)-coord(j,3))<rcut || abs(coord(i,3)-coord(j,3)) > zlength - rcut
                            %����PBC��i��j��λʸ
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
                            % ������һ����PBC������i��j�ľ���
                            dist = norm(distvec);
                            if dist < rcut %����Ҫ�ж�һ���ǲ����ڽضϰ뾶��
                                %��������ʱ�򻹿����ж�һ��ԭ�������پ���������������ʱ��ȥ�������
                                a = 1 ; %���ȡ�Ĳ���
                                D = 1 ; %���ȡ�Ĳ���
                                r0 = 1 ;
                                absforce = 2*a*d*(exp(-2*a*(dist-r0))-exp(-a*(dist-r0))); %���Ĵ�С����Ϊ�ų⣬��Ϊ����
                                force = distvec * absforce / dist; %i�ܵ�������j����ʸ��
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