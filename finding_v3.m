clear all
clc
close all
%in Matlab:
% y direction is vertical, downwards
% x direction is horizontal, rightwards

override=false;
%override=true;

if override==true
    %file_no
    start_no=2;    end_no=254;
    %extension
    ext='.jpg';
    %fancy save frames
    movie=false;    movie_plot_max=0.15;
    %LP filtering
    edges_LP=false;    spread_LP=false;
    %five point derivative?
    five=false;
else
    
    path=input('Path to files [current]: ');
    if isempty(path)
        path = pwd;
    end
    cd(path);
    
    ext=input('Path to files (in apostrophes!) [.jpeg]: ');
    if isempty(ext)
        ext = '.jpeg';
    end

    start_no=input('Number of first image [1]: ');
    if isempty(start_no)
        start_no = 1;
    end
    end_no=input('Number of last image [200]: ');
    if isempty(end_no)
        end_no = 200;
    end
    edge_img=input('Number of edge picking image [100]: ');
    if isempty(edge_img)
        edge_img = 100;
    end
    length=input('Length of the box in mm [100]: ');
    if isempty(length)
        length = 100;
    end
    width=input('Width of the box in mm [100]: ');
    if isempty(width)
        width = 100;
    end
    
    reply = input('Create movie frames? Y/N [N]: ', 's');
    if isempty(reply)
        reply = 'N';
    end
    if reply=='Y'
        movie=true;
        movie_plot_max=0.15;
    else
        movie=false;    
    end
    
    reply = input('Use five-point derivative? [N]: ', 's');
    if isempty(reply)
        reply = 'N';
    end
    if reply=='Y'
        five=true;
    else
        five=false;    
    end

end



I=imread([ num2str(edge_img) ext]);
%figure(1);
%imshow(imadjust(I(:,:,3)));
imshow(I);
%set(gcf, 'Position', get(0,'ScreenSize'));
%axishandle = gca;
clear x y
[x , y] = getcornersa(4);
close(1)


firel=round((x(1)+x(2))/2);
firer=round((x(3)+x(4))/2);
tl=round((y(1)+y(2))/2); %y coordinate also necessary for further rotation
tr=round((y(3)+y(4))/2);
deg=atan((round((x(2)+x(3))/2)-round((x(1)+x(4))/2))/(round((y(2)+y(3))/2)-round((y(1)+y(4))/2))); %calculate required rotation
if round((y(2)+y(3))/2)>round((y(1)+y(4))/2) %top is below bottom, then revert
    deg=deg+pi;
end

%I2=imrotate(I,-deg*180/pi,'bilinear','crop');
%figure(1);
%imshow(I2);
xc=round(size(I,2)/2);%middle x point of image X2
yc=round(size(I,1)/2);%middle y point of image X2
R=[cos(deg), -sin(deg); sin(deg), cos(deg)]; %rotation matrix about xc,yc


%corners
Vp=[x-xc;y-yc]; %vectors to points
ImagePoints=round(R*Vp); %rotation of vectors
ImagePoints(1,:)=ImagePoints(1,:)+xc; %x-component of vectors to points
ImagePoints(2,:)=ImagePoints(2,:)+yc; %y-component

%rotate fire edges as well
Vf=[firel-xc, firer-xc; tl-yc, tr-yc];
Vf=round(R*Vf); %rotation of vectors
firel=Vf(1)+xc; firer=Vf(3)+xc;
%tl=Vf(2)+yc; tr=Vf(4)+yc; %for check / plot


%create mapping for perspective correction
RealPoints = [0 0;0 length;width length;width 0]; 
T = cp2tform(ImagePoints',RealPoints,'projective');
xlu=min(ImagePoints(1,:));
ylu=min(ImagePoints(2,:));
xrb=max(ImagePoints(1,:));
yrb=max(ImagePoints(2,:));

firel=firel-xlu;
firer=firer-xlu;

%check if everything is inside the picture
if five
 if xlu<1
    xlu=1;
end
if ylu<3
    ylu=3;
end
if xrb>size(I,2)-2
    xrb=size(I,2)-2;
end
if yrb>size(I,1)-2
    yrb=size(I,1)-2;
end
if firel<1
    firel=1;
end
if firer>xrb-3
    firer=xrb-3;
end   
else
if xlu<1
    xlu=1;
end
if ylu<2
    ylu=2;
end
if xrb>size(I,2)-2
    xrb=size(I,2)-2;
end
if yrb>size(I,1)-2
    yrb=size(I,1)-2;
end
if firel<1
    firel=1;
end
if firer>xrb-2
    firer=xrb-2;
end
end

%initialize
DDav=zeros((yrb-ylu+1),(end_no-start_no));

for k=start_no:end_no
    %open image
    I=imread([num2str(k) ext]);    
    
    %improve quality for derivatives
    J=imrotate(I(:,:,1),-deg*180/pi,'bilinear','crop');
    mi=min(J(:));
    J=J-mi;
    mx = max(J(:));
    J=double(J)/double(mx);

    %calculate derrivative in x direction
    D=zeros(yrb-ylu,xrb-xlu);
    if five
    for y=ylu:yrb
        for x=xlu:xrb
            D(y-ylu+1,x-xlu+1)=(J(y-2,x)-8*J(y-1,x)+8*J(y+1,x)-J(y+2,x))/12;
        end
    end        
    else
    for y=ylu:yrb
        for x=xlu:xrb
            D(y-ylu+1,x-xlu+1)=0.5*(J(y+1,x)-J(y-1,x));
        end
    end        
    end

    
    %average for thicker strip - 2D -> 1D
    Dav=zeros(size(D,1),1);
    for i=firel:firer
        for x=1:size(D,1)
            Dav(x)=Dav(x)+D(x,i);
        end
    end
    DDav(:,k)=Dav/(abs(firer-firel));

    %create "movie"
    if movie     
        %create plot
        x=1:length(Dav);
        figure(1)
        plot(x,Dav);
        %add timer
        time = round(k);
        axis([0 length(x)+2 -movie_plot_max movie_plot_max]);
        if time < 10
            text(0,1.1*movie_plot_max, ['00' num2str(time) ' min'], 'color', 'black', 'fontsize', 18, 'horizontalalignment', 'left')
        elseif time < 100
            text(0,1.1*movie_plot_max, ['0' num2str(time) ' min'], 'color', 'black', 'fontsize', 18, 'horizontalalignment', 'left')
        else
            text(0,1.1*movie_plot_max, [num2str(time) ' min'], 'color', 'black', 'fontsize', 18, 'horizontalalignment', 'left') 
        end    
        %combine plot and image
        frame=getframe(gcf);
        [X,map]=frame2im(frame); 
        T=[ imresize(I,size(X,1)/size(I,1)) X ];
        %show / save
        if k==start_no
            mkdir('save'); cd('save'); mkdir('vid'); 
            %mkdir('img'); mkdir('plt'); 
            cd('..');
        end
        cd('save'); cd('vid'); imwrite(T,['vid' num2str(k) ext]); cd('..'); cd('..'); %save combined image
%       cd('save'); cd('img'); imwrite(I,['img' num2str(k) ext]); cd('..'); cd('..');  %save image only
%       cd('save'); cd('plt'); imwrite(X,['plt' num2str(k) ext]); cd('..'); cd('..');  %save plot only
        close 1
    end    
    
end

DDav=DDav-min(min(DDav));
DDav=DDav/max(max(DDav));
DDav=imadjust(DDav);

imshow(DDav)
%set(gcf, 'Position', get(0,'ScreenSize'));

%axishandle = gca;
%[xc , yc, linehandle] = getpoints(axishandle);

axishandle = gca;
%leading edge
[t , pos] = getpointsa();
pos=pos+ylu;
posle=zeros(1,size(DDav,2));
tle=zeros(1,size(DDav,2));
for i=1:size(DDav,2)
    if i<round(t(1))
        posle(i)=pos(1);
    elseif i>round(t(size(t,2)))
        posle(i)=pos(size(pos,2));
    else   
        posle(i)=spline(t,pos,i);
    end
    tle(i)=i;
end

%trailing edge
clear t pos
[t , pos] = getpointsa();
pos=pos+ylu;
poste=zeros(1,size(DDav,2));
tte=zeros(1,size(DDav,2));
for i=1:size(DDav,2)
    if i<round(t(1))
        poste(i)=pos(1);
    elseif i>round(t(size(t,2)))
        poste(i)=pos(size(pos,2));
    else   
        poste(i)=spline(t,pos,i);
    end
    tte(i)=i;
end


figure();
plot(tle,posle,'r');
hold on
plot(tte,poste,'b');
legend('Leading edge','Trailing edge');
xlabel('Time [min]');
ylabel('Position [pixel]');

%get real lengths
%leading
V=zeros(size(posle,2),2);
V(:,2)=posle;
V(:,1)=(firel+firer)/2;
V=tformfwd(T,V);
posleR=V(:,2);
%trailing
V=zeros(size(posle,2),2);
V(:,2)=poste;
V(:,1)=(firel+firer)/2;
V=tformfwd(T,V);
posteR=V(:,2);


figure();
plot(tle,posleR,'r');
hold on
plot(tte,posteR,'b');
legend('Leading edge','Trailing edge');
xlabel('Time [min]');
ylabel('Position [mm]');

%calculate spread rate
S_LE=zeros(size(posle,2)-2,1);
S_TE=S_LE;
S_LER=zeros(size(posle,2)-2,1);
S_TER=S_LER;

for i=2:round(t(size(t,2)))-1%length(tle)-1 %calculate derivatives 1D
    S_LER(i-1)=0.5*(posleR(i+1)-posleR(i-1));
    S_TER(i-1)=0.5*(posteR(i+1)-posteR(i-1));
    S_LE(i-1)=-0.5*(posle(i+1)-posle(i-1));
    S_TE(i-1)=-0.5*(poste(i+1)-poste(i-1));
end
tS=1:size(S_LE,1);
tS_R=1:size(S_LER,1);
S_LE(S_LE==0)=NaN;
S_TE(S_TE==0)=NaN;
S_LER(S_LER==0)=NaN;
S_TER(S_TER==0)=NaN;

figure()
plot(tS_R,S_LER,'-r');
hold on
plot(tS_R,S_TER,'-b');
legend('Leading edge','Trailing edge');
xlabel('Time [min]');
ylabel('Spread rate [mm/min]');
