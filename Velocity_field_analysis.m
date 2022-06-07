clc
clear
close all
 
tic
% testing
v=VideoReader('candle.MOV');
frame_rate = v.FrameRate;
delta_t = 1/frame_rate;  % time interval between two frames
 
no_of_pairs=125;
nof=no_of_pairs;
Frames=read(v,[260 (2*nof)+259]);
initial_frame = Frames(:,:,:,1);
imshow(initial_frame)
% final_frame=Frames(:,:,:,250);
% imshow(final_frame)
rect=[393.5 72.5 226 576];
F1=imcrop(Frames(:,:,:,1),rect);
 
%for calibration uncomment following commented line and comment above lines
%under testing section
 
% v=VideoReader('vid.MOV');
% frame_rate = v.FrameRate;
% delta_t = 1/frame_rate;  % time interval between two frames
 
% no_of_pairs=500;
% nof=no_of_pairs;
% Frames=read(v,[250 (2*nof)+249]);
% initial_frame = Frames(:,:,:,1);
% % rect=[391.5 101.5 380 500];  % cropped the images using first frame
% rect = [4.745100000000000e+02,54.510000000000000,2.079800000000000e+02,5.439800000000000e+02];
% F1=imcrop(Frames(:,:,:,1),rect);
 
%% Get the scale 
 
figure (1)
imshow(initial_frame)
[x, y]=ginput(2);
diameter_h = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
[x1, y1]=ginput(2);
diameter_v = sqrt((x1(1)-x1(2))^2+(y1(1)-y1(2))^2);
 
realdiameter_h = 8.1;  %diameter measured horizontally in cm
realdiameter_v = 8.0;  %diameter measured vertically in cm
 
diameter=(diameter_h+diameter_v)/2;
real_diameter= (realdiameter_h+realdiameter_v)*0.01/2; % in meters
 
Scale=real_diameter/diameter;
 
%% computing the displacements using iterative solution
 
tot_v=0;
tot_u=0;
for k=1:nof
    frame_1=Frames(:,:,:,(2*k-1));
    frame_2=Frames(:,:,:,(2*k));
    
    % converting them to gray scale and crop
    A1= rgb2gray(frame_1);     
    A2 = rgb2gray(frame_2);
    aa1=imcrop(A1,rect);
    aa2=imcrop(A2,rect);
    
    y = size(aa1);
    a1 = double(aa1);
    a2 = double(aa2);
    
    % Generating Ex Ey and Et
    rx = [-1 1;-1 1];
    ry = [1 1;-1 -1];
 
    Ex = zeros(y(1)-1,y(2)-1);
    Ey = zeros(y(1)-1,y(2)-1);
    Et = zeros(y(1)-1,y(2)-1);
    
    for i = 1:y(1)-1
        for j = 1:y(2)-1
            Ex(i,j) = 0.25*(sum(sum(a1(i:i+1,j:j+1).*rx)) + sum(sum(a2(i:i+1,j:j+1).*rx)));
            Ey(i,j) = 0.25*(sum(sum(a1(i:i+1,j:j+1).*ry)) + sum(sum(a2(i:i+1,j:j+1).*ry)));
            Et(i,j) = 0.25*sum(sum(a2(i:i+1,j:j+1) - a1(i:i+1,j:j+1)));
        end
    end
    
    u = zeros(y(1)-1,y(2)-1);
    v = zeros(y(1)-1,y(2)-1);
    u1 = zeros(y(1)-1,y(2)-1);
    v1 = zeros(y(1)-1,y(2)-1);
    
    ax= 4+2*sqrt(2);
    bx= 4+4*sqrt(2);
    mask = [1/bx 1/ax 1/bx;1/ax 0 1/ax;1/bx 1/ax 1/bx];
    lamda = 4;
    
    for i1 = 1:10
   
        for i = 2:y(1)-3
            for j = 2:y(2)-3
                uav = sum(sum(u(i-1:i+1,j-1:j+1).*mask));
                vav = sum(sum(v(i-1:i+1,j-1:j+1).*mask));
                P = Ex(i,j)*uav + Ey(i,j)*vav + Et(i,j);
                D = lamda^2 + (Ex(i,j))^2 + (Ey(i,j))^2;
                u1(i,j) = uav - Ex(i,j)*P/D;
                v1(i,j) = vav - Ey(i,j)*P/D;
            end
        end
    
    u = u1;
    v = v1;
    end
    u_frames(:,:,k)=u;
    v_frames(:,:,k)=v;
tot_v=tot_v+v;
tot_u=tot_u+u;
   
end
toc
%% calculating average displacements
avg_dx=tot_u/nof;
avg_dy=tot_v/nof;
 
% average velocities
avg_vx=avg_dx/delta_t;
avg_vy=avg_dy/delta_t;
 
uu=avg_vx;
vv=avg_vy;
 
% Ploting velocity vectors
figure (2)
imshow(rgb2gray(F1));
hold on
quiver(uu,vv,'g');
title('Velocity vector fields')
hold off
 
%% converting velocity to m/s from pixles/s
velocity_pixels= sqrt(uu.^2+vv.^2);
velocity = Scale*velocity_pixels;
 
% get the average velocity on several regions from the vent
% dis=[ 2.5 3.5 4.5]*0.01; 
dis=[3.5 4.5 5.5]*0.01;
size_vel=size(velocity);
P=size_vel(1)-floor(dis/Scale);
 
% Averaging vector magnitudes 
 
n=11;
kernel_size = n;
ii=floor(size_vel(1)/n)*n;
jj=floor(size_vel(2)/n)*n;
 
u_avg=NaN(ii,jj);
v_avg=NaN(ii,jj);
uu2=imresize(uu,[ii,jj]);
vv2=imresize(vv,[ii,jj]);
kernel = ones(kernel_size)/(kernel_size.^2);
 
for i2=1:ii-n
    for j2=1:jj-n
        if mod(i2,n)==0 && mod(j2,n)==0
        u_avg(i2+5,j2+5)=sum(sum(uu2(i2:i2+n-1,j2:j2+n-1).*kernel));
        v_avg(i2+5,j2+5)=sum(sum(vv2(i2:i2+n-1,j2:j2+n-1).*kernel));
        end
    end
end
 
figure (3)
imshow(rgb2gray(F1));
hold on
q=quiver(u_avg,v_avg,'g');
q.AutoScaleFactor=8 ;
q.LineWidth=1
q.MaxHeadSize=0.5
title('Averaged velocity vector fields');
 
figure (4)
imshow(rgb2gray(F1));
title('Converted image to gray scale ');
 
figure (5)
imshow(F1);
title('Cropped image of the candle flame');
 
%% measured velocity values
V_measured=zeros(4,size_vel(2));
V_measured(1,:)=2.7172;
V_measured(2,:)=2.6357;
V_measured(3,:)=2.5163;
V_measured(4,:)=(2.7172+2.6357+2.5153)/3;
 
%calculation of 
V=zeros(4,size_vel(2));
V(4,:)= (1:size_vel(2))*Scale*100;
for kk=1:3
    v_tot=0;
    for i=P(kk)-10:P(kk)+10
        v_tot=v_tot+velocity(i,:);
    end
    V(kk,:)=v_tot/20;
    mean_V(kk)=sum(V(kk,:))/size_vel(2);
     
    figure (kk+5)
    plot(V(4,:),V(kk,:),'-b');
    hold on
    xlim([0 2.55]);
    hl=yline(mean_V(kk),'r');
    hl.LineWidth = 1.5;
%     plot(V(4,:),V_measured(kk,:),'--b');
    title(['The graph of velocity vs Distance (at y = ',num2str(dis(kk)*100),' cm from vent)']);
    xlabel('Distance (cm)');
    ylabel('velocity (ms^-1)');
%     legend('calculated velocity');
    
end
hold off
 
figure (9)
plot(V(4,:),V(1,:));
hold on 
xlim([0 2.55]);
plot(V(4,:),V(2,:));
plot(V(4,:),V(3,:));
% plot(V(4,:),V_measured(4,:),'--');
title('The graph of velocity vs pixels' );
xlabel('Distance (cm)');
ylabel('velocity (ms^-1)');
legend('y= 3.5cm','y= 4.5cm','y= 5.5cm');%,'measured average velocity');  
 
%% error of the mean velocity
 
SD=0;% SD = Standard deviation
SEr=zeros(1,3);
for k3=1:3
    SD(k3)=sqrt(sum((V(k3,:)-mean_V(k3)).^2)/size_vel(2));     
    SEr(k3)=sqrt(sum((V(k3,:)-mean_V(k3)).^2)/size_vel(2))/sqrt(size_vel(2))
end
 
SE= round(SEr,5) 
mean_velocity=round(mean_V,5)
 
%% position with band 
 
figure (10)
imshow(rgb2gray(F1));
title('Converted image to gray scale ');
hold on
 
hl=yline(P(1),'y','y = 3.5 cm');
hl.LineWidth = 1.5;
h2=yline(P(2),'y','y = 4.5 cm');
h2.LineWidth = 1.5;
h3=yline(P(3),'y','y = 5.5 cm');
h3.LineWidth = 1.5;
 
for kk=1:3
    h2=yline(P(kk)-10,'--y');
    h3=yline(P(kk)+10,'--y');
end
 
%% finding the linear relationship
 
%after calibartion comment this section
 
x=mean_velocity;
y=[2.7172 2.6357 2.5163 ];
e=[0.007 0.007 0.007];
 
p=polyfit(x,y,1);
v=polyval(p,x);
plot(x,v,'-r');
hold on
errorbar(x,y,e,'o');
title('Graph of Calculated velocities using air flow meter vs Computed velocities using Schlieren images')
xlabel('Computed velocities using Schlieren images (m/s)')
ylabel('Calculated velocities using air flow meter (m/s)')
hold off
 
%% calibration y=mx+c
 
intercept=ones(size_vel(1),size_vel(2))*1.0769;
Velocity_New = 478.1923*velocity+intercept;
 
for kk=1:3
    v_tot=0;
    for i=P(kk)-10:P(kk)+10
        v_tot=v_tot+Velocity_New(i,:);
    end
    V(kk,:)=v_tot/20;
    mean_V(kk)=sum(V(kk,:))/size_vel(2);
    figure (kk+5)
    plot(V(4,:),V(kk,:),'-b');
    hold on
    xlim([0 2.55]);
    hl=yline(mean_V(kk),'r');
    hl.LineWidth = 1.5;
%     plot(V(4,:),V_measured(kk,:),'--k');
    title(['The graph of velocity vs Distance (at y = ',num2str(dis(kk)*100),' cm from candle)']);
    xlabel('Distance (cm)');
    ylabel('velocity (ms^-1)');
    legend('computed velocity','average of the computed velcoity');%,'calculated velocity using flow meter');
    
end
hold off
 
figure (9)
plot(V(4,:),V(1,:));
hold on 
xlim([0 2.55]);
plot(V(4,:),V(2,:));
plot(V(4,:),V(3,:));
% plot(V(4,:),V_measured(4,:),'--k');
title('The graph of velocity vs distance' );
xlabel('Distance (cm)');
ylabel('velocity (ms^-1)');
legend('y= 3.5cm','y= 4.5cm','y= 5.5cm');%,'measured average velocity'); 
 
% uncertainty of the velocity after calibration was done
Er_velocity = 478.1923*SE
uncertainty=round(Er_velocity,2)
mean_velocity=round(mean_V,2)
 
%% observing the velocity vectors in small region
 
% rect2=[155 420 20 30]; %for air flow
rect2=[95 225 20 30]; % for candle flame
F2=rgb2gray(F1);
 
figure (10)
imshow(F2);
hold on
quiver(uu,vv,'g');
rectangle('Position',rect2);
title('Velocity vector fields')
hold off
 
time=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
 
for i3=1:5:36
    u_crop=imcrop(u_frames(:,:,i3),rect2);
    v_crop=imcrop(v_frames(:,:,i3),rect2);
   
    figure (i3+10)
    imshow(imcrop(F2,rect2));
    hold on
    quiver(u_crop,v_crop,'g');
    title(['t = ',num2str(time((i3+4)/5)),' s']);
end


