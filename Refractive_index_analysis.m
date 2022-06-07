clc
clear
close all
 
% taking the average intensity value by considering two reference frames
ref = imread('C:\Users\mihir\OneDrive\Desktop\Research\Results\proper data set\organized data\refraction coe\ref.JPG');
rect=[1250.5 430.5 1956 2112];
reff = imcrop(rgb2gray(ref),rect);
meanInten= mean(reff(:));
 
ref2 = imread('C:\Users\mihir\OneDrive\Desktop\Research\Results\proper data set\organized data\refraction coe\ref1.JPG');
rect2=[1366.5 474.5 1904 2148];
reff2 = imcrop(rgb2gray(ref2),rect2);
meanInten2 = mean(reff2(:));
 
MeanIntensity = (meanInten+meanInten2)/2;
 
%% importing the video
 
v=VideoReader('vid.MOV');
frame_rate = v.FrameRate;
delta_t = 1/frame_rate;  % time interval between two frames
 
no_of_frames=1000;
nof=no_of_frames;
Frames=read(v,[250 nof+249]);
initial_frame = Frames(:,:,:,1);
rect=[5.725100000000000e+02,46.510000000000000,1.109800000000000e+02,5.499800000000000e+02];
F1=imcrop(Frames(:,:,:,1),rect);
 
%% Calculating average intensity
 
tt=zeros(size(F1,1),size(F1,2),nof);
resize=zeros(size(F1,1),size(F1,2),nof);
test_field=zeros(size(F1,1),size(F1,2),nof);
tot_inten=zeros(size(F1,1),size(F1,2),nof);
for k=1:nof
    frame=Frames(:,:,:,k);
    A1(:,:,:,k)= rgb2gray(frame);
    tt(:,:,k)=imcrop(A1(:,:,:,k),rect);
%     resized(:,:,k)=imresize(tt,[75 75]);
    test_field(:,:,k)=double(tt(:,:,k));
    tot_inten(:,:,k)= tot_inten(:,:,k)+test_field(:,:,k) ;
    
end
total=0;
for i=1:nof
    total=total+tot_inten(:,:,i);
    
end
avg_inten=total/nof;
%% plotting the variation of the refractive index
intensity=imresize(avg_inten,[75 75]);
refractive_index =(intensity-MeanIntensity)./MeanIntensity;
 
figure (1)
imshow(rgb2gray(F1));
 
figure (2)
surf(refractive_index');
hold on
zlim([-50 50])
title('Variation of the Refractive index')
zlabel('(n-n_{0}/n_{0})%')
 
%% Scale
 
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
 
%% positions
dis=[2.5 3.5 4.5]*0.01;
size_vel=size(F1);
P=size_vel(1)-floor(dis/Scale);
 
n0=1.0003; %refractive index of the steady air
 
n=((avg_inten-MeanIntensity)./MeanIntensity)*n0+n0
index=imresize(n,[75 75]);
figure (3)
surf(index);
hold on
% zlim([-50 50]
title('the variation of the Refractive index')
zlabel('Refractive index(n)')
 
 
%% finding mean refractive index
 
N=zeros(4,size_vel(2));
N(4,:)= (1:size_vel(2))*Scale*100;
for kk=1:3
    n_tot=0;
    for i=P(kk)-10:P(kk)+10
        n_tot=n_tot+n(i,:);
    end
    N(kk,:)=n_tot/20;
    mean_n(kk)=sum(N(kk,:))/size_vel(2);
     
    figure (kk+5)
    plot(N(4,:),N(kk,:),'-b');
    hold on
    xlim([0 1.4]);
    hl=yline(mean_n(kk),'r');
    hl.LineWidth = 1.5;
    title(['The graph of refractive index vs Distance (at y = ',num2str(dis(kk)*100),' cm from vent)']);
    xlabel('Distance (cm)');
    ylabel('Refractive index (n)');
end
hold off
 
figure (9)
plot(N(4,:),N(1,:));
hold on 
xlim([0 1.4]);
plot(N(4,:),N(2,:));
plot(N(4,:),N(3,:));
title('The graph of Refractive index vs pixels' );
xlabel('Distance (cm)');
ylabel('refractive index (n)');
legend('y= 3.5cm','y= 4.5cm','y= 5.5cm');
 
%% error of the mean Refractive index
 
SD=0; % SD = Standard deviation
SEr=zeros(1,3);
for k3=1:3
    SD(k3)=sqrt(sum((N(k3,:)-mean_n(k3)).^2)/size_vel(2));     
    SEr(k3)=sqrt(sum((N(k3,:)-mean_n(k3)).^2)/size_vel(2))/sqrt(size_vel(2))
end
 
SE= round(SEr,2) 
mean_Refractive=round(mean_n,2)
 
%%
figure (11)
x=[2.5 3.5 4.5];
plot(x,mean_Refractive,'-r');
hold on 
h2=yline(n0,'-k');
h2.LineWidth = 1;
errorbar(x,mean_Refractive,SE,'ob');
ylim([0.98 1.12]);
xlim([2.4 4.6])
title('the graph of refractive index vs distance from the vent');
xlabel('distance from the vent (cm)');
ylabel('refractive index (n)');
legend('n = computed refractive index ','n_{0} = 1.0003');


