function tool_cpm_plot( r_int,r_ext,phi_int,phi_ext,n_th,label_opt,cur_time_i,par_dt )
%TOOL_CPM_PLOT Summary of this function goes here
%   Detailed explanation goes here
[N_int,M]=size(phi_int);
[N_ext,~]=size(phi_ext);
tP_mx=zeros(M,16);
dth=2*pi/(n_th-1);
theta=0:dth:2*pi;
x_int=zeros(N_int,n_th);
y_int=zeros(N_int,n_th);
x_ext=zeros(N_ext,n_th);
y_ext=zeros(N_ext,n_th);
for j=1:n_th
    [tP_mx(:,j),~]...
        =Calc_LP(M-1,cos(theta(j)));
    for i=1:N_int
        x_int(i,j)=r_int(i)*cos(theta(j));
        y_int(i,j)=r_int(i)*sin(theta(j));
    end
    for i=1:N_ext
        x_ext(i,j)=r_ext(i)*cos(theta(j));
        y_ext(i,j)=r_ext(i)*sin(theta(j));
    end
end
data_int=phi_int*(tP_mx);
data_ext=phi_ext*(tP_mx);

figure;
hold on
%[C,h]=contour(x,y,data_r,[0.9991,0.9994,0.9997,0.9999,0.99999,1.00001,1.0001,1.0003,1.0006,1.0009],'b');
%[C1,h1]=contour(x_int,y_int,data_int,[-0.0001,-7e-5,-4e-5,4e-5,7e-5,0.0001],'b');
%[C3,h3]=contour(x_int,y_int,data_int_lt,[-0.0001,-7e-5,-4e-5,4e-5,7e-5,0.0001],'r--');
%[C2,h2]=contour(x_ext,y_ext,data_ext,[-0.0003,-0.0002,-0.0001,0.0001,0.0002,0.0003],'b');
%[C4,h4]=contour(x_ext,y_ext,data_ext_lt,[-0.0003,-0.0002,-0.0001,0.0001,0.0002,0.0003],'r--');
[C1,h1]=contour(x_int,y_int,data_int,[1.01,1.001,0.999,0.99],'b');
[C2,h2]=contour(x_ext,y_ext,data_ext,[1.01,1.001,0.999,0.99],'b');

plot(x_int(N_int,:),y_int(N_int,:),'b','LineWidth',1.0);
h1.LineWidth=1;
h2.LineWidth=1;

%contour(x,y,data_r,1100,'b');
xlabel('x');
ylabel('y');
title(['t=' num2str(cur_time_i*par_dt)]);
if label_opt==1
    %clabel(C1,h1,'FontSize',20);
    %clabel(C2,h2,'FontSize',20);
    %clabel(C3,h3,'FontSize',20);
    %clabel(C4,h4,'FontSize',20);
    clabel(C1,h1,'manual','FontSize',16);
    clabel(C2,h2,'manual','FontSize',16);
end
jpg_name=strcat('fig_t_c_',num2str(cur_time_i),'.fig');
saveas(gcf,jpg_name);
hold off
drawnow

%colormap jet;
%caxis([lr,ur]);
%view(0,0);
end

