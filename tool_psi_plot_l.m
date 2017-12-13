function tool_psi_plot_l( r_int,r_ext,psi_int,psi_ext,psi_int_lt,psi_ext_lt,n_th,label_opt,cur_time_i,par_dt )
%TOOL_CPM_PLOT Summary of this function goes here
%   Detailed explanation goes here
[N_int,M]=size(psi_int);
[N_ext,~]=size(psi_ext);
%psi_int_lt_m=zeros(N_int,M);
%psi_ext_lt_m=zeros(N_ext,M);
%psi_int_lt_m(:,3)=psi_int_lt;
%psi_ext_lt_m(:,3)=psi_ext_lt;
dtP_mx=zeros(M,n_th);
dth=2*pi/(n_th-1);
theta=0:dth:2*pi;
 x_int=zeros(n_th,N_int);
 y_int=zeros(n_th,N_int);
 x_ext=zeros(n_th,N_ext);
 y_ext=zeros(n_th,N_ext);
for j=1:n_th
    [~,dtP_mx(:,j)]...
        =Calc_LP(M-1,cos(theta(j)));
    for i=1:N_int
         x_int(j,i)=theta(j);
         y_int(j,i)=r_int(i);
     end
     for i=1:N_ext
         x_ext(j,i)=theta(j);
         y_ext(j,i)=r_ext(i);
     end
end
data_int=psi_int*(dtP_mx*diag(-dtP_mx(2,:)));
data_ext=psi_ext*(dtP_mx*diag(-dtP_mx(2,:)));
%data_int_lt=psi_int_lt_m*(dtP_mx*diag(-dtP_mx(2,:)));
%data_ext_lt=psi_ext_lt_m*(dtP_mx*diag(-dtP_mx(2,:)));

figure;
hold on
%[C,h]=contour(x,y,data_r,[0.9991,0.9994,0.9997,0.9999,0.99999,1.00001,1.0001,1.0003,1.0006,1.0009],'b');
[C1,h1]=contour(x_int,y_int,data_int',[0.0005,-0.0005,0.001,-0.001],'b');
%[C3,h3]=contour(x_int,y_int,data_int_lt,[-0.001,0.001],'r--');
[C2,h2]=contour(x_ext,y_ext,data_ext',[0.0005,-0.0005,0.001,-0.001],'b');
%[C4,h4]=contour(x_ext,y_ext,data_ext_lt,[0.005,0.01,-0.005,-0.01],'r--');
plot(x_int(:,N_int),y_int(:,N_int),'b','LineWidth',1.0);
h1.LineWidth=1;
h2.LineWidth=1;
%h3.LineWidth=1;
%h4.LineWidth=1;
%contour(x,y,data_r,1100,'b');
xlabel('\theta');
ylabel('r');
title(['t=' num2str(cur_time_i*par_dt)]);
if label_opt==1
    %clabel(C1,h1,'FontSize',20);
    %clabel(C2,h2,'FontSize',20);
    %clabel(C3,h3,'FontSize',20);
    %clabel(C4,h4,'FontSize',20);
    clabel(C1,h1,'manual','FontSize',16);
    clabel(C2,h2,'manual','FontSize',16);
    %clabel(C3,h3,'manual','FontSize',16);
    %clabel(C4,h4,'manual','FontSize',16);
end
jpg_name=strcat('l_fig_t_psi_',num2str(cur_time_i),'.fig');
saveas(gcf,jpg_name);
hold off
drawnow

%colormap jet;
%caxis([lr,ur]);
%view(0,0);
end

