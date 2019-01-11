%%%%%%%%%%created by Qi Yan%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot results
%%
% clear
% load CL_data.mat;
%% Generate random colors
rng(0,'twister');
color_setting = cell(1,N);
for i = 1:N
    color_setting{i} = rand(1,3);
end
%% Plot aggregated trajectories
% figure;
% box on;
% hold on;
% for i = 1:N
%     plot(XX(1,:,i),XX(2,:,i),'--','color',color_setting{i},'linewidth',1);
% end
% for i = 1:N
%     plot(XX_hat_DR(1,:,i),XX_hat_DR(2,:,i),'-+','color',color_setting{i},'linewidth',0.3,'markersize',4,'MarkerIndices',1:0.1*k_f:k_f);
%     plot(XX_hat_CEKF(1,:,i),XX_hat_CEKF(2,:,i),'color',color_setting{i},'linewidth',1);
%     plot(XX_hat_CEKF(1,1,i),XX_hat_CEKF(2,1,i),'x','color',color_setting{i},'markersize',12);
% %     plot(XX_hat_CEKF(1,:,i),XX_hat_TEKF(2,:,i),'-*','color',color_setting{i},'linewidth',1,'markersize',1,'MarkerIndices',1:0.1*k_f:k_f);
% %     plot(XX_hat_CUKF(1,:,i),XX_hat_CUKF(2,:,i),'-*','color',color_setting{i},'linewidth',1,'markersize',4,'MarkerIndices',1:0.1*k_f:k_f);
% %     plot(XX_hat_ICI(1,:,i),XX_hat_ICI_EKF(2,:,i),'--+','color',color_setting{i},'linewidth',1,'markersize',4,'MarkerIndices',1:0.12*k_f:k_f);
% end
% 
% str_lgd = cell(1,N);
% for i = 1:N
%     str_lgd{i} = ['Robot ',num2str(i)];
% end
% legend(str_lgd,'location','best');
% xlabel('x direction /m');
% ylabel('y direction /m');

%% Plot separate trajectories
% figure;
% row_subplot = ceil(N/4);
% for i = 1:N
%     subplot(row_subplot,4,i);
%     plot(XX(1,:,i),XX(2,:,i),'color',color_setting{i},'LineWidth',1)
%     hold on
%     plot(XX_hat_CD(1,:,i),XX_hat_CD(2,:,i),'--','color',color_setting{i},'LineWidth',1);
%     plot(XX_hat_SD(1,:,i),XX_hat_SD(2,:,i),'--o','color',color_setting{i},'linewidth',1,'markersize',2,'MarkerIndices',1:0.05*k_f:k_f);
%     plot(XX_hat_CD(1,1,i),XX_hat_CD(2,1,i),'x','color',color_setting{i},'MarkerSize',12);
%     legend(['Robot ',num2str(i)]);
% end

%% NEES result
figure; box on;
hold on;
plot(Tk,NEES_DR,'linewidth',1,'color','black');
% plot(Tk,NEES_CEKF,'linewidth',1);
plot(Tk,NEES_TEKF,'linewidth',1);
plot(Tk,NEES_SATEKF,'linewidth',1);
line([0 max(xlim)],[1 1],'linewidth',1,'color','red','linestyle','--');
ylim([0 max(NEES_CEKF)]);
% legend('DR','CEKF','SA-TEKF','TEKF','NEES = 1');
legend('DR','Dense-CL','OPT-CL','NEES = 1');
xlabel('Time/s');
ylabel('NEES');
% title('NEES');
% ylim([0 25]);
%% Absolute Error
% figure; box on;
% for i = 1:N
%     subplot(ceil(N/2),2,i);
%     box on;
%     hold on;
%     plot(Tk,abs(XX(1,:,i)-XX_hat_DR(1,:,i)),'LineWidth',1,'color','black');
%     plot(Tk,abs(XX(1,:,i)-XX_hat_CEKF(1,:,i)),'LineWidth',1);
%     plot(Tk,abs(XX(1,:,i)-XX_hat_CUKF(1,:,i)),'LineWidth',1);
%     xlabel('time /s');
%     ylabel('Abs error in X /m');
%     xlim([0 t_f]);
%         legend('DR','CEKF','CUKF','location','best');
%     legend('DR','CEKF','location','best');
%     title(['Robot ',num2str(i)]);
% end
%% RMSE
figure; box on;
for i = 1:N
    subplot(ceil(N/3),3,i);
    box on;
    hold on;
    plot(Tk,RMSE_DR(:,i),'linewidth',1,'color','black');
%     plot(Tk,RMSE_CEKF(:,i),'linewidth',1);
    plot(Tk,RMSE_TEKF(:,i),'linewidth',1);
    plot(Tk,RMSE_SATEKF(:,i),'linewidth',1);
%     legend('DR','CEKF','SA-TEKF','TEKF','location','best');
    xlabel('time/s');
    ylabel('RMS error/m');
    title(['Robot ',num2str(i)]);
end
legend('DR','Dense-CL','OPT-CL','location','best');
%% Error with 3-sigma boundaries
figure; box on;
for i = 1:N
    subplot(ceil(N/3),3,i);
    box on;
    hold on;
%     plot(Tk,+Three_sigma_CEKF(1,:,i),[color_plot,'-.'],'LineWidth',1);
%     plot(Tk,-Three_sigma_CEKF(1,:,i),[color_plot,'-.'],'LineWidth',1);
%     plot(Tk,XX(1,:,i)-XX_hat_CEKF(1,:,i),'red','LineWidth',1);
    title(['Robot ',num2str(i)]);
    
    plot(Tk,+Three_sigma_TEKF(1,:,i),['blue','--*'],'LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
    plot(Tk,-Three_sigma_TEKF(1,:,i),['blue','--*'],'LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
    plot(Tk,XX(1,:,i)-XX_hat_TEKF(1,:,i),'b-*','LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);

    plot(Tk,+Three_sigma_SATEKF(1,:,i),['red','--^'],'LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
    plot(Tk,-Three_sigma_SATEKF(1,:,i),['red','--^'],'LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
    plot(Tk,XX(1,:,i)-XX_hat_SATEKF(1,:,i),'r-*','LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
    
%     plot(Tk,+Three_sigma_CUKF(1,:,i),['blue','--*'],'LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
%     plot(Tk,-Three_sigma_CUKF(1,:,i),['blue','--*'],'LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
%     plot(Tk,XX(1,:,i)-XX_hat_CUKF(1,:,i),'b-*','LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
    
    plot(Tk,+Three_sigma_DR(1,:,i),['black','--+'],'LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
    plot(Tk,-Three_sigma_DR(1,:,i),['black','--+'],'LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
%     plot(Tk,XX(1,:,i)-XX_hat_DR(1,:,i),'black-+','LineWidth',1,'markersize',2,'MarkerIndices',1:0.1*k_f:k_f);
    
    min_y = min(ylim); max_y = max(ylim);
    if Missed_on_off==0
        for i_mea=1:length(RelMea_Table)
            if(RelMea_Table{i_mea}(1,1) > 0) % if it is a reasonable measurement
                mea_start = RelMea_Table{i_mea}(1,1)*delta;
                mea_end = RelMea_Table{i_mea}(1,2)*delta;
                line([floor(mea_start),floor(mea_start)],[min_y,max_y],'color','green','linestyle','--')
                line([floor(mea_end),floor(mea_end)],[min_y,max_y],'color','blue','linestyle','--')
            end
        end
%         axis([0 250 -10 10])
        xlim([0 t_f]);
    else
%         axis([0 250 -30 30])
    end
    xlabel('Time/s');
    ylabel('Error/m');
end

%% plot trace of the aggregated covariance matrix
figure; box on;
hold on;
% plot(Tk,Tr_DR,'linewidth',1,'color','black');
% plot(Tk,Tr_CEKF,'linewidth',1);
plot(Tk,Tr_TEKF,'linewidth',1);
plot(Tk,Tr_SATEKF,'linewidth',1);
[V_1,~] = TrajGen(1);
% bound_CEKF = tr_bound(P_CEKF{ceil(mea_start/delta)+1},V_1(1),delta);
% line([mea_start mea_end],[bound_CEKF bound_CEKF],'linewidth',2,'linestyle','--');
% line([mea_start mea_end],[trace(Pi_bar_opt) trace(Pi_bar_opt)],'linewidth',2,'linestyle','--','color','blue');
% line([mea_start mea_end],[trace(Pi_bar_init) trace(Pi_bar_init)],'linewidth',2,'linestyle','-.','color','red');
% plot(Tr_ICI,'linewidth',1);
% plot(Tr_ICI_EKF,'linewidth',1);
% xlim([max(Tk)*49/60 max(Tk)*55/60]);
% ylim([0 trace(Pi_bar_opt)*1.2]);
xlabel('Time/s');
ylabel('Trace of covariance/m^2');
% ylim([0 0.02]);
% legend('DR','CEKF','SA-TEKF','TEKF','Upper bound-OPT','Upper bound-Init','location','best');
% legend('DR','CEKF');
% legend('DR','Dense-CL','OPT-CL');
legend('Dense-CL','OPT-CL');
%% Communication times evolution
% figure;
% subplot(2,1,1);
% hold on;
% box on;
% plot(Tk,count_com_TEKF(:,1),'linewidth',1);
% plot(Tk,count_com_SATEKF(:,1),'linewidth',1);
% plot(Tk,count_com_TEKF(:,2),'linewidth',1);
% ylabel('# Communications');
% legend('TEKF','SA-TEKF','CEKF','location','best');
% subplot(2,1,2);
% hold on; box on;
% plot(Tk,count_com_TEKF(:,1)./count_com_TEKF(:,2),'linewidth',1);
% ylabel('Ratio');
% xlabel('Time/s');
%% Topology changing
G = digraph();
for i=1:N
    for j=1:N
        if i~=j
            G = addedge(G,i,j);
        end
    end
end
figure;
plot(G);

G = cell(1,5);
for i=1:5
    G{i} = digraph();
    G{i} = addnode(G{i},N);
end
cur_edge = zeros(3,2); ii = 1;
time_edge = zeros(5,1);
for i=1:length(mea_schedule)
    if ~isempty(mea_schedule{i})
        if ~isequal(mea_schedule{i}(:,2:3),cur_edge)
            cur_edge = mea_schedule{i}(:,2:3);
            G{ii} = digraph(table(cur_edge,'VariableNames',{'EndNodes'}));
            time_edge(ii) = floor(i*delta);
%             if G{ii}.numnodes~=10
%                 disp('bad!');
%             end
            ii = ii + 1;
        end
    end
end
figure;
for i=1:4
    j = 3*i;
    subplot(2,2,i);
    plot(G{j});
    xlabel(['t = ',num2str(time_edge(j)),'s, # Edges = ',num2str(G{j}.numedges)]);
end