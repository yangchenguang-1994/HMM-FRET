function [FRET_trace,seq_show] = HMM(photo_intensity,photo_time,sta_num)
    figure
    plot(photo_time,photo_intensity,'*','Color',[0 1 0],'MarkerSize',5);
    hold on;
    plot(photo_time,photo_intensity,'-','Color',[0 1 0],'LineWidth',1.4);
    hold on;
    zoom on;
    pause;
    zoom off;
    abc=ginput(sta_num);
%     top = abc(1,2);
%     bottom = abc(2,2);
%     T1 = length(photo_time);
    % HMM fit
    FRET_trace=cell(1,1);
    FRET_trace{1}=[photo_time photo_intensity];
    trace_num=1; 
    [timetags,sd,cent] = kmean(photo_intensity,sta_num,1000);
    cent= abc(1:end,2);%可以设置可能的初始值来迭代
    T=length(FRET_trace{trace_num});
    %设置A，B，pi初始值,包含几个未知自由变量由sta_num决定
    %创建A，B，Pi向量
    A_matrix_0 = zeros(sta_num,sta_num)+ 1/sta_num;%设置A初始值

    B_matrix_para_0=zeros(sta_num,2);  %设置B迭代初始值
    for i=1:sta_num
        B_matrix_para_0(i,1)=cent(i);             %状态大体值
        B_matrix_para_0(i,2)=sd(i);                    %每个态的标准差
    end
    Pi_matrix_0=zeros(1,sta_num)+1/sta_num;          %设置Pi初始值
    %EM算法求极值
    max_iter=150;   %迭代次数
    A_matrix_1=A_matrix_0;
    B_matrix_para_1=B_matrix_para_0;
    Pi_matrix_1=Pi_matrix_0;
    B_calcu=zeros(sta_num,T); %根据中状态序列计算B_matrix所用到的工具矩阵，每次迭代要初始化
    B_calcu_1=zeros(sta_num,1);
    P_lase=zeros(2,max_iter);
    %EM算法E步
    for i=1:max_iter
      %由观测序列trace和初始参数得到状态序列（维特比算法）
      [seq,P]=viterbi_FRET(A_matrix_0,B_matrix_para_0,Pi_matrix_0,FRET_trace{trace_num}(:,2)');
      P_lase(i)=P;
      %由挂测序列trace和状态序列得到参数
      %求B_matrix
      B_calcu=zeros(sta_num,T); %根据中状态序列计算B_matrix所用到的工具矩阵，每次迭代要初始化
      B_calcu_1=zeros(sta_num,1);
      for j=1:T
          state_j=seq(j);
          B_calcu_1(state_j)=B_calcu_1(state_j)+1;
          B_calcu(state_j,B_calcu_1(state_j))=FRET_trace{trace_num}(j,2)'; 
      end
      for k=1:sta_num 
          B_matrix_para_1(k,1)=mean(B_calcu(k,1:B_calcu_1(k)),2);
          B_matrix_para_1(k,2)=std(B_calcu(k,1:B_calcu_1(k)),0,2);
      end
      %B_matrix_para_0
      %求A_matrix
      A_matrix_1=zeros(sta_num,sta_num);
      for m=1:T-1
          A_matrix_1(seq(m),seq(m+1))=A_matrix_1(seq(m),seq(m+1))+1;
      end
    %   A_matrix_1
      %归一化A_matrix
      for n=1:sta_num
          A_matrix_1(n,:)=A_matrix_1(n,:)/(sum(A_matrix_1(n,:),2)+0.01);
      end

      B_matrix_para_0=B_matrix_para_1;
      A_matrix_0=A_matrix_1;

    end

    %seq_show是状态序列，不过状态值不是1 2 3，而是对应的FRET值
    P_lase(2,1:max_iter-1)=P_lase(1,2:max_iter)-P_lase(1,1:max_iter-1);
    seq_show=seq;
    for o=1:T
        seq_show(o)=B_matrix_para_0(seq(o),1);
    end
    figure
    plot(FRET_trace{trace_num}(:,1),seq_show,'LineWidth',2)
    hold on
    plot(FRET_trace{trace_num}(:,1),FRET_trace{trace_num}(:,2))
end