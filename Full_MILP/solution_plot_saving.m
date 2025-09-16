%% SOLUTION PLOTTING:
function [max_UC_betta_0,max_EC_betta_1,max_TTS_betta_1]=solution_plot_saving(scenario,run,betta,normalization_UC,normalization_EC,normalization_time,data,solutions,c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range)

    global PWD;
    
    % Weighting parameter estimation:
    weight_traffic_EC=(1-betta)/2;
    weight_traffic_UC=betta;
    weight_timing=(1-betta)/2;
    
    cd (strcat(PWD,'/FoM_Estimation/'))
    OBJETIVO=zeros(1,size(solutions,1));
    O1=zeros(1,size(solutions,1));
    O2=zeros(1,size(solutions,1));
    O3=zeros(1,size(solutions,1));
    UC_perc_M=zeros(1,size(solutions,1));
    EC_perc_M=zeros(1,size(solutions,1));
    TTS_M=zeros(1,size(solutions,1));
    
    for i=1:size(solutions,1)
        [UC_perc,EC_perc,TTS,UC_perc_db,EC_perc_db,TTS_db, o1, o2, o3, Ill_out, B_out, P_out, MS, GS, FS, CS, ZS, XS, PS, DS]=FOM_Estimation(c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,solutions(i,:),i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range);
        O1(i)=o1;
        O2(i)=o2;
        O3(i)=o3;
        OBJETIVO(i)=(weight_traffic_UC.*O1(i))./normalization_UC + (weight_traffic_EC.*O2(i))./normalization_EC +(weight_timing.*O3(i))./normalization_time;
        UC_perc_M(i)=UC_perc;
        EC_perc_M(i)=EC_perc;
        TTS_M(i)=TTS;
    end
    
    Filename=strcat('[',num2str(scenario+1),'s_',num2str(run),'r]_',num2str(betta),'betta_',num2str(frame),'frames_',num2str(n_users),'u','.mat');
    b_slots=N;
    save([fullfile(fileparts(PWD), 'Simulation_Results'),'/',Filename],"UC_perc_M","EC_perc_M","TTS_M","UC_perc_db","EC_perc_db","TTS_db","O1","O2","O3","OBJETIVO","normalization_UC","normalization_EC","normalization_time","c_scenario","B_T","b_slots","P_T","Adj_c","Adj_u","UpC","n_users","rings","beams","theta","colours","frame","frame_dur","freq","Ill_out","B_out","P_out","MS","GS","FS","CS","ZS","XS","PS","DS","data","solutions","betta")
    
    
    % BETTA DATA:
    max_UC_betta_0=sum(GS,"all"); %o1
    max_EC_betta_1=sum(FS,"all"); %o2
    max_TTS_betta_1=sum(CS,"all"); %o3
    
    % % PLOT:
    % figure
    % hold on
    % plot(data(4:end,1),data(4:end,2))
    % plot(data(4:end,1),data(4:end,3))
    % plot(data(4:end,1),data(end,2)*ones(1,length(data(4:end,1))),'black --') % Objective Value
    % legend('Incumbent Bound','Best Bound')
    % title(strcat('Incumbent and Best Bound approximation (',num2str(n_users),'u) (',num2str(frame),'frames) (betta=',num2str(betta),')'))
    % hold off
    % 
    % figure
    % hold on
    % %yyaxis left
    % plot(data(2:end,1),O1(2:end)./normalization_UC,"red -")
    % plot(data(2:end,1),O2(2:end)./normalization_EC,"blue -")
    % plot(data(2:end,1),O3(2:end)./normalization_time,"magenta -")
    % %plot(data(2:end,1),((weight_traffic_UC.*O1(2:end))./normalization_UC + (weight_traffic_EC.*O2(2:end))./normalization_EC +(weight_timing.*O3(2:end))./normalization_time),"black -")
    % plot(data(2:end,1),OBJETIVO(2:end),"black -")
    % ylabel ('FOM_{normalized}')
    % %yyaxis right
    % %ylabel ('FOM_{normalized} O3')
    % hold off
    % title(strcat('Convergence evaluation (',num2str(n_users),'u) (',num2str(frame),'frames) (betta=',num2str(betta),')'))
    % xlabel('Analytical simulation time (s)')
    % legend({'O1','O2','O3','O'})
    % 
    % 
    % figure
    % hold on
    % yyaxis left
    % plot(data(2:end,1),(UC_perc_M(2:end)),"red -")
    % plot(data(2:end,1),(UC_perc_db)*ones(1,length(data(2:end,1))),"red--")
    % plot(data(2:end,1),(EC_perc_M(2:end)),"blue -")
    % plot(data(2:end,1),(EC_perc_db)*ones(1,length(data(2:end,1))),"blue--")
    % ylabel ('FOM (%)')
    % yyaxis right
    % plot(data(2:end,1),(TTS_M(2:end)),"magenta -")
    % plot(data(2:end,1),(TTS_db)*ones(1,length(data(2:end,1))),"magenta--")
    % hold off
    % title(strcat('Convergence evaluation (',num2str(n_users),'u) (',num2str(frame),'frames) (betta=',num2str(betta),')'))
    % xlabel('Analytical simulation time (s)')
    % ylabel ('FOM (s)')
    % legend({'UC','UC -- DB','EC','EC -- DB','TTS','TTS -- DB'})

end

