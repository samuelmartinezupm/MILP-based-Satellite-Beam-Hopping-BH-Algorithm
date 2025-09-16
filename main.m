clc 
clear all
close all

global PWD;

% USE CASE DEFINITION VARIABLES:
scenario=1;
run=1;

freq=20e9; % Ka Band (DL)
h_sat=1100; % Starklink [550km] ; OneWeb [1100km] 
el_min=25; % Minimum Elevation Angle [º] 
Re=6378; % Earth Radius [km]
frame=30; %100ms/frame = 10s
frame_dur=0.1; %s
colours=1;
beams=4;
P_T=18; %W
n_users=5; 
traffic_model='uniform'; %'linear' %'hotspot'
B_T=250/2; %MHz per colour
N=B_T/25; %Bandwidth bins
rings=3; %Fixed grid 7 cellsP
number_cells=0;
for i=1:rings
    number_cells=number_cells+6*(i-1);
end
number_cells=number_cells+1;
    
% SELECTRESOLUTION METHODOLOGY: FULL MILP vs TIME-SPLIT MILP
method='time-split'; % 'full';



if strcmp(method,'full') % FULL MILP 

    disp('Loading Full MILP...');

    cd ./Full_MILP
    PWD=pwd;
    cd (strcat(PWD,'/Analytical/'))
    
    [UpC,Adj_c,Adj_u,D,P,R,c_scenario,theta,M,N,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range,row,col,val,total_range,total_constraints,b, vtype,mip_start,g_aux,f_aux,c_aux]=BH(scenario,freq,h_sat,el_min,Re,frame,frame_dur,colours,beams,P_T,n_users,traffic_model,B_T,N,rings,number_cells);    
    
    for betta=[0,1,0.7]
        if betta==0 % NORMALIZATION: max_UC_betta_0
            MIPGap=0.4;
            normalization_UC=1;
            normalization_EC=1;
            normalization_time=1;
            cd (PWD)
            [result,data,solutions]=gurobi_execution(n_users,number_cells,frame,betta,normalization_UC,normalization_EC,normalization_time,total_range,total_constraints,val,row,col,b,vtype,MIPGap,mip_start,g_aux,f_aux,c_aux);
            cd (PWD)
            [max_UC_betta_0,~,~]=solution_plot_saving(scenario,run,betta,normalization_UC,normalization_EC,normalization_time,data,solutions,c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range);
    
        elseif betta==1 % NORMALIZATION: max_EC_betta_1, max_TTS_betta_1
            MIPGap=0.4;
            normalization_UC=1;
            normalization_EC=1;
            normalization_time=1;
            cd (PWD)
            [result,data,solutions]=gurobi_execution(n_users,number_cells,frame,betta,normalization_UC,normalization_EC,normalization_time,total_range,total_constraints,val,row,col,b,vtype,MIPGap,mip_start,g_aux,f_aux,c_aux);
            cd (PWD)
            [~,max_EC_betta_1, max_TTS_betta_1]=solution_plot_saving(scenario,run,betta,normalization_UC,normalization_EC,normalization_time,data,solutions,c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range);
    
        else % FULL MILP EXECUTION
            MIPGap=0.2;
            normalization_UC=max_UC_betta_0;
            if max_EC_betta_1==0
                normalization_EC=1;
            else
                normalization_EC=max_EC_betta_1;
            end
            if max_TTS_betta_1==0
                normalization_time=1;
            else
                normalization_time=max_TTS_betta_1;
            end
            cd (PWD)
            [result,data,solutions]=gurobi_execution(n_users,number_cells,frame,betta,normalization_UC,normalization_EC,normalization_time,total_range,total_constraints,val,row,col,b,vtype,MIPGap,mip_start,g_aux,f_aux,c_aux);
            cd (PWD)
            [~,~,~]=solution_plot_saving(scenario,run,betta,normalization_UC,normalization_EC,normalization_time,data,solutions,c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range);
    
        end  
    end

elseif strcmp(method,'time-split') %TIME-SPLIT MILP

   disp('Loading Time-Split MILP...');

    cd ./Time_Split_MILP
    PWD=pwd;

    cd (strcat(PWD,'/Analytical/'))
    
    [c_scenario,Adj_c,Adj_u,UpC,theta,M,N,P,D,R,A,total_range,total_constraints,b, vtype,mip_start,i_index,m_index,b_index,x_index,z_index,g_index,f_index,c_index,a_p_index,a_n_index,i_range,m_range,b_range,x_range,z_range,g_range,f_range,c_range,a_p_range,a_n_range,i_aux,m_aux,b_aux,x_aux,z_aux,g_aux,f_aux,c_aux,a_p_aux,a_n_aux,t_contraint_label]=BH(scenario,freq,h_sat,el_min,Re,frame,frame_dur,colours,beams,P_T,n_users,traffic_model,B_T,N,rings,number_cells);

    % Restrict ALL Variables:
    lb=zeros(1,total_range);
    ub=zeros(1,total_range);
    t_range_boolean=zeros(1,total_range); % boolean variable that identifies the variables that must be accounted (1), always the ones lower than t_idx+9.
        
    for t_idx=1:10:frame
        % Unrestric the variables within the time-frame:
        %lb([reshape(i_aux(:,t:t+9),1,[]), reshape(m_aux(:,:,t:t+9),1,[]), reshape(b_aux(:,:,t:t+9),1,[]), reshape(x_aux(:,:,t:t+9),1,[]), reshape(z_aux(:,:,:,t:t+9),1,[]), reshape(a_p_aux(:,:,t:t+9),1,[]), reshape(a_n_aux(:,:,t:t+9),1,[])])=0;
        ub([reshape(i_aux(:,t_idx:t_idx+9),1,[]), reshape(m_aux(:,:,t_idx:t_idx+9),1,[]), reshape(b_aux(:,:,t_idx:t_idx+9),1,[]), reshape(x_aux(:,:,t_idx:t_idx+9),1,[]), reshape(z_aux(:,:,:,t_idx:t_idx+9),1,[]), reshape(a_p_aux(:,:,t_idx:t_idx+9),1,[]), reshape(a_n_aux(:,:,t_idx:t_idx+9),1,[])])=1;
        t_range_boolean([reshape(i_aux(:,t_idx:t_idx+9),1,[]), reshape(m_aux(:,:,t_idx:t_idx+9),1,[]), reshape(b_aux(:,:,t_idx:t_idx+9),1,[]), reshape(x_aux(:,:,t_idx:t_idx+9),1,[]), reshape(z_aux(:,:,:,t_idx:t_idx+9),1,[]), reshape(a_p_aux(:,:,t_idx:t_idx+9),1,[]), reshape(a_n_aux(:,:,t_idx:t_idx+9),1,[])])=1;
        %lb([reshape(g_aux(:,t:t+9),1,[]), reshape(f_aux(:,t:t+9),1,[])])=0;
        ub([reshape(g_aux(:,t_idx:t_idx+9),1,[]), reshape(f_aux(:,t_idx:t_idx+9),1,[])])=max(sum(R(:,1:t_idx+9),2))+1;%2*max(sum(R(:,1:t_idx+9),2));%max(sum(R(:,1:t_idx+9),2))+1;
        t_range_boolean([reshape(g_aux(:,t_idx:t_idx+9),1,[]), reshape(f_aux(:,t_idx:t_idx+9),1,[])])=1;
        %lb(reshape(c_aux(:,t:t+9),1,[]))=0;
        ub(reshape(c_aux(:,t_idx:t_idx+9),1,[]))=(t_idx+9)+1;%2*(t_idx+9);%(t_idx+9)+1;
        t_range_boolean(reshape(c_aux(:,t_idx:t_idx+9),1,[]))=1;%(t_idx+9)+1;

        for betta=[0,1,0.7]
            if betta==0
                MIPGap=0.40;
                normalization_UC=1;
                normalization_EC=1;
                normalization_time=1;
                cd (PWD)
                [result,data,solutions]=gurobi_execution(n_users,number_cells,frame,betta,normalization_UC,normalization_EC,normalization_time,total_range,total_constraints,A,b,vtype,MIPGap,mip_start,g_aux,f_aux,c_aux,lb,ub,t_idx,t_contraint_label,t_range_boolean);
                cd (PWD)
                [max_UC_betta_0,~,~]=solution_plot_saving(scenario,run,betta,normalization_UC,normalization_EC,normalization_time,data,solutions,c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range,t_idx);
        
            elseif betta==1
                MIPGap=0.40;
                normalization_UC=1;
                normalization_EC=1;
                normalization_time=1;
                cd (PWD)
                [result,data,solutions]=gurobi_execution(n_users,number_cells,frame,betta,normalization_UC,normalization_EC,normalization_time,total_range,total_constraints,A,b,vtype,MIPGap,mip_start,g_aux,f_aux,c_aux,lb,ub,t_idx,t_contraint_label,t_range_boolean);
                cd (PWD)
                [~,max_EC_betta_1, max_TTS_betta_1]=solution_plot_saving(scenario,run,betta,normalization_UC,normalization_EC,normalization_time,data,solutions,c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range,t_idx);
        
            else
                MIPGap=0.2;
                normalization_UC=max_UC_betta_0;
                if max_EC_betta_1==0
                    normalization_EC=1;
                else
                    normalization_EC=max_EC_betta_1;
                end
                if max_TTS_betta_1==0
                    normalization_time=1;
                else
                    normalization_time=max_TTS_betta_1;
                end
                cd (PWD)
                [result,data,solutions]=gurobi_execution(n_users,number_cells,frame,betta,normalization_UC,normalization_EC,normalization_time,total_range,total_constraints,A,b,vtype,MIPGap,mip_start,g_aux,f_aux,c_aux,lb,ub,t_idx,t_contraint_label,t_range_boolean);
                cd (PWD)
                [~,~,~]=solution_plot_saving(scenario,run,betta,normalization_UC,normalization_EC,normalization_time,data,solutions,c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range,t_idx);
                % Fix computed result for the given time frame:
                lb([reshape(i_aux(:,t_idx:t_idx+9),1,[]), reshape(m_aux(:,:,t_idx:t_idx+9),1,[]), reshape(b_aux(:,:,t_idx:t_idx+9),1,[]), reshape(x_aux(:,:,t_idx:t_idx+9),1,[]), reshape(z_aux(:,:,:,t_idx:t_idx+9),1,[]), reshape(g_aux(:,t_idx:t_idx+9),1,[]), reshape(f_aux(:,t_idx:t_idx+9),1,[]), reshape(c_aux(:,t_idx:t_idx+9),1,[]), reshape(a_p_aux(:,:,t_idx:t_idx+9),1,[]), reshape(a_n_aux(:,:,t_idx:t_idx+9),1,[])])=result([reshape(i_aux(:,t_idx:t_idx+9),1,[]), reshape(m_aux(:,:,t_idx:t_idx+9),1,[]), reshape(b_aux(:,:,t_idx:t_idx+9),1,[]), reshape(x_aux(:,:,t_idx:t_idx+9),1,[]), reshape(z_aux(:,:,:,t_idx:t_idx+9),1,[]), reshape(g_aux(:,t_idx:t_idx+9),1,[]), reshape(f_aux(:,t_idx:t_idx+9),1,[]), reshape(c_aux(:,t_idx:t_idx+9),1,[]), reshape(a_p_aux(:,:,t_idx:t_idx+9),1,[]), reshape(a_n_aux(:,:,t_idx:t_idx+9),1,[])]);
                ub([reshape(i_aux(:,t_idx:t_idx+9),1,[]), reshape(m_aux(:,:,t_idx:t_idx+9),1,[]), reshape(b_aux(:,:,t_idx:t_idx+9),1,[]), reshape(x_aux(:,:,t_idx:t_idx+9),1,[]), reshape(z_aux(:,:,:,t_idx:t_idx+9),1,[]), reshape(g_aux(:,t_idx:t_idx+9),1,[]), reshape(f_aux(:,t_idx:t_idx+9),1,[]), reshape(c_aux(:,t_idx:t_idx+9),1,[]), reshape(a_p_aux(:,:,t_idx:t_idx+9),1,[]), reshape(a_n_aux(:,:,t_idx:t_idx+9),1,[])])=result([reshape(i_aux(:,t_idx:t_idx+9),1,[]), reshape(m_aux(:,:,t_idx:t_idx+9),1,[]), reshape(b_aux(:,:,t_idx:t_idx+9),1,[]), reshape(x_aux(:,:,t_idx:t_idx+9),1,[]), reshape(z_aux(:,:,:,t_idx:t_idx+9),1,[]), reshape(g_aux(:,t_idx:t_idx+9),1,[]), reshape(f_aux(:,t_idx:t_idx+9),1,[]), reshape(c_aux(:,t_idx:t_idx+9),1,[]), reshape(a_p_aux(:,:,t_idx:t_idx+9),1,[]), reshape(a_n_aux(:,:,t_idx:t_idx+9),1,[])]);
                
                % Adjust and recompute DB solution:
                if (t_idx+9)<frame
                    TTL=5;
                    [~,~,~,~,~, Ill_db,B_db,P_db]=FOM_calculation_demand_based_band_slots_fixed_MODCODs_SPLITTED(B_T,N,M, P_T, zeros(number_cells,frame), zeros(n_users,N,frame), zeros(n_users,frame), c_scenario, UpC, n_users, number_cells, beams, theta, colours, frame, frame_dur, TTL, freq, P, D, result, i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range,t_idx);
                    i_start=Ill_db(1:number_cells*frame);
                    b_start=B_db(1:n_users*N*frame);
                    
                    X_db=zeros(n_users,N,frame);
                    M_db=zeros(n_users,M,frame);
                    Z_db=zeros(n_users,M,N,frame);
                    for time_idx=1:frame
                        for u_idx=1:n_users
                           if P_db(u_idx,time_idx)>0 && sum(B_db(u_idx,:,time_idx))>0
                                X_db(u_idx,uint8(sum(B_db(u_idx,:,time_idx))),time_idx)=1;
                    
                                p_idx=find(P_db(u_idx,time_idx)==P(:,uint8(sum(B_db(u_idx,:,time_idx))),u_idx));
                                if p_idx==0
                                disp()
                                end
                                M_db(u_idx,p_idx,time_idx)=1;
                                Z_db(u_idx,p_idx,uint8(sum(B_db(u_idx,:,time_idx))),time_idx)=1;
                            end
                        end
                    end
                    
                    x_start=X_db(1:n_users*N*frame);
                    m_start=M_db(1:n_users*M*frame);
                    z_start=Z_db(1:n_users*M*N*frame);
        
                    mip_start=[i_start,m_start,b_start,x_start,z_start,NaN(1,total_range-length([i_start,m_start,b_start,x_start,z_start]))];
                end
    
            end
        
        end
    end


else 
    disp('Wrong MILP method selection.')
end