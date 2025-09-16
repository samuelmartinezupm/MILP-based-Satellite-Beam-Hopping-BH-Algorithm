function [UC_perc,EC_perc,TTS,UC_perc_db,EC_perc_db,TTS_db, o1, o2, o3, Ill_out, B_out, P_out, MS, GS, FS, CS, ZS, XS, PS, DS]=FOM_Estimation(c_scenario,B_T,P_T,Adj_c,Adj_u,UpC,n_users,number_cells,rings,beams,theta,colours,frame,frame_dur,freq, M,N,P,D,solutions,i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range, t_current)

% It's importat to consider just till the current framing (frame_current), not to alter the FoM:
frame_current=t_current+9;

%% 1) ANALYTICAL SOLUTION
result=solutions;
% Reshape result back: [Ill], [B] and [P].
Ill_out=reshape(result(1:i_range),[number_cells,frame]);
B_out=reshape(result(b_index:(b_index+b_range-1)),[n_users,N,frame]);
MS=reshape(result(m_index:(m_index+m_range-1)),[n_users,M,frame]);
GS=reshape(result(g_index:(g_index+g_range-1)),[n_users,frame]);
FS=reshape(result(f_index:(f_index+f_range-1)),[n_users,frame]);
CS=reshape(result(c_index:(c_index+c_range-1)),[number_cells,frame]);
ZS=reshape(result(z_index:(z_index+z_range-1)),[n_users,M,N,frame]);
XS=reshape(result(x_index:(x_index+x_range-1)),[n_users,N,frame]);

P_out=zeros(n_users,frame);
for u=1:n_users
    for t=1:frame
        [m,n]=find(round(squeeze(ZS(u,:,:,t))));
        if isempty(m) || isempty(n)
             P_out(u,t)=0;
        else
            P_out(u,t)=P(m,n,u);
        end
    end
end

PS=zeros(n_users,frame);
for u=1:n_users    
    for t=1:frame
        for m=1:M
            for n=1:N
                if ZS(u,m,n,t)*P(m,n,u)>0
                    PS(u,t)=PS(u,t)+ZS(u,m,n,t)*P(m,n,u); % iT SHOULD ENTER HERE ONLY FOR A SINGLE PAIR M,N
                end
            end
        end
    end
end

DS=zeros(n_users,frame);
for u=1:n_users    
    for t=1:frame
        for m=1:M
            for n=1:N
                if ZS(u,m,n,t)*D(m,n)>0
                    DS(u,t)=DS(u,t)+ZS(u,m,n,t)*D(m,n); % iT SHOULD ENTER HERE ONLY FOR A SINGLE PAIR M,N
                end
            end
        end
    end
end

b_slots=N;

% Matrix cleaning till frame_current:
Ill_out=Ill_out(:,1:frame_current);
B_out=B_out(:,:,1:frame_current);
MS=MS(:,:,1:frame_current);
GS=GS(:,1:frame_current);
FS=FS(:,1:frame_current);
CS=CS(:,1:frame_current);
ZS=ZS(:,:,:,1:frame_current);
XS=XS(:,:,1:frame_current);
PS=PS(:,1:frame_current);
DS=DS(:,1:frame_current);
P_out=P_out(:,1:frame_current);

%save(strcat('[',num2str(params.MIPGap),'MIPGap',num2str(weight_traffic_UC),'_',num2str(weight_traffic_EC),'_',num2str(weight_timing),'_',num2str(n_users),'u_', num2str(frame),'frames_', num2str(b_slots),'bslots_normalized]_GUROBI_c_scenario_',num2str(scenario),'_analytical','.mat'), "c_scenario","B_T","b_slots", "P_T", "Adj_c", "Adj_u", "UpC", "n_users", "rings","beams", "theta", "colours", "frame", "frame_dur", "freq","Ill_out","B_out","P_out");

%objetivo=(sum(GS,"all")*weight_traffic_UC+sum(FS,"all")*weight_traffic_EC)+sum(CS,"all")*weight_timing;

o1=sum(GS,"all");
o2=sum(FS,"all");
o3=sum(CS,"all");

%% 2) DEMAND BASED:
    % NO Band Segregation:
    %[RC,SC,UC,EC,TTS, Ill,B,P]=FOM_calculation_demand_based(B_T,b_slots, P_T, Ill, B, P,c_scenario,n_users, beams, theta, colours, frame, frame_dur, TTL, freq);
    %save(strcat('[demand_based_NO_band_segregation]data_',string(beams),'beams_',string(colours),'colours'),"SC","UC","EC","RC","TTS","Ill", "B", "P", "Adj_c","Adj_u","UpC")
    
    %Tentative Outputs:
    Ill=zeros(number_cells,frame_current);
    B=zeros(n_users,b_slots,frame_current);
    P=zeros(n_users,frame_current);

    TTL=5;

    %Band Segregation:
    [RC_db,SC_db,UC_db,EC_db,TTS_db, Ill_db,B_db,P_db]=FOM_calculation_demand_based_band_slots(B_T,b_slots, P_T, Ill, B, P,c_scenario,n_users, beams, theta, colours, frame_current, frame_dur, TTL, freq);
    % save(strcat('[EXAMPLE]_demand_based_FOM_and_Ill_B_P'),"SC","UC","EC","RC","TTS","Ill", "B", "P")

%% 3) OPTIMIZATION-BASED:
    % load("[demand_based_band_segregation]data_4beams_1colours.mat","Ill","B","P") % Let's start with demand based results!
    % load("[demand_based_NO_band_segregation]data_4beams_1colours.mat","Adj_c","Adj_u","UpC") 
    [RC,SC,UC,EC,TTS]=FOM_calculation_optimization_based(B_T,b_slots, P_T,  Ill_out, B_out, P_out, UpC, c_scenario, n_users, theta, frame_current, frame_dur, freq);
    % save(strcat('[optimization_based]data_',string(beams),'beams_',string(colours),'colours'),"SC","UC","EC","RC","TTS","Ill", "B", "P")

 UC_perc=(UC/RC)*100;
 EC_perc=(EC/RC)*100;
 TTS;


 UC_perc_db=(UC_db/RC_db)*100;
 EC_perc_db=(EC_db/RC_db)*100;
 TTS_db;

%PAINT
% plot(sum(P_out,1))
% hold on
% plot(sum(P,1))
% legend('analytical','demand-based')
% xlabel('Time instants (0.1s)')
% ylabel('P [W]')
% title('Power variation')
% hold off