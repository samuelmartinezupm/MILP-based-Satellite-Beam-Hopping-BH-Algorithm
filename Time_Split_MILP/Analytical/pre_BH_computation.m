 function [UpC,Adj_c,Adj_u,D,P,R,c_scenario,theta,M,N]=pre_BH_computation(scenario,freq,h_sat,el_min,Re,frame,frame_dur,colours,beams,P_T,n_users,traffic_model,B_T,N,rings,number_cells); % Precompute scenery based matrixes.

global PWD;

% Helping Matrixes:
UpC=zeros(number_cells,n_users);
Adj_c=zeros(number_cells,number_cells);
Adj_u=zeros(n_users,n_users);

SE=[0.1523, 0.3770, 0.8770, 1.4766, 1.9141, 2.4063, 2.7305, 3.3223, 3.9023, 4.5234, 5.1152, 5.5547, 6.2266, 6.9141, 7.4063]; % MCS index table 2 for PDSCH
M=length(SE);

for CQI=1:length(SE)
    C_N(CQI)=(CQI-2.9142)/0.2789;
    if C_N(CQI)>=4
        C_N(CQI)=(CQI-1.8749)/0.5246;
        if C_N(CQI)>=14
            C_N(CQI)=(CQI-1.4946)/0.5369;
        end
    end
end

%  PAINT
% C/N vs SE plot:
% figure
% plot(C_N, SE,'bx')
% hold on
% stairs(C_N, SE)
% xlabel('C/N [dB]')
% ylabel('SE [bit/s/Hz]')
% title('MCS PDSCH 5G')
% hold off


% 1) Ring-based cell generation
% Trigonometric raltionships:
betta_max=(asin((Re/(Re+h_sat))*cos(el_min*pi/180)))*180/pi;
gamma_max=90-betta_max-el_min;
h_prime_max=Re*(1-cos(gamma_max*pi/180));
A_footprint=2*pi*Re*h_prime_max;
D_footprint=((2*gamma_max)/360)*(2*pi*Re);

% Cell Definition:
[l,h,theta,D_footprint,number_cells,centers, interfering,Adj_c]=Cell_Scenario(h_sat,el_min,rings,Adj_c);


%% PAINT
% % Satellite Footprint:
% figure (2)
% pgon = nsidedpoly(1000, 'Center',  [0,0], 'Radius', D_footprint/2);
% plot(pgon, 'FaceColor', 'y');
% axis equal
% hold on

% Cell Footprint 
c_scenario=[];
for i=1:number_cells
    c_scenario=[c_scenario c(i,centers(i,:),l,nonzeros(interfering(i,:)))]; %constructor
    c_scenario(i).compute_betta_to_sat(h_sat); % Compute betta to the cell center so that the gain loss due to beam scanning (moving ftom boresight) can be then accounted: non-ideal isotropic behavior of the embedded element gain.
   % % PAINT
   % c_scenario(i).draw(1);
   % hold on
end
% % PAINT
% title(strcat('Cell Scenario:  ',num2str(rings),' rings'))

%2) User and traffic generation:
%User Generation

[x_y,demand,type,g_rx,T_noise_rx]=Traffic_Distribution(traffic_model,n_users,D_footprint, freq);

%User Aggregation into cells:
for i=1:n_users
    %Euclidean distance to cell centers to determine the cell:
    [dist_min,cell_num]=min(sqrt((x_y(i,1)-centers(:,1)).^2+(x_y(i,2)-centers(:,2)).^2)); 
    %Assign user to corresponding cell:
    c_scenario(cell_num).adduser(u(i,[x_y(i,1),x_y(i,2)],type,demand(i),g_rx,T_noise_rx));  %constructor
    c_scenario(cell_num).users(length(c_scenario(cell_num).users)).compute_distance_elevation_betta_to_sat(h_sat); % Compute distance, elevation and betta angle to the satellite of the the last added user
    c_scenario(cell_num).users(length(c_scenario(cell_num).users)).compute_betta_to_cell_center(h_sat,c_scenario(cell_num)); % Compute the betta angle with respect to the cell center so that the gain loss can be accounted by the fact of not being at the cell center
    UpC(cell_num,i)=1;
end
Adj_u=UpC'*Adj_c*UpC;


%3) 

% Scenery based MATRIX generation initialization:
P=zeros(M,N,n_users);
D=zeros(M,N);
R=zeros(n_users,frame);

% Link Budget Caculation for each M (MODCOD and therefore C/N->SE), N (numer of assigned BAND) and u (user, depending on the scenario):
% The Demand matrix (D) does not change, depending on the user. It is the
% power the one varies from user to user depending on the particular link
% conditions. 

% Input Kink Budget data:
% CANCELAR EXTRA LOSSES:
Extra_losses.At=0;
G_t=10*log10(0.65*48360/theta^2); % It needs to be corrected by user's position.
k=-228.6012; %dBW/(HzK^-1)
T_ant=30; %K
alpha=0.1;

% Antenna pattern
q=1.3;
np=1001;
theta_scanning=linspace(-pi/2,pi/2,np);

for i=1:M
    for j=1:N
        D(i,j)=SE(i)*((j*(B_T/N))/(1+alpha))*frame_dur; % [Mbits/s]*frame_dur
        for k_c=1:number_cells
            for k_u=1:length(c_scenario(k_c).users)
                %R(c_scenario(k_c).users(k_u).id,:)=c_scenario(k_c).users(k_u).traffic_demand*frame_dur; % [Mbits/s]*frame_dur
                % C_N=(G_tx+P_tx-Lfsl-Extra_losses.At+G_rx)-(k+T_eq+BW_rx)->Determine:P_tx
                 % Gain TX:
                 G_tx=G_t-12*(c_scenario(k_c).users(k_u).betta_to_center/theta)^2; % Gain loss approximation when moving away from the cell center.
                 [minScanning,closestScanningIndex]=min(abs(theta_scanning(1:((np-1)/2+1))+deg2rad(c_scenario(k_c).betta_to_sat))); % Gain loss due to scanning
                 G_tx=G_tx+20*log10(cos(abs(theta_scanning(closestScanningIndex)))^q);
                 % Gain RX:
                 G_rx=c_scenario(k_c).users(k_u).gain; %dB           
                 % T_eq RX:
                 T_eq=10*log10(T_ant+c_scenario(k_c).users(k_u).T_noise+280*(1-10^(-Extra_losses.At/10))); %[dBK] -> Contributors: captured by antenna + rx equivalent + attenuation due to extra losses
                 Lfsl=20*log10(freq)+20*log10(c_scenario(k_c).users(k_u).distace_to_sat*10^3)-147.55;
                 BW_rx=10*log10(j*(B_T/N)*10^6); 
                 P(i,j,c_scenario(k_c).users(k_u).id)=10^((C_N(i)-G_tx+Lfsl+Extra_losses.At-G_rx+(k+T_eq+BW_rx))/10); % [W]
            end
        end
    end
end

%save('input_BH_data','UpC','Adj_c','Adj_u','D','P');
% 
for t=1:frame %1frame=10 subframe by 10subframe (beamswithching is performed each frame, for 10 slots constant illumination constant)
    % Traffic Pending Adjustment based on frame length: from second to second add the total demand/s to the pending:
    if mod(t-1,1/frame_dur)==0
        for i=1:length(c_scenario)
            for j=1:length(c_scenario(i).users)
                R(c_scenario(i).users(j).id,t)=c_scenario(i).users(j).traffic_pending+c_scenario(i).users(j).traffic_demand;
            end
        end
    end
end


end
