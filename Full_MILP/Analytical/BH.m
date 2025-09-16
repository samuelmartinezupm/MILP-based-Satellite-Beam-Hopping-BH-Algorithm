function [UpC,Adj_c,Adj_u,D,P,R,c_scenario,theta,M,N, i_range, b_index, b_range, m_index, m_range, g_index, g_range, f_index, f_range, c_index, c_range, z_index, z_range, x_index, x_range,row,col,val,total_range,total_constraints,b, vtype,mip_start,g_aux,f_aux,c_aux]=BH(scenario,freq,h_sat,el_min,Re,frame,frame_dur,colours,beams,P_T,n_users,traffic_model,B_T,N,rings,number_cells)

global PWD;

seeds = [2086611235, 134546460, 372556643, 174984732, 1701384495, 1849275000, 319285460, 814873822, 1149458473, 1037906549, 1936202636, 752251260, 331834558, 1073404773, 473651765, 495752086, 1890790170, 30767557, 149610889, 549432576, 1899222341, 1305877330, 965789845, 136830328, 692993820, 173302510, 1936751275, 435756074, 1429130986, 1202487946, 1330904806, 1561538049, 1584558490, 1773010390, 475077922, 704819280, 1101917058, 916813798, 896500662, 1209693945, 1103696678, 1351622864, 1395343572, 978851367, 739895583, 295136723, 1425639839, 770762499, 119047332, 618439756];
rng(seeds(scenario+1));
[UpC,Adj_c,Adj_u,D,P,R,c_scenario,theta,M,N]=pre_BH_computation(scenario,freq,h_sat,el_min,Re,frame,frame_dur,colours,beams,P_T,n_users,traffic_model,B_T,N,rings,number_cells); % Precompute scenery based matrixes.
rng('shuffle');

% MIP Start:
cd (strcat(PWD,'/FoM_Estimation/'))
TTL=5;

[~,~,~,~,~, Ill_db,B_db,P_db]=FOM_calculation_demand_based_band_slots_fixed_MODCODs(B_T,N, P_T, zeros(number_cells,frame), zeros(n_users,N,frame), zeros(n_users,frame),c_scenario,n_users, beams, theta, colours, frame, frame_dur, TTL, freq,P);

%[RC,SC,UC,EC,TTS, Ill,B,P]=FOM_calculation_demand_based_band_slots_fixed_MODCODs(B_T,N, P_T, zeros(number_cells,frame), zeros(n_users,N,frame), zeros(n_users,frame),c_scenario,n_users, beams, theta, colours, frame, frame_dur, TTL, freq,P);

%[RC,SC,UC,EC,TTS, Ill,B,P]=FOM_calculation_demand_based_band_slots(B_T,25, P_T, zeros(number_cells,frame), zeros(n_users,N,frame), zeros(n_users,frame),c_scenario,n_users, 8, theta, colours, frame, frame_dur, TTL, freq);

i_start=Ill_db(1:number_cells*frame);
b_start=B_db(1:n_users*N*frame);

X_db=zeros(n_users,N,frame);
M_db=zeros(n_users,M,frame);
Z_db=zeros(n_users,M,N,frame);
for t_idx=1:frame
    for u_idx=1:n_users
       if P_db(u_idx,t_idx)>0 && sum(B_db(u_idx,:,t_idx))>0
            X_db(u_idx,sum(B_db(u_idx,:,t_idx)),t_idx)=1;

            p_idx=find(P_db(u_idx,t_idx)==P(:,sum(B_db(u_idx,:,t_idx)),u_idx));
            if p_idx==0
            disp()
            end
            M_db(u_idx,p_idx,t_idx)=1;
            Z_db(u_idx,p_idx,sum(B_db(u_idx,:,t_idx)),t_idx)=1;
        end
    end
end

x_start=X_db(1:n_users*N*frame);
m_start=M_db(1:n_users*M*frame);
z_start=Z_db(1:n_users*M*N*frame);

cd (strcat(PWD,'/Analytical/'))

% Number of interference users: 
interference=0;
for u_fila=1:size(Adj_u,1)
    for u_columna=1:size(Adj_u,2)
        if Adj_u(u_fila,u_columna)==1 && u_fila<u_columna
            interference=interference+1;
        end
    end
end

% Set of variables: x=[i,m,b,x,z,g,c]
i_index=1;
i_range=number_cells*frame;
i_aux=reshape(1:i_range,[number_cells,frame]);
vtype(1:i_range)='B';

m_index=i_index+i_range;
m_range=n_users*M*frame;
m_aux=reshape(m_index:(m_index+m_range-1),[n_users,M,frame]);
vtype(m_index:(m_index+m_range-1))='B';

b_index=m_index+m_range;
b_range=n_users*N*frame;
b_aux=reshape(b_index:(b_index+b_range-1),[n_users,N,frame]);
vtype(b_index:(b_index+b_range-1))='B';

x_index=b_index+b_range;
x_range=n_users*N*frame;
x_aux=reshape(x_index:(x_index+x_range-1),[n_users,N,frame]);
vtype(x_index:(x_index+x_range-1))='B';

z_index=x_index+x_range;
z_range=n_users*M*N*frame;
z_aux=reshape(z_index:(z_index+z_range-1),[n_users,M,N,frame]);
vtype(z_index:(z_index+z_range-1))='B';

g_index=z_index+z_range;
g_range=n_users*frame;
g_aux=reshape(g_index:(g_index+g_range-1),[n_users,frame]);
vtype(g_index:(g_index+g_range-1))='S';

f_index=g_index+g_range;
f_range=n_users*frame;
f_aux=reshape(f_index:(f_index+f_range-1),[n_users,frame]);
vtype(f_index:(f_index+f_range-1))='S';


c_index=f_index+f_range;
c_range=number_cells*frame;
c_aux=reshape(c_index:(c_index+c_range-1),[number_cells,frame]);
vtype(c_index:(c_index+c_range-1))='I';

a_p_index=c_index+c_range;
a_p_range=n_users*N*frame;
a_p_aux=reshape(a_p_index:(a_p_index+a_p_range-1),[n_users,N,frame]);
vtype(a_p_index:(a_p_index+a_p_range-1))='B';

a_n_index=a_p_index+a_p_range;
a_n_range=n_users*N*frame;
a_n_aux=reshape(a_n_index:(a_n_index+a_n_range-1),[n_users,N,frame]);
vtype(a_n_index:(a_n_index+a_n_range-1))='B';

vtype=vtype';
total_range=i_range+m_range+b_range+x_range+z_range+g_range+f_range+c_range+a_p_range+a_n_range %total number of variables

% Set of constraints (Ax<=b): A, b
total_constraints=n_users*M*N*frame+n_users*M*N*frame+n_users*M*N*frame+n_users*N*frame+n_users*N*frame+n_users*N*frame+n_users*N*frame+n_users*N*frame+n_users*N*frame+frame+frame+interference*N*frame+n_users*frame+n_users*frame+n_users*frame+number_cells*frame+n_users*frame+n_users*frame; %1)+2)+3)+4)+5)+6)+7)+8)+9)+10)+11)+12)+13)+14)+15)+16)+17)
A=spalloc(total_constraints,total_range,total_constraints*5); % number of non-zero elements: total_constraints*5 (aprox)
b=zeros(total_constraints,1);

% vf=[]
% vc=[]
% vd=[]
%A=sparse(vf,vc,vd)

nf=1;

%1) x(u,n,t)+m(u,m,t)-z(u,m,n,t)<=1 : n_users*M*N*frame
%2) z(u,m,n,t)-m(u,m,t)<=0 :  n_users*M*N*frame
%3) z(u,m,n,t)-x(u,n,t)<=0 :  n_users*M*N*frame

for u=1:n_users
    for m=1:M
        for n=1:N
            for t=1:frame
                %1)
                A(nf,x_aux(u,n,t))=1;
                A(nf,m_aux(u,m,t))=1;
                A(nf,z_aux(u,m,n,t))=-1;
                b(nf)=1;
                nf=nf+1;

                %2)
                A(nf,z_aux(u,m,n,t))=1;
                A(nf,m_aux(u,m,t))=-1;
                b(nf)=0;
                nf=nf+1;

                %3)
                A(nf,z_aux(u,m,n,t))=1;
                A(nf,x_aux(u,n,t))=-1;
                b(nf)=0;
                nf=nf+1;

            end
        end
    end
end

%4) ∑w b(u,w,t)-N*a_p(u,n,t)<=n : n_users*N*frame
%5) -∑w b(u,w,t)+(N+1)*a_p(u,n,t)<=N+1-n-1 : n_users*N*frame
%6) -∑w b(u,w,t)-N*a_n(u,n,t)<=-n : n_users*N*frame
%7) ∑w b(u,w,t)+N*a_n(u,n,t)<=N+n-1 : n_users*N*frame
%8) x(u,n,t)+a_p(u,n,t)<=1 : n_users*N*frame
%9) x(u,n,t)+a_n(u,n,t)<=1 : n_users*N*frame

for u=1:n_users
        for n=1:N
            for t=1:frame
                %4)
                for w=1:N
                    A(nf,b_aux(u,w,t))=1;
                end
                A(nf,a_p_aux(u,n,t))=-N;
                b(nf)=n;
                nf=nf+1;

                %5)
                for w=1:N
                    A(nf,b_aux(u,w,t))=-1;
                end
                A(nf,a_p_aux(u,n,t))=N+1;
                b(nf)=N-n;
                nf=nf+1;

                %6) 
                for w=1:N
                    A(nf,b_aux(u,w,t))=-1;
                end
                A(nf,a_n_aux(u,n,t))=-N;
                b(nf)=-n;
                nf=nf+1;

                %7)
                for w=1:N
                    A(nf,b_aux(u,w,t))=1;
                end
                A(nf,a_n_aux(u,n,t))=N;
                b(nf)=N-1+n;
                nf=nf+1;

                %8)
                A(nf,a_p_aux(u,n,t))=1;
                A(nf,x_aux(u,n,t))=1;
                b(nf)=1;
                nf=nf+1;

                %9)
                A(nf,a_n_aux(u,n,t))=1;
                A(nf,x_aux(u,n,t))=1;
                b(nf)=1;
                nf=nf+1;

            end
        end
end

%10) ∑c i(c,t)<=beams : frame
%11) ∑u p(u,t)=∑u ∑m,n p(u,m,n,t)=∑u ∑m,n P(m,n,u)*z(u,m,n,t)<=P_T : frame

for t=1:frame
    %10)
    for c=1:number_cells
        A(nf,i_aux(c,t))=1;
    end
    b(nf)=beams;
    nf=nf+1;

    %11)
    for m=1:M
        for n=1:N
            for u=1:n_users
                A(nf,z_aux(u,m,n,t))=P(m,n,u);
            end
        end
    end
    b(nf)=P_T;
    nf=nf+1;

end


%12) b(u(α),w,t) + b(u(β),w,t)<=1: interference*N*frames
for u_fila=1:size(Adj_u,1)
    for u_columna=1:size(Adj_u,2)
        if Adj_u(u_fila,u_columna)==1 && u_fila<u_columna
            %disp(strcat(num2str(fila),',',num2str(columna)))
            for w=1:N
                for t=1:frame
                    A(nf,b_aux(u_fila,w,t))=1;
                    A(nf,b_aux(u_columna,w,t))=1;
                    b(nf)=1;
                    nf=nf+1;
                end
            end
        end
    end
end


%13)  ∑m m(u,m,t)-i(c,t)<=0 : n_users*frame
%14)  ∑n x(u,n,t)-i(c,t)<=0 : n_users*frame  -> old: ∑w b(u,w,t)-i(c,t)<=0
%14v)  ∑n b(u,n,t)-N*i(c,t)<=0 : n_users*frame  
%16)  g(u,t-1)-d(u,t)-g(u,t)= g(u,t-1)-∑m,n D(m,n)*z(u,m,n,t)-g(u,t)<=-R(u,t) : n_users*frame
%17)  -g(u,t-1)+d(u,t)-f(u,t)= -g(u,t-1)+∑m,n D(m,n)*z(u,m,n,t)-f(u,t)<=R(u,t) : n_users*frame

for u=1:n_users
    for t=1:frame
        %13)
        for m=1:M
            A(nf,m_aux(u,m,t))=1;
        end
        A(nf,i_aux(find(UpC(:,u)),t))=-1;
        b(nf)=0;
        nf=nf+1;

        %14)
        for n=1:N
            A(nf,x_aux(u,n,t))=1;
        end
        A(nf,i_aux(find(UpC(:,u)),t))=-1;
        b(nf)=0;
        nf=nf+1;

        %14v)
        for n=1:N
            A(nf,b_aux(u,n,t))=1;
        end
        A(nf,i_aux(find(UpC(:,u)),t))=-N;
        b(nf)=0;
        nf=nf+1;

        %16)
        for m=1:M
            for n=1:N
                A(nf,z_aux(u,m,n,t))=-D(m,n);
            end
        end
        A(nf,g_aux(u,t))=-1;
        if t~=1
            A(nf,g_aux(u,t-1))=1;
        end
        b(nf)=-R(u,t);
        nf=nf+1;

        %17)
        for m=1:M
            for n=1:N
                A(nf,z_aux(u,m,n,t))=+D(m,n);
            end
        end
        A(nf,f_aux(u,t))=-1;
        if t~=1
            A(nf,g_aux(u,t-1))=-1;
        end
        b(nf)=R(u,t);
        nf=nf+1;

    end
end



%15)  c(c,t-1)-frame*i(c,t)-c(c,t)+g(u,t)/∑t (Ru,t)<=0 : number_cells*frame
for c=1:number_cells
    for t=1:frame
        if t==1
            A(nf,i_aux(c,t))=-frame;
            A(nf,c_aux(c,t))=-1;
            userpercell=find(UpC(c,:));
            for u=1:length(userpercell)
                A(nf,g_aux(userpercell(u),t))=1/sum(R(userpercell(u),1:t));
            end
            b(nf)=0;
            nf=nf+1;   
        else
            A(nf,c_aux(c,t-1))=1;
            A(nf,i_aux(c,t))=-frame;
            A(nf,c_aux(c,t))=-1;
            userpercell=find(UpC(c,:));
            for u=1:length(userpercell)
                A(nf,g_aux(userpercell(u),t))=1/sum(R(userpercell(u),:));
            end            
            b(nf)=0;
            nf=nf+1;
        end
    end
end

%MIP Start
mip_start=[i_start,m_start,b_start,x_start,z_start];
%mip_start=[i_start,m_start,b_start,x_start,z_start,NaN(1,total_range-length([i_start,m_start,b_start,x_start,z_start]))];

[row,col,val] = find(A);

toc
end
