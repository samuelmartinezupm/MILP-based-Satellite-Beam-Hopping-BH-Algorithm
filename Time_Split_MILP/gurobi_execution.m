%% GUROBI EXECUTION:
function [result,data,solutions]=gurobi_execution(n_users,number_cells,frame,betta,normalization_UC,normalization_EC,normalization_time,total_range,total_constraints,A,b,vtype,MIPGap,mip_start,g_aux,f_aux,c_aux,lb,ub,t_current,t_contraint_label,t_range_boolean)
    
    global PWD;

    % Weighting parameter estimation:
    weight_traffic_EC=(1-betta)/2;
    weight_traffic_UC=betta;
    weight_timing=(1-betta)/2;
    
   % Objective function: 
    obj=zeros(1,total_range);

   for u=1:n_users
        for t=1:frame
            obj(1,g_aux(u,t))=weight_traffic_UC/normalization_UC;
            obj(1,f_aux(u,t))=weight_traffic_EC/normalization_EC;
        end
    end
    
    for c=1:number_cells
        for t=1:frame
            obj(1,c_aux(c,t))=weight_timing/normalization_time;
        end
    end
    
    % % Without MIP Start
    % save('BH_input_data.mat', "total_range","total_constraints","val","row","col","b","obj","vtype","MIPGap");
    % [result,data_,solutions_]=(pyrunfile("BH.py",["result" "data" "solutions"],total_range=total_range,total_constraints=total_constraints,val=val,row=row,col=col, b=b, obj=obj, vtype=vtype', MIPGap=MIPGap));
    % %result.x=double(pyrunfile("BH.py","result","data","solutions",total_range=total_range,total_constraints=total_constraints,val=val,row=row,col=col, b=b, obj=obj, vtype=vtype', MIPGap=MIPGap))
    
    % A=A((t_contraint_label<=(t_current+9)),:); % Filter the contraints that must be considered in the execution.
    % [row,col,val] = find(A);
    % 
    % b=b((t_contraint_label<=(t_current+9))); % Filter the contraints that must be considered in the execution.
    % 
    % total_constraints_to_consider=sum(t_contraint_label<=(t_current+9));

    % TEST INIT
    t_constraint_boolean=(t_contraint_label<=(t_current+9)).*(t_contraint_label>=(t_current));

    A=A(t_constraint_boolean==1,:);
    [row,col,val] = find(A);
   
    b=b(t_constraint_boolean==1,:);

    total_constraints_to_consider=sum(t_constraint_boolean);
    % TEST FINISH

    % With MIP Start
    mip_start(t_range_boolean~=1)=NaN; % Cancel the variables above the study!
    cd (strcat(PWD,'/Analytical/')) 
    [result,data_,solutions_]=(pyrunfile("BH_MIPStart_lb_ub.py",["result" "data" "solutions"],total_range=total_range,total_constraints=total_constraints_to_consider,val=val,row=row,col=col, b=b, obj=obj, vtype=vtype', MIPGap=MIPGap, mip_start=mip_start, lb=lb, ub=ub));
    
    % % % 
    % % % Direct in MATLAB:
    % addpath 'C:\gurobi1100\win64\matlab'
    % model.vtype=vtype';
    % %A=sparse(row,col,val);   
    % model.A           = A;
    % model.rhs         = b;
    % model.sense       = '<';
    % model.obj = obj';
    % model.modelsense  = 'min';
    % model.lb=lb;
    % model.ub=ub;
    % params.outputflag = 1;
    % params.MIPGap= 0.4;
    % %iis = gurobi_iis(model,params); 
    % model.start=mip_start;
    % result = gurobi(model, params);
    % % disp(result)
    
    result=double(result);
    data=[];
    solutions=[];
    for idx=1:length(data_)
        data=[data; double(data_{idx})];
        solutions=[solutions; double(solutions_{idx})];    
    end
end
