%% GUROBI EXECUTION:
function [result,data,solutions]=gurobi_execution(n_users,number_cells,frame,betta,normalization_UC,normalization_EC,normalization_time,total_range,total_constraints,val,row,col,b,vtype,MIPGap,mip_start,g_aux,f_aux,c_aux)
    
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
    
    % With MIP Start
    cd (strcat(PWD,'/Analytical/'))
    [result,data_,solutions_]=(pyrunfile("BH_MIPStart.py",["result" "data" "solutions"],total_range=total_range,total_constraints=total_constraints,val=val,row=row,col=col, b=b, obj=obj, vtype=vtype', MIPGap=MIPGap, mip_start=mip_start));
    
    result=double(result);
    data=[];
    solutions=[];
    for idx=1:length(data_)
        data=[data; double(data_{idx})];
        solutions=[solutions; double(solutions_{idx})];    
    end
end
