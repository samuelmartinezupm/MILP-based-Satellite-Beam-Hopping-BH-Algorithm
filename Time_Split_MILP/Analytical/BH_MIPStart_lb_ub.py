#python3.10

import numpy as np
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
import scipy.io

def BH(total_range,total_constraints,val,row,col, b, obj, vtype, MIPGap, mip_start,lb,ub):

    # Create a new model
    model = gp.Model("named_model")

    val=np.array(val,dtype=float)
    row=np.array(row,dtype=int)
    col=np.array(col,dtype=int)
    # Create variables
    A = sp.csr_matrix((val, (row-1, col-1)), shape=(int(total_constraints), int(total_range)))

    vtype=[char for char in vtype]
    #vtype=np.array(vtype,dtype=str)

    lb=np.array(lb,dtype=float)
    ub=np.array(ub,dtype=float)

    x=model.addMVar(shape=int(total_range), vtype=vtype, lb=lb, ub=ub, name="x")
    
    obj=np.array(obj,dtype=float)

    mip_start=np.array(mip_start,dtype=float)

    # Set objective
    model.setObjective(obj @ x, gp.GRB.MINIMIZE)

    b=np.array(b,dtype=float)

    # Add all constraints
    model.addConstr(A @ x <= b)

    model.setParam('OutputFlag', 1)
    model.setParam('MIPGap', MIPGap)
   # Adding MIPStart
    if mip_start is not None:
       for i in range(len(mip_start)):
           if np.isnan(mip_start[i]):
                x[i].start = GRB.UNDEFINED
           else:
               x[i].start =mip_start[i] 
    # model.computeIIS()

    # Optimize model
    model._vars = x
    model._data = []
    model._solutions = []
    model.optimize(cb)

    print(x.X)
    print('Obj: %g' % model.objVal)

    return x.X, model._data, model._solutions


def cb(model, where):
    #if where == gp.GRB.Callback.MIP:
    #    cur_time = model.cbGet(gp.GRB.Callback.RUNTIME)
    #    cur_obj = model.cbGet(gp.GRB.Callback.MIP_OBJBST)
    #    cur_bd = model.cbGet(gp.GRB.Callback.MIP_OBJBND)
    #    model._data.append([cur_time, cur_obj, cur_bd])

    if where == GRB.Callback.MIPSOL:

        model._solutions.append(model.cbGetSolution(model._vars))

        cur_time = model.cbGet(gp.GRB.Callback.RUNTIME)
        cur_obj = model.cbGet(gp.GRB.Callback.MIPSOL_OBJBST)
        cur_bd = model.cbGet(gp.GRB.Callback.MIPSOL_OBJBND)
        model._data.append([cur_time, cur_obj, cur_bd])
       

    #if where == gp.GRB.Callback.MIPNODE:
    #    cur_time = model.cbGet(gp.GRB.Callback.RUNTIME)
    #    cur_obj = model.cbGet(gp.GRB.Callback.MIPNODE_OBJBST)
    #    cur_bd = model.cbGet(gp.GRB.Callback.MIPNODE_OBJBND)
    #    model._data.append([cur_time, cur_obj, cur_bd])


#mat = scipy.io.loadmat('BH_input_data.mat', squeeze_me=True)
#result,data,solutions=BH(mat['total_range'],mat['total_constraints'],mat['val'],mat['row'],mat['col'], mat['b'], mat['obj'], mat['vtype'], mat['MIPGap'],mat['mip_start'],mat['lb'],mat['ub'],mat['t_contraint_label'],mat['t_current'])
result,data,solutions=BH(total_range,total_constraints,val,row,col, b, obj, vtype, MIPGap,mip_start,lb,ub)


