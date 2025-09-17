# MILP-based-Satellite-Beam-Hopping-BH-Algorithm

This code is © Samuel M. Zamacola, 2025, and it is made available under the GPL license enclosed with the software.

Over and above the legal restrictions imposed by this license, if you use this software for an academic publication then you are obliged to provide proper attribution to the paper that describes it:
+ S. M. Zamacola, N. Pachler,  R. M. Rodríguez-Osorio and B. G. Cameron, ‘Joint Illumination, Power and Band Allocation for Multi-beam LEO Satellites with Beam-Hopping Using Mixed Integer Linear Programming’, .....

DEF.: "The MILP based BH algorithm is an analytical algorithm where illumination scheduling (BH) and resource allocation (bandwindth+power) is jointly performed in response to users' distribution and traffic demand. The overall flowchart is presented in Simulation_Flow.pdf. The diagram is divided in two parts: (1) pre-computation modules, where input system and scenario related variables are computed; (2) MILP formulation module, where decision variables, constraints and objective function are defined; and (3) the MILP solver, where the BH algorithm dynamically allocates illumination [Ill], power [P] and bandwidth [B] resources for each time slot, xomputing the parameters of interest (PoI): UC/EC/TTS. For more information, check the reference paper."

Two resolution methodologies are presented: Full MILP and Time-Split MILP. The full MILP, considers the high dimensional problem, optimizing resource allocation: illumination + resource allocation (bandwidth and power) for all the time instants (frames). To reduce the computational complexity, an alternative time-split MILP is proposed by partitioning the problem into smaller time-based subproblems, balancing execution time and overall performance.

* INPUT: satellite altitude (h_sat), minimum elevation angle (el_min), number of rings within the satellite's FoV (rings), frequency (freq), total ilumination slots (frame), slot duration (frame_dur), number of colours (colours), simultaneous beams (beams), frequency bins(N), number of users (n_users), distribution model of users (traffic_model): random/linear/hotspot, total RF power (P_T), total bandwidth per colour (B_T).

* OUTPUT: Illumination [Ill] + Resource Allocation (power [P] and bandwidth [B]) + FoM (Unserved Capacity (UC), Extra Served Capacity (EC), Time To Serve (TTS)).
 
Included Files:
 
+ User class (u.m): Each user is characterized by a given location, a type of station, a given traffic demand and
counts at UE level with some gain and noise characteristics. These are the static attributes of
the class, but there are other attributes that are going to be dynamically computed through
methods based on these prior attributes. Depending on the subsatellite point location and the
position of the satellite the following parameters are calculated: the slant range, the elevation
angle from the UE to the satellite, and the nadir angles, with respect to the satellite, and relative
to its cell center. Dynamic attributes, those that are going to be filled in each frame iteration and
therefore are time dependent, are used to store the intermediate and output variables of the
iteration. These are: assigned power and bandwidth resources for link budget computation, C/N
which is link budget’s output, and then served, pending and extra traffic counters.
 
+ Cell class (c.m): Cells are labelled by a given identifier or cell number, count with a given center and have a
radius that is going to be dependant of the total footprint area and the number of chosen rings.
Cell objects will include a position based computed attribute in which the user objects presented
above are stored. This is where the composition property comes into action, as cells count with
users located within their area, and the organisation of these is much simpler by incorporating
a list of user objects, if any, in each of the cells. Additionally, the nadir angle from the cell centre
and the list of interfering cells is computed based on the initial scheme. During each frame
iteration the aggregated traffic requirement is calculated based on the user list within the given
cell. Based on the cells that are selected for illumination, these are going to be classified in each
time instant as active or inactive. If that is the case, the aggregated power and bandwidth
resources assigned at cell level are computed. In case the cell is active, the assigned colour is
also stored as a dynamic attribute.
 
+ main.m: the main executable file where BH algorithm is executed. 1) Instance Generation. 2) BH Pre-computation. 3) Demand-based (DB) Initialization 4) MILP BH Definition. 5) MILP BH Calculation. 6) Instance Result Saving.

  +1) main.m: input parameters are defined for each type of scenario, and the resolution methodology is selected, method='full' for Full MILP and method='time-split' for Time-split MILP.
  
  +2) pre_BH_computation.m: parameter pre-computation is performed, (traffic_model->Traffic_Distribution, cell_scenario_model->Cell_Scenario) that determines (user per cell allocation (UpC) together with adjacency matrixes (AdjU and AdjC)) allowing to compute requested capacity (R) potential offered capacity (D), required power (P) matrixes.
 
    + Traffic_Distribution.m: based on the selected user distribution type, the generation of users is performed in the file, by defining UE related specifications.
 
    + Cell_Scenario.m: based on the number of rings that are intended to be allocated within the satellite's FoV, cells are generated in the file.

  +3) FOM_calculation_demand_based_band_slots_fixed_MODCODs.m: Demand-based BH algorithm is executed to include the solution as the first solution in MILP  initialization. 
 
  +4) BH.m: initialization of MILP decision variables, constraints and objective function.
 
  +5) gurobi_execution.m: call to Gurobi solver through its Python API (BH_MIPStart_lb_ub.py) where initialization (mip_start), objective function (obj) and constraints (A->(row,col,val) and b) are passed. . NOTE: For its execution Gurobu solver is required to be installed.
  
  +6) solution_plot_saving.m: obtained [Ill], [B] and [P] decision variables togwrhwe with main PoIs are saved -> 'NSGA_II_[', instance,'_', num2str(scenario), ']_result_','MAXFEs', num2str(maxFEs), '.P',num2str(popSize),'.i',init,'.cR', num2str(c_rate_LEO),'.mType', m_type_LEO,'.mRate', num2str(m_rate_LEO),'.lsType', ls_type, '.lsRate',num2str(ls_rate),'.users', num2str(n_users),'.r',num2str(run),'.mat': Instance result saving (NSGA-II + DB).
		+ Full MILP: '[',num2str(scenario+1),'s_',num2str(run),'r]_',num2str(betta),'betta_',num2str(frame),'frames_',num2str(n_users),'u','.mat' <- ("UC_perc_M","EC_perc_M","TTS_M","UC_perc_db","EC_perc_db","TTS_db","O1","O2","O3","OBJETIVO","normalization_UC","normalization_EC","normalization_time","c_scenario","B_T","b_slots","P_T","Adj_c","Adj_u","UpC","n_users","rings","beams","theta","colours","frame","frame_dur","freq","Ill_out","B_out","P_out","MS","GS","FS","CS","ZS","XS","PS","DS","data","solutions","betta")
		+ Time-split MILP:'[',num2str(scenario+1),'s_',num2str(run),'r_split]_',num2str(betta),'betta_',num2str(t_current+9),'frame_',num2str(n_users),'u','.mat' <- ("UC_perc_M","EC_perc_M","TTS_M","UC_perc_db","EC_perc_db","TTS_db","O1","O2","O3","OBJETIVO","normalization_UC","normalization_EC","normalization_time","c_scenario","B_T","b_slots","P_T","Adj_c","Adj_u","UpC","n_users","rings","beams","theta","colours","frame","frame_dur","freq","Ill_out","B_out","P_out","MS","GS","FS","CS","ZS","XS","PS","DS","data","solutions","betta")
		
		
NOTE: Demand-based (DB) BH algorithm is available at: https://github.com/samuelmartinezupm/Demand-based-Satellite-Beam-Hopping-BH-Algorithm 
NOTE: NSGA-II BH algorithm is available at: https://github.com/samuelmartinezupm/NSGA-II-based-Satellite-Beam-Hopping-BH-Algorithm

## Contact
For questions, issues, or contributions, please contact Samuel M. Zamacola at samuel.martinez@upm.es.
