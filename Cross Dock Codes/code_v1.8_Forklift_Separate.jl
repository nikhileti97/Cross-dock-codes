#updated objective to makespan

using JuMP, Gurobi;

CDSSP = Model(optimizer_with_attributes(Gurobi.Optimizer, "TIME_LIMIT" => 4.95*3600));


@variable(CDSSP, 0 <= Entry1Inbound[i = 1:n_inbound] <= HorizEnd); # e_I => Entry1Inbound
@variable(CDSSP, 0 <= Leave1Inbound[i = 1:n_inbound] <= HorizEnd); # l_I => Leave1Inbound
@variable(CDSSP, 0 <= Entry1Outbound[j = 1:n_outbound] <= HorizEnd); # e_O => Entry1Outbound
@variable(CDSSP, 0 <= Leave1Outbound[j = 1:n_outbound] <= HorizEnd); # l_O => Leave1Outbound

#@variable(CDSSP, 0 <= ExchLeave1Outbound[j = 1:n_outbound, g = 1:n_inbound] <= HorizEnd); # l_O => Leave1Outbound
@variable(CDSSP, 0 <= Tardiness1[j = 1:n_outbound]); # gamma => Tardiness1
@variable(CDSSP, 0 <= MaxPwl1[j = 1:n_outbound]);

@variable(CDSSP, 0 <= FirstLoad[i = 1:n_inbound, j = 1:n_outbound] <= 1, Bin); # OutAssign => y
@variable(CDSSP, 0 <= LastLoad[i = 1:n_inbound, j = 1:n_outbound] <= 1, Bin); # OutAssign => y
@variable(CDSSP, 0 <= FirstUnload[i = 1:n_inbound] <= 1, Bin); # OutAssign => y
@variable(CDSSP, 0 <= LastUnload[i = 1:n_inbound] <= 1, Bin); # OutAssign => y

@variable(CDSSP, 0 <= UnFork[i = 1:n_inbound, f = 1:n_forklifts] <= 1, Bin); # Unloading Forklift_Assignment => x
@variable(CDSSP, 0 <= LFork[i = 1:n_inbound, j = 1:n_outbound, f = 1:n_forklifts] <= 1, Bin); # Forklift_Assignment => x
@variable(CDSSP, 0 <= gamma[i = 1:n_inbound, q = 1:n_inbound] <= 1, Bin);  # alpha => InbPrecedence
@variable(CDSSP, 0 <= alphaij[i = 1:n_inbound, j = 1:n_outbound,  i1 = 1:n_inbound] <= 1, Bin); # Foklift Usage => Unloading fisrt, loading later
@variable(CDSSP, 0 <= betaij[i = 1:n_inbound, j = 1:n_outbound,  i1 = 1:n_inbound] <= 1, Bin); # Forklift Usage => Loading First, unloading later
@variable(CDSSP, 0 <= deltajj[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound,  j1 = 1:n_outbound] <= 1, Bin); # Foklift Usage => loading fisrt, loading later

@variable(CDSSP, 0 <= St[i = 1:n_inbound] ); # St_if => Forklift Unloading Start
#@variable(CDSSP, 0 <= Lt[i = 1:n_inbound, f = 1:n_forklifts] <= HorizEnd); # Lt_if => Forklift Unloading Finish
@variable(CDSSP, 0 <= Sl[i = 1:n_inbound, j = 1:n_outbound] ); # Sl_ijf => Forklift loading Starts
@variable(CDSSP, 0 <= Ll[i = 1:n_inbound, j = 1:n_outbound] ); # Sl_ijf => Forklift loading Starts





@constraint(CDSSP, InbAssigFork[i = 1:n_inbound], sum(UnFork[i, f] for f = 1:n_forklifts) == 1); # Unloading Forklift assignment
#@constraint(CDSSP, OutAssigCon[j = 1:n_outbound], sum(OutAssign[j, g] for g = 1:n_stack) == 1);

@constraint(CDSSP, EnterStripCon[i = 1:n_inbound], Entry1Inbound[i] >= Arrival[i]);
@constraint(CDSSP, EnterUnFork[i = 1:n_inbound], St[i] >= Entry1Inbound[i]); # Starting time of forklift unloading

#@constraint(CDSSP, LeaveStripCon[i = 1:n_inbound, k = 1:n_strip], Leave1Inbound[i, k] - Entry1Inbound[i, k] >= InbAssign[i, k]*sum(Supply[i, p]*LoadTime[p] for p = 1:n_products));
#@constraint(CDSSP, LeaveUnFork[i = 1:n_inbound, f = 1:n_forklifts], Lt[i, f] - St[i, f] >= UnFork[i, f]*sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products)); # End Time for forklift unloading

#@constraint(CDSSP, LeaveStripCon1[i = 1:n_inbound], sum(Leave1Inbound[i, k] for k = 1:n_strip) >= sum(Lt[i,f] for f = 1:n_forklifts));

@constraint(CDSSP, LeaveStripCon1[i = 1:n_inbound], Leave1Inbound[i] >= St[i]+sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products));

@constraint(CDSSP, EntrStripCon4[j = 1:n_inbound, q = 1:n_inbound], Entry1Inbound[j] >= (Leave1Inbound[q]+Change_time)*Inb_Precedence[q,j]);                                        

@constraint(CDSSP, LoadFirst[i = 1:n_inbound, j = 1:n_outbound], FirstLoad[i,j] + sum(alphaij[i, j, q] for q = 1:n_inbound) + sum(deltajj[i1,j1,i,j] for i1 = 1:n_inbound for j1 = 1:n_outbound if i1 != i || j1 != j) == 1);     
                                                                                                
@constraint(CDSSP, LoadLast[i = 1:n_inbound, j = 1:n_outbound], LastLoad[i,j] + sum(betaij[i, j, q] for q = 1:n_inbound) + sum(deltajj[i,j,i1,j1] for i1 = 1:n_inbound for j1 = 1:n_outbound if i1 != i || j1 != j) == 1);          
                                                                                                                                
@constraint(CDSSP, UnLoadFirst[i = 1:n_inbound], FirstUnload[i] + sum(betaij[i1, j, i] for i1 = 1:n_inbound for j = 1:n_outbound) + sum(gamma[i1,i] for i1 = 1:n_inbound if i1 != i) == 1);     
                                                                                                
@constraint(CDSSP, UnLoadLast[i = 1:n_inbound], LastUnload[i] + sum(alphaij[i1, j, i] for i1 = 1:n_inbound for j = 1:n_outbound) + sum(gamma[i,i1] for i1 = 1:n_inbound if i1 != i) == 1);                                                                                                                                          
                                                                                                                                
                                                                                                                                
@constraint(CDSSP, First1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound, f = 1:n_forklifts], FirstLoad[i, j] + FirstUnload[q] <= 3 - LFork[i,j,f] - UnFork[q,f]);

@constraint(CDSSP, First2[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound,  f = 1:n_forklifts; i1 != i || j1 != j], FirstLoad[i, j] + FirstLoad[i1,j1] <= 3 - LFork[i,j,f] - LFork[i1,j1,f]);  
                                                                                                                                
@constraint(CDSSP, First3[i = 1:n_inbound, q = 1:n_inbound, f = 1:n_forklifts; q != i], FirstUnload[i] + FirstUnload[q] <= 3 - UnFork[i,f] - UnFork[q,f]);                                                                                                                                
@constraint(CDSSP, Last1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound, f = 1:n_forklifts], LastLoad[i, j] + LastUnload[q] <= 3 - LFork[i,j,f] - UnFork[q,f]);

@constraint(CDSSP, Last2[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound,  f = 1:n_forklifts; i1 != i || j1 != j], LastLoad[i, j] + LastLoad[i1,j1] <= 3 - LFork[i,j,f] - LFork[i1,j1,f]);  
                                                                                                                                
@constraint(CDSSP, Last3[i = 1:n_inbound, q = 1:n_inbound, f = 1:n_forklifts; q != i], LastUnload[i] + LastUnload[q] <= 3 - UnFork[i,f] - UnFork[q,f]);          

#@constraint(CDSSP, UnLoadLoad1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound], Sl[i,j] >= C_time + St[q] + sum(Supply[q, p]*UnLoadTime[p] for p = 1:n_products)- HorizEnd*(1 - alphaij[i, j, q])); # Loading followed by unloading
@constraint(CDSSP, UnLoadLoad1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound], alphaij[i,j,q] => {Sl[i,j] >= C_time + St[q] + sum(Supply[q, p]*UnLoadTime[p] for p = 1:n_products)}); # Loading followed by unloading

@constraint(CDSSP, UnLoadLoad1a[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound], betaij[i,j,q] => {St[q] >= C_time + Ll[i,j]}); # Loading followed by unloading

@constraint(CDSSP, UnLoadLoad2[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound,f = 1:n_forklifts], alphaij[i,j,q] - 1 <=  LFork[i, j, f] - UnFork[q,f]); # Loading followed by unloading
                                                                                                                                
@constraint(CDSSP, UnLoadLoad3[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound,f = 1:n_forklifts], 1 - alphaij[i,j,q] >=  LFork[i, j, f] - UnFork[q,f]); # Loading followed by unloading                                                                                             
                                                                                                                                
@constraint(CDSSP, UnLoadLoad4[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound,f = 1:n_forklifts], betaij[i,j,q] - 1 <=  LFork[i, j, f] - UnFork[q,f]); # Loading followed by unloading

@constraint(CDSSP, UnLoadLoad5[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound,f = 1:n_forklifts], 1 - betaij[i,j,q] >=  LFork[i, j, f] - UnFork[q,f]); # Loading followed by unloading                                                                                               
                                                                                                                                
#@constraint(CDSSP, LoadLoad1[i = 1:n_inbound, i1 = 1:n_inbound, j = 1:n_outbound, j1 = 1:n_outbound, f = 1:n_forklifts; i!=i1 || j !=j1], Sl[i1,j1] >= C_time + Ll[i, j] - HorizEnd*(1 - deltajj[i, j, i1, j1])); # Loading followed by loading
@constraint(CDSSP, LoadLoad1[i = 1:n_inbound, i1 = 1:n_inbound, j = 1:n_outbound, j1 = 1:n_outbound, f = 1:n_forklifts; i!=i1 || j !=j1], deltajj[i,j,i1,j1] => {Sl[i1,j1] >= C_time + Ll[i, j]}); # Loading followed by loading

@constraint(CDSSP, LoadLoad2[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound, f = 1:n_forklifts; i!=i1 || j !=j1], deltajj[i,j,i1,j1] - 1 <=  LFork[i, j, f] - LFork[i1,j1,f]); # Loading followed by unloading
                                                                                                                                
@constraint(CDSSP, LoadLoad3[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound, f = 1:n_forklifts; i!=i1 || j !=j1], 1 - deltajj[i,j,i1,j1] >=  LFork[i, j, f] - LFork[i1,j1,f]); # Loading followed by unloading                                                                                             

@constraint(CDSSP, LoadLoad4[i = 1:n_inbound, j = 1:n_outbound], deltajj[i,j,i,j] == 0); # UnLoading followed by Loading

@constraint(CDSSP, AuxUnLoadUnLoad1[i = 1:n_inbound, q = 1:n_inbound; i != q], 1 - Inb_Precedence[i,q] >= gamma[q,i]); # Unloading followed by unloading
@constraint(CDSSP, AuxUnLoadLoad1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound; i != q], 1 - Inb_Precedence[i,q] >= alphaij[i,j,q]);
@constraint(CDSSP, AuxUnLoadLoad2[i = 1:n_inbound, j = 1:n_outbound], betaij[i,j,i] == 0);
@constraint(CDSSP, AuxLoadLoad1[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound; i!=i1], deltajj[i,j,i1,j1]  <=  1 - Inb_Precedence[i1,i]);# Loading followed by unloading
@constraint(CDSSP, AuxLoadLoad2[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound; j!=j1], deltajj[i,j,i1,j1]  <=  1 - Out_Precedence[j1,j]);# Loading followed by unloading

#@constraint(CDSSP, UnLoadUnLoad1[i = 1:n_inbound, q = 1:n_inbound, f = 1:n_forklifts; i != q], St[q] >= C_time + St[i] + sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products) - HorizEnd*(1 - gamma[i, q])); # Unloading followed by unloading
@constraint(CDSSP, UnLoadUnLoad1[i = 1:n_inbound, q = 1:n_inbound, f = 1:n_forklifts; i != q], gamma[i,q] => {St[q] >= C_time + St[i] + sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products)}); # Unloading followed by unloading

@constraint(CDSSP, UnLoadUnLoad2[i = 1:n_inbound, q = 1:n_inbound,f = 1:n_forklifts; q != i], gamma[i,q] - 1 <=  UnFork[i, f] - UnFork[q,f]); # Loading followed by unloading
                                                                                                                                
@constraint(CDSSP, UnLoadUnLoad3[i = 1:n_inbound, q = 1:n_inbound,f = 1:n_forklifts; q != i], 1 - gamma[i,q] >=  UnFork[i, f] - UnFork[q,f]); # Loading followed by unloading        

@constraint(CDSSP, UnLoadUnLoad4[i = 1:n_inbound], gamma[i,i] == 0); # UnLoading followed by unloading

@constraint(CDSSP, LoadAssign[i = 1:n_inbound, j = 1:n_outbound], sum(LFork[i, j, f] for f = 1:n_forklifts) == Exch[i,j]); # Assign one loading task to one forklift, otherwise divide loads of 'rijp' between multiple forklifts, making the problem complicated  


@constraint(CDSSP, StartLoadFork[i = 1:n_inbound, j = 1:n_outbound], Sl[i, j] >= St[i]*Exch[i,j] + sum(Supply[i,p]*UnLoadTime[p] for p in 1:n_products)/2*Exch[i,j]); # Starting time of forklift unloading

@constraint(CDSSP, EndLoadFork[i = 1:n_inbound, j = 1:n_outbound], Ll[i, j] >= Sl[i, j] + sum(Load_Tran[i,j,p]*LoadTime[p] for p in 1:n_products)); #End time of forklift unloading

@constraint(CDSSP, LeaveStackCon1[i = 1:n_inbound, j = 1:n_outbound], Leave1Outbound[j] >= Ll[i, j]);
# ensures minimum stay of outbound truck;
# truck should be docked before starting any loading job; truck latest docking time given by precedence constraints;
@constraint(CDSSP, LeaveStackCon2[i = 1:n_inbound, j = 1:n_outbound], Leave1Outbound[j] >= Entry1Outbound[j] + sum(Demand[j, p]*LoadTime[p] for p = 1:n_products));

@constraint(CDSSP, LeaveStackCon3[i = 1:n_inbound, j = 1:n_outbound], Sl[i, j] >= Entry1Outbound[j]);
                                                                                                
@constraint(CDSSP, LeaveStackCon4[j = 1:n_outbound, q = 1:n_outbound], Entry1Outbound[j] >= (Leave1Outbound[q]+Change_time)*Out_Precedence[q,j]);
                                                                                                
@constraint(CDSSP, DeadCon[j = 1:n_outbound], Leave1Outbound[j] <= DeadLines[j] + Tardiness1[j]);
#@constraint(CDSSP, DeadCon[j = 1:n_outbound], Leave1Outbound[j] <= Leav_OT[j] + Tardiness1[j]);
#@@constraint(CDSSP, DeadCon1[j = 1:n_outbound], LeaveOutbound[j] <= value(LeaveOutbound[j]) + Earliness[j]);                                                                                                
                                                                                                
                                                                                                
@constraint(CDSSP, LinearizeCon11[j = 1:n_outbound], MaxPwl1[j] >= Multiplier[1]*Penalty[j]*(Tardiness1[j]));
@constraint(CDSSP, LinearizeCon22[j = 1:n_outbound], MaxPwl1[j] >= Multiplier[2]*Penalty[j]*(Tardiness1[j]) + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]));
@constraint(CDSSP, LinearizeCon33[j = 1:n_outbound], MaxPwl1[j] >= Multiplier[3]*Penalty[j]*(Tardiness1[j]) +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]));
                                                                                                
#@constraint(CDSSP, InbLastEntry[i = 1:n_inbound, k = 1:n_strip], Entry1Inbound[i, k] <= InbAssign[i,k]*HorizEnd);
#@constraint(CDSSP, MakespanCon[j = 1:n_outbound, g = 1:n_stack], Makespan[j] >= Leave1Outbound[j,g]);

@objective(CDSSP, Min, sum(MaxPwl1[j] for j = 1:n_outbound));
                                                                                                
#optimize!(CDSSP) 