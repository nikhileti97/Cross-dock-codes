#updated objective to makespan

using JuMP, Gurobi;

CDSP = Model(optimizer_with_attributes(Gurobi.Optimizer, "TIME_LIMIT" => 4.5*60*60));

@variable(CDSP, 0 <= InbAssign[i = 1:n_inbound, k = 1:n_strip] <= 1, Bin); # InbAssign => x
@variable(CDSP, 0 <= OutAssign[j = 1:n_outbound, g = 1:n_stack] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= FirstInAssign[j = 1:n_inbound] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= LastInAssign[j = 1:n_inbound] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= FirstOutAssign[j = 1:n_outbound] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= LastOutAssign[j = 1:n_outbound] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= InbPrecedence[i = 1:n_inbound, q = 1:n_inbound] <=1, Bin);  # alpha => InbPrecedence
@variable(CDSP, 0 <= OutPrecedence[j = 1:n_outbound, r = 1:n_outbound]<=1, Bin); # beta => OutPrecedence
@variable(CDSP, 0 <= Exch[i = 1:n_inbound, j = 1:n_outbound] <= 1, Bin); # v => TransRoute
@variable(CDSP, 0 <= LoadTrans[i = 1:n_inbound, j = 1:n_outbound, p = 1:n_products], Int); # r => LoadTrans
@variable(CDSP, 0 <= EntryInbound[i = 1:n_inbound] <= HorizEnd); # e_I => EntryInbound
@variable(CDSP, 0 <= LeaveInbound[i = 1:n_inbound] <= HorizEnd); # l_I => LeaveInbound
@variable(CDSP, 0 <= EntryOutbound[j = 1:n_outbound] <= HorizEnd); # e_O => EntryOutbound
@variable(CDSP, 0 <= LeaveOutbound[j = 1:n_outbound] <= HorizEnd); # l_O => LeaveOutbound

#@variable(CDSP, 0 <= ExchLeaveOutbound[j = 1:n_outbound, g = 1:n_inbound] <= HorizEnd); # l_O => LeaveOutbound
@variable(CDSP, 0 <= Tardiness[j = 1:n_outbound] <= HorizEnd); # gamma => Tardiness
@variable(CDSP, 0 <= MaxPwl[j = 1:n_outbound]); # z => MaxPwl


@variable(CDSP, 0 <= FirstLoad[i = 1:n_inbound, j = 1:n_outbound] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= LastLoad[i = 1:n_inbound, j = 1:n_outbound] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= FirstUnload[i = 1:n_inbound] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= LastUnload[i = 1:n_inbound] <= 1, Bin); # OutAssign => y

@variable(CDSP, 0 <= UnFork[i = 1:n_inbound, f = 1:n_forklifts] <= 1, Bin); # Unloading Forklift_Assignment => x
@variable(CDSP, 0 <= LFork[i = 1:n_inbound, j = 1:n_outbound, f = 1:n_forklifts] <= 1, Bin); # Forklift_Assignment => x
@variable(CDSP, 0 <= gamma[i = 1:n_inbound, q = 1:n_inbound] <= 1, Bin);  # alpha => InbPrecedence
@variable(CDSP, 0 <= alphaij[i = 1:n_inbound, j = 1:n_outbound,  i1 = 1:n_inbound] <= 1, Bin); # Foklift Usage => Unloading fisrt, loading later
@variable(CDSP, 0 <= betaij[i = 1:n_inbound, j = 1:n_outbound,  i1 = 1:n_inbound] <= 1, Bin); # Forklift Usage => Loading First, unloading later
@variable(CDSP, 0 <= deltajj[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound,  j1 = 1:n_outbound] <= 1, Bin); # Foklift Usage => loading fisrt, loading later

@variable(CDSP, 0 <= St[i = 1:n_inbound] ); # St_if => Forklift Unloading Start
#@variable(CDSP, 0 <= Lt[i = 1:n_inbound, f = 1:n_forklifts] <= HorizEnd); # Lt_if => Forklift Unloading Finish
@variable(CDSP, 0 <= Sl[i = 1:n_inbound, j = 1:n_outbound] ); # Sl_ijf => Forklift loading Starts
@variable(CDSP, 0 <= Ll[i = 1:n_inbound, j = 1:n_outbound] ); # Sl_ijf => Forklift loading Starts



@constraint(CDSP, LinearizeCon1[j = 1:n_outbound], MaxPwl[j] >= Multiplier[1]*Penalty[j]*Tardiness[j]);
@constraint(CDSP, LinearizeCon2[j = 1:n_outbound], MaxPwl[j] >= Multiplier[2]*Penalty[j]*Tardiness[j] + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]));
@constraint(CDSP, LinearizeCon3[j = 1:n_outbound], MaxPwl[j] >= Multiplier[3]*Penalty[j]*Tardiness[j] +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]));
@constraint(CDSP, InbAssigCon[i = 1:n_inbound], sum(InbAssign[i, k] for k = 1:n_strip) == 1);
@constraint(CDSP, InbAssigFork[i = 1:n_inbound], sum(UnFork[i, f] for f = 1:n_forklifts) == 1); # Unloading Forklift assignment
@constraint(CDSP, OutAssigCon[j = 1:n_outbound], sum(OutAssign[j, g] for g = 1:n_stack) == 1);
@constraint(CDSP, EnterStripCon[i = 1:n_inbound], EntryInbound[i] >= Arrival[i]);
@constraint(CDSP, EnterUnFork[i = 1:n_inbound], St[i] >= EntryInbound[i]); # Starting time of forklift unloading

#@constraint(CDSP, LeaveStripCon[i = 1:n_inbound, k = 1:n_strip], LeaveInbound[i, k] - EntryInbound[i, k] >= InbAssign[i, k]*sum(Supply[i, p]*LoadTime[p] for p = 1:n_products));
#@constraint(CDSP, LeaveUnFork[i = 1:n_inbound, f = 1:n_forklifts], Lt[i, f] - St[i, f] >= UnFork[i, f]*sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products)); # End Time for forklift unloading

#@constraint(CDSP, LeaveStripCon1[i = 1:n_inbound], sum(LeaveInbound[i, k] for k = 1:n_strip) >= sum(Lt[i,f] for f = 1:n_forklifts));

@constraint(CDSP, LeaveStripCon1[i = 1:n_inbound], LeaveInbound[i] >= St[i]+sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products));

@constraint(CDSP, InbFirstCon1[i = 1:n_inbound, q = 1:n_inbound, k = 1:n_strip; i != q], FirstInAssign[q] + FirstInAssign[i] <= 3 - InbAssign[i,k] - InbAssign[q,k]);
@constraint(CDSP, InbLastCon1[i = 1:n_inbound, q = 1:n_inbound, k = 1:n_strip; i != q], LastInAssign[q] + LastInAssign[i] <= 3 - InbAssign[i,k] - InbAssign[q,k]);
@constraint(CDSP, OutFirstCon1[i = 1:n_outbound, q = 1:n_outbound, k = 1:n_stack; i != q], FirstOutAssign[q] + FirstOutAssign[i] <= 3 - OutAssign[i,k] - OutAssign[q,k]);
@constraint(CDSP, OutLastCon1[i = 1:n_outbound, q = 1:n_outbound, k = 1:n_stack; i != q], LastOutAssign[q] + LastOutAssign[i] <= 3 - OutAssign[i,k] - OutAssign[q,k]);


# pick a sequence according to first and last trucks
@constraint(CDSP, InbFirstCon2[q = 1:n_inbound], FirstInAssign[q] + sum(InbPrecedence[i,q] for i = 1:n_inbound if  i != q) == 1);
@constraint(CDSP, InbFirstCon3[i = 1:n_inbound], LastInAssign[i] + sum(InbPrecedence[i,q] for q = 1:n_inbound if  q != i) == 1);
@constraint(CDSP, OutFirstCon2[q = 1:n_outbound], FirstOutAssign[q] + sum(OutPrecedence[i,q] for i = 1:n_outbound if  i != q) == 1);
@constraint(CDSP, OutFirstCon3[i = 1:n_outbound], LastOutAssign[i] + sum(OutPrecedence[i,q] for q = 1:n_outbound if  q != i) == 1);


@constraint(CDSP, InbPrecCon1[i = 1:n_inbound, q = 1:n_inbound; i != q], InbPrecedence[i, q] => {EntryInbound[q] >= Change_time + LeaveInbound[i]});                                                                
#@constraint(CDSP, InbPrecCon1[i = 1:n_inbound, q = 1:n_inbound; i != q], EntryInbound[q] >= Change_time + LeaveInbound[i] - HorizEnd*(1 - InbPrecedence[i, q]));

@constraint(CDSP, InbPrecCon2[i = 1:n_inbound, q = 1:n_inbound, k = 1:n_strip; i != q], 1 - InbPrecedence[i, q] >= InbAssign[i, k] - InbAssign[q, k])
@constraint(CDSP, InbPrecCon3[i = 1:n_inbound, q = 1:n_inbound, k = 1:n_strip; i != q], InbPrecedence[i, q] - 1 <= InbAssign[i, k] - InbAssign[q, k]);
@constraint(CDSP, InbPrecCon4[i = 1:n_inbound], InbPrecedence[i, i] == 0);

@constraint(CDSP, OutPrecCon1[j = 1:n_outbound, r = 1:n_outbound; j != r], OutPrecedence[j,r] => {EntryOutbound[r] >= Change_time + LeaveOutbound[j]});
#@constraint(CDSP, OutPrecCon1[j = 1:n_outbound, r = 1:n_outbound; j != r], EntryOutbound[r] >= Change_time + LeaveOutbound[j] - HorizEnd*(1 - OutPrecedence[j,r]));
                                                                
@constraint(CDSP, OutPrecCon2[i = 1:n_outbound, q = 1:n_outbound, k = 1:n_stack; i != q], 1 - OutPrecedence[i, q]  >= OutAssign[i, k] - OutAssign[q, k]);
@constraint(CDSP, OutPrecCon3[i = 1:n_outbound, q = 1:n_outbound, k = 1:n_stack; i != q], OutPrecedence[i, q] - 1  <= OutAssign[i, k] - OutAssign[q, k]);
                                                                
@constraint(CDSP, OutPrecCon4[i = 1:n_outbound], OutPrecedence[i, i] == 0);

@constraint(CDSP, TransCons1[i = 1:n_inbound, j = 1:n_outbound], sum(LoadTrans[i, j, p] for p = 1:n_products) <= min(sum(Supply[i,p] for p = 1:n_products),Demand_max[j])*Exch[i,j]);

@constraint(CDSP, SupplyCon[i = 1:n_inbound, p = 1:n_products], sum(LoadTrans[i, j, p] for j = 1:n_outbound) <= Supply[i, p]);
@constraint(CDSP, DemandCon[j = 1:n_outbound, p = 1:n_products], sum(LoadTrans[i, j, p] for i = 1:n_inbound) >= Demand[j, p]);
                                                                
@constraint(CDSP, LoadFirst[i = 1:n_inbound, j = 1:n_outbound], FirstLoad[i,j] + sum(alphaij[i, j, q] for q = 1:n_inbound) + sum(deltajj[i1,j1,i,j] for i1 = 1:n_inbound for j1 = 1:n_outbound if i1 != i || j1 != j) == 1);     
                                                                                                
@constraint(CDSP, LoadLast[i = 1:n_inbound, j = 1:n_outbound], LastLoad[i,j] + sum(betaij[i, j, q] for q = 1:n_inbound) + sum(deltajj[i,j,i1,j1] for i1 = 1:n_inbound for j1 = 1:n_outbound if i1 != i || j1 != j) == 1);          
                                                                                                                                
@constraint(CDSP, UnLoadFirst[i = 1:n_inbound], FirstUnload[i] + sum(betaij[i1, j, i] for i1 = 1:n_inbound for j = 1:n_outbound) + sum(gamma[i1,i] for i1 = 1:n_inbound if i1 != i) == 1);     
                                                                                                
@constraint(CDSP, UnLoadLast[i = 1:n_inbound], LastUnload[i] + sum(alphaij[i1, j, i] for i1 = 1:n_inbound for j = 1:n_outbound) + sum(gamma[i,i1] for i1 = 1:n_inbound if i1 != i) == 1);                                                                                                                                          
                                                                                                                                
                                                                                                                                
@constraint(CDSP, First1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound, f = 1:n_forklifts], FirstLoad[i, j] + FirstUnload[q] <= 3 - LFork[i,j,f] - UnFork[q,f]);

@constraint(CDSP, First2[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound,  f = 1:n_forklifts; i1 != i || j1 != j], FirstLoad[i, j] + FirstLoad[i1,j1] <= 3 - LFork[i,j,f] - LFork[i1,j1,f]);  
                                                                                                                                
@constraint(CDSP, First3[i = 1:n_inbound, q = 1:n_inbound, f = 1:n_forklifts; q != i], FirstUnload[i] + FirstUnload[q] <= 3 - UnFork[i,f] - UnFork[q,f]);                                                                                                                                
@constraint(CDSP, Last1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound, f = 1:n_forklifts], LastLoad[i, j] + LastUnload[q] <= 3 - LFork[i,j,f] - UnFork[q,f]);

@constraint(CDSP, Last2[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound,  f = 1:n_forklifts; i1 != i || j1 != j], LastLoad[i, j] + LastLoad[i1,j1] <= 3 - LFork[i,j,f] - LFork[i1,j1,f]);  
                                                                                                                                
@constraint(CDSP, Last3[i = 1:n_inbound, q = 1:n_inbound, f = 1:n_forklifts; q != i], LastUnload[i] + LastUnload[q] <= 3 - UnFork[i,f] - UnFork[q,f]);                                                                                                                                                                                                                                                                

@constraint(CDSP, UnLoadLoad1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound], alphaij[i,j,q] => {Sl[i,j] >= C_time + St[q] + sum(Supply[q, p]*UnLoadTime[p] for p = 1:n_products)}); # Loading followed by unloading
#@constraint(CDSP, UnLoadLoad1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound], Sl[i,j] >= C_time + St[q] + sum(Supply[q, p]*UnLoadTime[p] for p = 1:n_products)- HorizEnd*(1 - alphaij[i, j, q])); # Loading followed by unloading

@constraint(CDSP, UnLoadLoad1a[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound], betaij[i,j,q] => {St[q] >= C_time + Ll[i,j]}); # Loading followed by unloading
#@constraint(CDSP, UnLoadLoad1a[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound], St[q] >= C_time + Ll[i,j]- HorizEnd*(1 - betaij[i, j, q])); # Loading followed by unloading

@constraint(CDSP, UnLoadLoad2[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound,f = 1:n_forklifts], alphaij[i,j,q] - 1 <=  LFork[i, j, f] - UnFork[q,f]); # Loading followed by unloading
                                                                                                                                
@constraint(CDSP, UnLoadLoad3[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound,f = 1:n_forklifts], 1 - alphaij[i,j,q] >=  LFork[i, j, f] - UnFork[q,f]); # Loading followed by unloading                                                                                             
                                                                                                                                
@constraint(CDSP, UnLoadLoad4[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound,f = 1:n_forklifts], betaij[i,j,q] - 1 <=  LFork[i, j, f] - UnFork[q,f]); # Loading followed by unloading

@constraint(CDSP, UnLoadLoad5[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound,f = 1:n_forklifts], 1 - betaij[i,j,q] >=  LFork[i, j, f] - UnFork[q,f]); # Loading followed by unloading                                                                                               
                                                                                                                                

@constraint(CDSP, LoadLoad1[i = 1:n_inbound, i1 = 1:n_inbound, j = 1:n_outbound, j1 = 1:n_outbound, f = 1:n_forklifts; i!=i1 || j !=j1], deltajj[i,j,i1,j1] => {Sl[i1,j1] >= C_time + Ll[i, j]}); 
#@constraint(CDSP, LoadLoad1[i = 1:n_inbound, i1 = 1:n_inbound, j = 1:n_outbound, j1 = 1:n_outbound, f = 1:n_forklifts; i!=i1 || j !=j1], Sl[i1,j1] >= C_time + Ll[i, j] - HorizEnd*(1 - deltajj[i, j, i1, j1])); # Loading followed by loading

@constraint(CDSP, LoadLoad2[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound, f = 1:n_forklifts; i!=i1 || j !=j1], deltajj[i,j,i1,j1] - 1 <=  LFork[i, j, f] - LFork[i1,j1,f]); # Loading followed by unloading
                                                                                                                                
@constraint(CDSP, LoadLoad3[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound, f = 1:n_forklifts; i!=i1 || j !=j1], 1 - deltajj[i,j,i1,j1] >=  LFork[i, j, f] - LFork[i1,j1,f]); # Loading followed by unloading                                                                                             

@constraint(CDSP, LoadLoad4[i = 1:n_inbound, j = 1:n_outbound], deltajj[i,j,i,j] == 0); # UnLoading followed by Loading

@constraint(CDSP, AuxUnLoadUnLoad1[i = 1:n_inbound, q = 1:n_inbound; i != q], 1 - InbPrecedence[i,q] >= gamma[q,i]); # Unloading followed by unloading
@constraint(CDSP, AuxUnLoadLoad1[i = 1:n_inbound, j = 1:n_outbound, q = 1:n_inbound; i != q], 1 - InbPrecedence[i,q] >= alphaij[i,j,q]);
@constraint(CDSP, AuxUnLoadLoad2[i = 1:n_inbound, j = 1:n_outbound], betaij[i,j,i] == 0);
@constraint(CDSP, AuxLoadLoad1[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound; i!=i1], deltajj[i,j,i1,j1]  <=  1 - InbPrecedence[i1,i]);# Loading followed by unloading
@constraint(CDSP, AuxLoadLoad2[i = 1:n_inbound, j = 1:n_outbound, i1 = 1:n_inbound, j1 = 1:n_outbound; j!=j1], deltajj[i,j,i1,j1]  <=  1 - OutPrecedence[j1,j]);# Loading followed by unloading

@constraint(CDSP, UnLoadUnLoad1[i = 1:n_inbound, q = 1:n_inbound, f = 1:n_forklifts; i != q], gamma[i,q] => {St[q] >= C_time + St[i] + sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products)}); # Unloading followed by unloading
#@constraint(CDSP, UnLoadUnLoad1[i = 1:n_inbound, q = 1:n_inbound, f = 1:n_forklifts; i != q], St[q] >= C_time + St[i] + sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products) - HorizEnd*(1 - gamma[i, q])); # Unloading followed by unloading

@constraint(CDSP, UnLoadUnLoad2[i = 1:n_inbound, q = 1:n_inbound,f = 1:n_forklifts; q != i], gamma[i,q] - 1 <=  UnFork[i, f] - UnFork[q,f]); # Loading followed by unloading
                                                                                                                                
@constraint(CDSP, UnLoadUnLoad3[i = 1:n_inbound, q = 1:n_inbound,f = 1:n_forklifts; q != i], 1 - gamma[i,q] >=  UnFork[i, f] - UnFork[q,f]); # Loading followed by unloading        

@constraint(CDSP, UnLoadUnLoad4[i = 1:n_inbound], gamma[i,i] == 0); # UnLoading followed by unloading

@constraint(CDSP, LoadAssign[i = 1:n_inbound, j = 1:n_outbound], sum(LFork[i, j, f] for f = 1:n_forklifts) == Exch[i,j]); # Assign one loading task to one forklift, otherwise divide loads of 'rijp' between multiple forklifts, making the problem complicated  

@constraint(CDSP, StartLoadFork[i = 1:n_inbound, j = 1:n_outbound], Exch[i,j] => {Sl[i, j] >= St[i] + sum(Supply[i,p]*UnLoadTime[p] for p in 1:n_products)/2}); # Starting time of forklift unloading

@constraint(CDSP, EndLoadFork[i = 1:n_inbound, j = 1:n_outbound], Ll[i, j] >= Sl[i, j] + sum(LoadTrans[i,j,p]*LoadTime[p] for p in 1:n_products)); #End time of forklift unloading

@constraint(CDSP, LeaveStackCon1[i = 1:n_inbound, j = 1:n_outbound], LeaveOutbound[j] >= Ll[i, j]);
# ensures minimum stay of outbound truck;
# truck should be docked before starting any loading job; truck latest docking time given by precedence constraints;
@constraint(CDSP, LeaveStackCon2[j = 1:n_outbound], LeaveOutbound[j] >= EntryOutbound[j] + sum(Demand[j, p]*LoadTime[p] for p = 1:n_products));

@constraint(CDSP, LeaveStackCon3[i = 1:n_inbound, j = 1:n_outbound], Sl[i, j] >= EntryOutbound[j]);


@constraint(CDSP, DeadCon[j = 1:n_outbound], LeaveOutbound[j] <= DeadLines[j] + Tardiness[j]);
#@constraint(CDSP, InbLastEntry[i = 1:n_inbound, k = 1:n_strip], EntryInbound[i, k] <= InbAssign[i,k]*HorizEnd);
#@constraint(CDSP, MakespanCon[j = 1:n_outbound, g = 1:n_stack], Makespan[j] >= LeaveOutbound[j,g]);

@objective(CDSP, Min, sum(MaxPwl[j] for j = 1:n_outbound));
