#updated objective to makespan

using JuMP, Gurobi;

CDSP = Model(optimizer_with_attributes(Gurobi.Optimizer, "TIME_LIMIT" => 1.95*3600));

@variable(CDSP, 0 <= InbAssign[i = 1:n_inbound, k = 1:n_doors] <= 1, Bin); # InbAssign => x
@variable(CDSP, 0 <= OutAssign[j = 1:n_outbound, g = 1:n_doors] <= 1, Bin); # OutAssign => y
@variable(CDSP, 0 <= FirstInAssign[j = 1:n_inbound] <= 1, Bin ); # OutAssign => y
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

@variable(CDSP, 0 <= FirstInEntry[i = 1:n_inbound] <= HorizEnd); # e_I => EntryInbound
@variable(CDSP, 0 <= LastInLeave[i = 1:n_inbound] <= HorizEnd); # l_I => LeaveInbound
@variable(CDSP, 0 <= FirstOutEntry[j = 1:n_outbound] <= HorizEnd); # e_O => EntryOutbound
@variable(CDSP, 0 <= LastOutLeave[j = 1:n_outbound] <= HorizEnd); # l_O => LeaveOutbound

#@variable(CDSP, 0 <= ExchLeaveOutbound[j = 1:n_outbound, g = 1:n_inbound] <= HorizEnd); # l_O => LeaveOutbound
@variable(CDSP, 0 <= Tardiness[j = 1:n_outbound] <= HorizEnd); # gamma => Tardiness
@variable(CDSP, 0 <= MaxPwl[j = 1:n_outbound]); # z => MaxPwl
@variable(CDSP, 0 <= Makespan[j = 1:n_outbound]);



@variable(CDSP, MixPrecedence[i = 1:n_inbound, r = 1:n_outbound], Bin); # gamma => MixPrecedence
@variable(CDSP, MixPrecedence1[i = 1:n_outbound, r = 1:n_inbound], Bin); # gamma => MixPrecedence


@constraint(CDSP, LinearizeCon1[j = 1:n_outbound], MaxPwl[j] >= Multiplier[1]*Penalty[j]*Tardiness[j]);
@constraint(CDSP, LinearizeCon2[j = 1:n_outbound], MaxPwl[j] >= Multiplier[2]*Penalty[j]*Tardiness[j] + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]));
@constraint(CDSP, LinearizeCon3[j = 1:n_outbound], MaxPwl[j] >= Multiplier[3]*Penalty[j]*Tardiness[j] +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]));
@constraint(CDSP, InbAssigCon[i = 1:n_inbound], sum(InbAssign[i, k] for k = 1:n_doors) == 1);
@constraint(CDSP, OutAssigCon[j = 1:n_outbound], sum(OutAssign[j, g] for g = 1:n_doors) == 1);

@constraint(CDSP, EnterStripCon[i = 1:n_inbound], EntryInbound[i] >= Arrival[i]);
@constraint(CDSP, LeaveStripCon[i = 1:n_inbound], LeaveInbound[i] - EntryInbound[i] >= sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products));



@constraint(CDSP, InbFirstCon1[i = 1:n_inbound, q = 1:n_inbound, k = 1:n_doors; i != q], FirstInAssign[q] + FirstInAssign[i] <= 3 - InbAssign[i,k] - InbAssign[q,k]);

@constraint(CDSP, InbLastCon1[i = 1:n_inbound, q = 1:n_inbound, k = 1:n_doors; i != q], LastInAssign[q] + LastInAssign[i] <= 3 - InbAssign[i,k] - InbAssign[q,k]);
@constraint(CDSP, OutFirstCon1[i = 1:n_outbound, q = 1:n_outbound, k = 1:n_doors; i != q], FirstOutAssign[q] + FirstOutAssign[i] <= 3 - OutAssign[i,k] - OutAssign[q,k]);
@constraint(CDSP, OutLastCon1[i = 1:n_outbound, q = 1:n_outbound, k = 1:n_doors; i != q], LastOutAssign[q] + LastOutAssign[i] <= 3 - OutAssign[i,k] - OutAssign[q,k]);

@constraint(CDSP, MixFirstCon1[i = 1:n_inbound, q = 1:n_outbound, k = 1:n_doors], FirstInAssign[i] + FirstOutAssign[q] <= 3 - InbAssign[i,k] - OutAssign[q,k]);
@constraint(CDSP, MixLastCon1[i = 1:n_inbound, q = 1:n_outbound, k = 1:n_doors], LastInAssign[i] + LastOutAssign[q] <= 3 - InbAssign[i,k] - OutAssign[q,k]);


# constraints include sequencing between inbound and outbound trucks
@constraint(CDSP, InbFirstCon2[q = 1:n_inbound], FirstInAssign[q] + sum(InbPrecedence[i,q] for i = 1:n_inbound if  i != q) + sum(MixPrecedence1[i,q] for i = 1:n_outbound) == 1);
@constraint(CDSP, InbFirstCon3[i = 1:n_inbound], LastInAssign[i] + sum(InbPrecedence[i,q] for q = 1:n_inbound if  q != i)  + sum(MixPrecedence[i,q] for q = 1:n_outbound) == 1);
@constraint(CDSP, OutFirstCon2[q = 1:n_outbound], FirstOutAssign[q] + sum(OutPrecedence[i,q] for i = 1:n_outbound if  i != q) + sum(MixPrecedence[i,q] for i = 1:n_inbound) == 1);
@constraint(CDSP, OutFirstCon3[i = 1:n_outbound], LastOutAssign[i] + sum(OutPrecedence[i,q] for q = 1:n_outbound if  q != i)  + sum(MixPrecedence1[i,q] for q = 1:n_inbound) == 1);

@constraint(CDSP, InbPrecCon1[i = 1:n_inbound, q = 1:n_inbound; i != q], EntryInbound[q] >= Change_time + LeaveInbound[i] - M1*(1 - InbPrecedence[i, q]));

@constraint(CDSP, InbPrecCon2[i = 1:n_inbound, q = 1:n_inbound, k = 1:n_doors; i != q], 1 - InbPrecedence[i, q] >= InbAssign[i, k] - InbAssign[q, k])
@constraint(CDSP, InbPrecCon3[i = 1:n_inbound, q = 1:n_inbound, k = 1:n_doors; i != q], InbPrecedence[i, q] - 1 <= InbAssign[i, k] - InbAssign[q, k]);
@constraint(CDSP, InbPrecCon4[i = 1:n_inbound], InbPrecedence[i, i] == 0);


@constraint(CDSP, OutPrecCon1[j = 1:n_outbound, r = 1:n_outbound; j != r], EntryOutbound[r] >= Change_time + LeaveOutbound[j] - M2*(1 - OutPrecedence[j,r]));
                                                                
@constraint(CDSP, OutPrecCon2[i = 1:n_outbound, q = 1:n_outbound, k = 1:n_doors; i != q], 1 - OutPrecedence[i, q]  >= OutAssign[i, k] - OutAssign[q, k]);
@constraint(CDSP, OutPrecCon3[i = 1:n_outbound, q = 1:n_outbound, k = 1:n_doors; i != q], OutPrecedence[i, q] - 1  <= OutAssign[i, k] - OutAssign[q, k]);
                                                                
@constraint(CDSP, OutPrecCon4[i = 1:n_outbound], OutPrecedence[i, i] == 0);


@constraint(CDSP, MixPrecCon1[i = 1:n_inbound, q = 1:n_outbound], EntryOutbound[q] >= Change_time-M1*(1-MixPrecedence[i,q]) + LeaveInbound[i]);
@constraint(CDSP, MixPrecCon1a[i = 1:n_outbound, q = 1:n_inbound], EntryInbound[q] >= Change_time-M2*(1-MixPrecedence1[i,q]) + LeaveOutbound[i]);
                                
@constraint(CDSP, MixPrecCon2[i = 1:n_inbound, q = 1:n_outbound, k = 1:n_doors], 1 - MixPrecedence[i, q] >= InbAssign[i, k] - OutAssign[q, k])
@constraint(CDSP, MixPrecCon2a[i = 1:n_outbound, q = 1:n_inbound, k = 1:n_doors], 1 - MixPrecedence1[i, q] >= OutAssign[i, k] - InbAssign[q, k])                                                                                                
                                                                                                
@constraint(CDSP, MixPrecCon3[i = 1:n_inbound, q = 1:n_outbound, k = 1:n_doors], MixPrecedence[i, q] - 1 <= InbAssign[i, k] - OutAssign[q, k]);
@constraint(CDSP, MixPrecCon3a[i = 1:n_outbound, q = 1:n_inbound, k = 1:n_doors], MixPrecedence1[i, q] - 1 <= OutAssign[i, k] - InbAssign[q, k]) 
                                                                                                
                                                                                                
@constraint(CDSP, TransCons1[i = 1:n_inbound, j = 1:n_outbound], sum(LoadTrans[i, j, p] for p = 1:n_products) <= Demand_max[j]*Exch[i,j]);
#@constraint(CDSP, TransCons2[i = 1:n_inbound, j = 1:n_outbound], sum(TransRoute[i, j, k, g] for k = 1:n_doors, g = 1:n_doors) <= 1);
#@constraint(CDSP, TransCons3[i = 1:n_inbound, j = 1:n_outbound, k = 1:n_doors, g = 1:n_doors], TransRoute[i, j, k, g] <= InbAssign[i, k]);
#@constraint(CDSP, TransCons4[i = 1:n_inbound, j = 1:n_outbound, k = 1:n_doors, g = 1:n_doors], TransRoute[i, j, k, g] <= OutAssign[j, g]);
@constraint(CDSP, SupplyCon[i = 1:n_inbound, p = 1:n_products], sum(LoadTrans[i, j, p] for j = 1:n_outbound) <= Supply[i, p]);
@constraint(CDSP, DemandCon[j = 1:n_outbound, p = 1:n_products], sum(LoadTrans[i, j, p] for i = 1:n_inbound) >= Demand[j, p]);
                                                                                                                                
@constraint(CDSP, LeaveStackCon1[i = 1:n_inbound, j = 1:n_outbound], LeaveOutbound[j] >= EntryInbound[i] + sum(Supply[i, p]*UnLoadTime[p] for p = 1:n_products)/2 + sum(LoadTrans[i, j, p]*LoadTime[p] for p = 1:n_products) - M1*(1-Exch[i, j]));
@constraint(CDSP, LeaveStackCon2[j = 1:n_outbound], LeaveOutbound[j] >= EntryOutbound[j] + sum(Demand[j, p]*LoadTime[p] for p = 1:n_products));
#@constraint(CDSP, LeaveStackCon3[j = 1:n_outbound, g = 1:n_outbound; j!=g], EntryOutbound[j] >= PrecEntryOutbound[g,j]);
#@constraint(CDSP, LeaveStackCon4[j = 1:n_outbound, g = 1:n_outbound; j!=g], LeaveOutbound[j] <= PrecLeaveOutbound[j,g]);

@constraint(CDSP, DeadCon[j = 1:n_outbound], LeaveOutbound[j] <= DeadLines[j] + Tardiness[j]);                                                                                                                                

@objective(CDSP, Min, sum(MaxPwl[j] for j = 1:n_outbound));
