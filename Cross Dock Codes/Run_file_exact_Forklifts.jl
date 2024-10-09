include("data_v1.2.jl");

using CSV;

resultfile = open("Supply_Demand_exact.csv", "w");
cols = vcat(["Trailer_num", ", Trailer_type", ", Supply/Demand/Loadtime"], [", P$(p)" for p in 1:n_products]);
cols1 = join(cols)
println(resultfile, cols1);
ltime = vcat(["NA", ", NA", ", Loadtime"], [", $(LoadTime[p])" for p in 1:n_products]);
ltime1 = join(ltime);
println(resultfile, ltime1);
for i in 1:n_inbound
    entry = vcat(["IT$(i)", ", Inbound trailer", ", Supply"], [", $(Supply[i, p])" for p in 1:n_products]);
    entry1 = join(entry)
    println(resultfile, entry1);
end;
for j in 1:n_outbound
    entry2 = vcat(["OT$(j)", ", Outbound trailer", ", Demand"], [", $(Demand[j, p])" for p in 1:n_products]);
    entry3 = join(entry2)
    println(resultfile, entry3);
end;
close(resultfile);

n_strip = Int(round(mean(UnLoadTime)/(mean(UnLoadTime)+mean(LoadTime))*n_doors));;
n_stack = n_doors-n_strip;
#=
Entr_IT = zeros(Float64,n_inbound,n_strip);
Leav_IT = zeros(Float64,n_inbound,n_strip);
Entr_OT = zeros(Float64,n_outbound,n_stack);
Leav_OT = zeros(Float64,n_outbound,n_stack);
Load_Tran = zeros(Int32,n_inbound,n_outbound,n_products);
Out_Entr = Dict();
Out_Leav = Dict();
Inf_Entr = zeros(Float64,n_inbound,n_forklifts); # Latest times available for forklifts
Inf_Leav = zeros(Float64,n_inbound,n_forklifts);
Outf_Entr = zeros(Float64,n_outbound,n_forklifts);
Outf_Leav = zeros(Float64,n_outbound,n_forklifts);
for i in 1:n_outbound
    for j in 1:n_forklifts
        Out_Entr[(i,j)] = Array{Float64}(undef,0);
        Out_Leav[(i,j)] = Array{Float64}(undef,0);
    end;
end;
=#

println("Method Used - Forklifts - Integrated");
println("Inbound - $(n_inbound), Outbound - $(n_outbound), Doors - $(n_strip)_$(n_stack), Forklifts - $(n_forklifts)");
include("code_v1.8_Forklift_Modified.jl");

sol_time = @elapsed begin
optimize!(CDSP);
end;

resultfile = open("Solution_file_exact_Forklifts_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "w");
println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_strip, n_stack, n_inbound, n_outbound, n_products, Change_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, Assigned_DD, DD_entry_time, DD_leave_time, Tardiness, Objective value, Termination Status , Solution_time");
for i in 1:n_inbound
    door_assign = sum(k for k = 1:n_strip if value(InbAssign[i, k]) > 0);
    println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), Strip Door $(door_assign), $(value(EntryInbound[i])), $(value(LeaveInbound[i])), NA, $(objective_value(CDSP)), $(termination_status(CDSP)), $(sol_time)");
end;
            
for j in 1:n_outbound
        door_assign = sum(g for g = 1:n_stack if value(OutAssign[j, g]) > 0);
        println(resultfile, "OT$(j), Outbound trailer, $(0), $(DeadLines[j]), $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), NA, $(sum(Demand[j,p] for p = 1:n_products)), $(Penalty[j]), $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), Stack Door $(door_assign), $(value(EntryOutbound[j])), $(value(LeaveOutbound[j])), $(value(Tardiness[j])), $(objective_value(CDSP)), $(termination_status(CDSP)), $(sol_time)");
end;
close(resultfile);

df = DataFrame(
Door= cat(["Strip Door $k" for i in 1:n_inbound, k = 1:n_strip if value(InbAssign[i,k]) == 1], ["Stack Door $g" for j in 1:n_outbound, g = 1:n_stack if value(OutAssign[j,g]) == 1], dims = 1),
Start_time = cat([value(EntryInbound[i]) for i in 1:n_inbound, k = 1:n_strip if value(InbAssign[i,k]) == 1], [value(EntryOutbound[j]) for j =1:n_outbound, g = 1:n_stack if value(OutAssign[j,g]) == 1], dims = 1) ,
Stop_time = cat([value(LeaveInbound[i]) for i in 1:n_inbound, k = 1:n_strip if value(InbAssign[i,k]) == 1], [value(LeaveOutbound[j]) for j =1:n_outbound, g = 1:n_stack if value(OutAssign[j,g]) == 1], dims = 1),
Trailer_type = cat(["Inbound Trailer" for i in 1:n_inbound, k = 1:n_strip if value(InbAssign[i,k]) == 1], ["Outbound Trailer" for j =1:n_outbound, g = 1:n_stack if value(OutAssign[j,g]) == 1], dims = 1),
Labels = cat(["IT-$i" for i in 1:n_inbound, k = 1:n_strip if value(InbAssign[i,k]) == 1], ["OT-$j" for j =1:n_outbound, g = 1:n_stack if value(OutAssign[j,g]) == 1], dims = 1),
Label_pos = (cat([value(LeaveInbound[i]) for i in 1:n_inbound, k = 1:n_strip if value(InbAssign[i,k]) == 1], [value(LeaveOutbound[j]) for j =1:n_outbound, g = 1:n_stack if value(OutAssign[j,g]) == 1], dims = 1)).+ 2
);                        
                        
CSV.write("forklift_output.csv", df, header=true)

resultfile = open("Solution_file_Forklift_Swap_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "w");
println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_strip, n_stack, n_inbound, n_outbound, n_forklifts, n_products, C_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, Assigned_forklift, forklift_start_time, forklift_leave_time, Tardiness, Objective value");
for i in 1:n_inbound
    forklift_assign = sum(k for k = 1:n_forklifts if value(UnFork[i, k] > 0.1);
    println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_forklifts), $(n_products), $(C_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(forklift_assign), $(value(St[i,forklift_assign])), $(value(St[i,forklift_assign])+sum(Supply[i,p] for p = 1:n_products)), NA, $(objective_value(CDSP)) ");
end;

for i in 1:n_outbound
    forklift_assign = sum(k for k = 1:n_forklifts if value(LFork[i, k] > 0.1);
    println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_forklifts), $(n_products), $(C_time), $(sum(Demand[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(forklift_assign), $(value(Sl[i,forklift_assign])), $(value(Ll[i,forklift_assign])), NA, $(objective_value(CDSP)) ");
end;
close(resultfile);

#=
for i in 1:n_inbound
    for j in 1:n_strip
        if value(EntryInbound[i,j]) > 0
            Entr_IT[i,j] = value(EntryInbound[i,j]);
            Leav_IT[i,j] = value(LeaveInbound[i,j]);
        end;
    end;
end;

for i in 1:n_outbound
    for j in 1:n_stack
        if value(OutAssign[i,j]) == 1
            Entr_OT[i,j] = value(EntryOutbound[i,j]);
            Leav_OT[i,j] = value(LeaveOutbound[i,j]);
        end;
    end;
end;

for j in 1:n_outbound
    for i in 1:n_inbound
        for k in 1:n_forklifts
            if value(LFork[i,j,k]) == 1
                push!(Out_Entr[(j,k)],value(Sl[i,j,k]));
                push!(Out_Leav[(j,k)],value(Ll[i,j,k]));
            end;
        end;
    end;
end;


for i in 1:n_inbound
    for j in 1:n_outbound
        if sum([value(LoadTrans[i,j,p]) for p in 1:n_products]) > 0
            for p in 1:n_products
                Load_Tran[i,j,p] = value(LoadTrans[i,j,p])
            end;
        end;
    end;
end;
=#
#include("schedule_gantt_exact_Forklifts.jl");
#=
if (termination_status(CDSP) == MOI.OPTIMAL) || ((termination_status(CDSP) == MOI.TIME_LIMIT) && has_values(CDSP))

    include("schedule_gantt_exact_Forklifts.jl");

    resultfile = open("Solution_file_exact_Forklifts.csv", "w");
    println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_strip, n_stack, n_inbound, n_outbound, n_products, Change_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, Assigned_DD, DD_entry_time, DD_leave_time, Tardiness, Objective value, Termination Status , Solution_time");
    for i in 1:n_inbound
        door_assign = sum(k for k = 1:n_strip if value(InbAssign[i, k]) == 1);
        println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), Strip Door $(door_assign), $(value(EntryInbound[i,door_assign])), $(value(LeaveInbound[i,door_assign])), NA, $(objective_value(CDSP)), $(termination_status(CDSP)), $(sol_time)");
    end;
    for j in 1:n_outbound
        door_assign = sum(g for g = 1:n_stack if value(OutAssign[j, g]) == 1);
        println(resultfile, "OT$(j), Outbound trailer, $(0), $(DeadLines[j]), $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), NA, $(sum(Demand[j,p] for p = 1:n_products)), $(Penalty[j]), $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), Stack Door $(door_assign), $(value(EntryOutbound[j, door_assign])), $(value(LeaveOutbound[j, door_assign])), $(value(Tardiness[j])), $(objective_value(CDSP)), $(termination_status(CDSP)), $(sol_time)");
    end;
    close(resultfile);

    resultfile = open("Flows_exact_Forklifts.csv", "w");
    println(resultfile, "Inbound Trailer, Outbound Trailer, Product, Quantity");
    for i in 1:n_inbound
        str_dd = sum(k for k = 1:n_strip if value(InbAssign[i, k]) == 1);
        for j in 1:n_outbound
            sta_dd = sum(g for g = 1:n_stack if value(OutAssign[j, g]) == 1);
            for p in 1:n_products
                if value(LoadTrans[i, j, p]) != 0
                    println(resultfile, "IT$(i), OT$(j), P$(p), $(value(LoadTrans[i, j, p]))");
                end;
            end;
        end;
    end;
    close(resultfile);
    
else

    println("No Feasible solution found in 24 hours time limit")

end;
=#