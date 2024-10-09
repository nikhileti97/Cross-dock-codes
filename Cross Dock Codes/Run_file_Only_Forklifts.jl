# find forklift scheduling from given truck door assignment

include("data_v1.2.jl");

using CSV,DataFrames;

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

println("Method Used - Forklifts - Sequential");
println("Inbound - $(n_inbound), Outbound - $(n_outbound), Doors - $(n_strip)_$(n_stack), Forklifts - $(n_forklifts)");

df_read = CSV.read("Solution_file_exact_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack).csv", DataFrame);

Entr_IT = zeros(Float64,n_inbound);
door_IT = zeros(Int8,n_inbound);
Entry_IT = zeros(Float64,n_inbound,n_strip);
Leave_IT = zeros(Float64,n_inbound,n_strip);
Leav_IT = zeros(Float64,n_inbound);
Entr_OT = zeros(Float64,n_outbound);
door_OT = zeros(Int8,n_outbound);
Leav_OT = zeros(Float64,n_outbound);
Entry_OT = zeros(Float64,n_outbound,n_stack);
Leave_OT = zeros(Float64,n_outbound,n_stack);
Load_Tran = zeros(Int32,n_inbound,n_outbound,n_products);
Start_Load = zeros(Float64,n_inbound,n_outbound);
End_Load = zeros(Float64,n_inbound,n_outbound);
Inb_Precedence = zeros(Int8,n_inbound,n_inbound);
Out_Precedence = zeros(Int8,n_outbound,n_outbound);
Exch = zeros(Int8,n_inbound,n_outbound);
Lateness = zeros(Float64,n_outbound);


for i in 1:n_inbound
    Entr_IT[i] = df_read[i,20];
    expression = df_read[i,19];
    words = split(expression);
    last_element = last(words)

    # Attempt to parse the last element as an integer
    try
        integer_value = parse(Int, last_element)
        door_IT[i] = integer_value
    catch
        println("No integer found in the expression.")
    end
    #door_IT[i] = df_read[i,19];
    Entry_IT[i,door_IT[i]] = Entr_IT[i];
    Leav_IT[i] = df_read[i,20]+sum(Supply[i,p]*UnLoadTime[p] for p in 1:n_products);
    Leave_IT[i,door_IT[i]] = Leav_IT[i];
end;

for i in 1:n_outbound
    Entr_OT[i] = df_read[n_inbound+i,20];
    expression = df_read[n_inbound+i,19];
    words = split(expression);
    last_element = last(words)

    # Attempt to parse the last element as an integer
    try
        integer_value = parse(Int, last_element)
        door_OT[i] = integer_value
    catch
        println("No integer found in the expression.")
    end
    
    #door_OT[i] = df_read[n_inbound+i,19];
    Entry_OT[i,door_OT[i]] = Entr_OT[i];
    Leav_OT[i] = df_read[n_inbound+i,21];
    Leave_OT[i,door_OT[i]] = Leav_OT[i];
    Lateness[i] = parse(Float64,df_read[n_inbound+i,22]);
end;

df_flow = CSV.read("Flows_exact_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack).csv", DataFrame);

for i in 1:length(df_flow[:,1])
    # Input string
    expression = df_flow[i,1];
    expression1 = df_flow[i,2];
    expression2 = df_flow[i,3];
    
    integer_pattern = r"\d+"

# Search for the pattern in the input string
    match_result = match(integer_pattern, expression)
    match_result1 = match(integer_pattern, expression1)
    match_result2 = match(integer_pattern, expression2)

    # Check if a match was found and extract the integer
    if match_result !== nothing
        integer_value = parse(Int, match_result.match)
        integer_value1 = parse(Int, match_result1.match)
        integer_value2 = parse(Int, match_result2.match)
        Load_Tran[integer_value,integer_value1,integer_value2] = df_flow[i,4];
        Exch[integer_value,integer_value1] = 1;
    end;
end;

df_in = CSV.read("In_Precedence_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack).csv", DataFrame);

for i in 1:length(df_in[:,1])
    Inb_Precedence[df_in[i,1],df_in[i,2]] = df_in[i,3];
end;

df_out = CSV.read("Out_Precedence_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack).csv", DataFrame);

for i in 1:length(df_out[:,1])
    Out_Precedence[df_out[i,1],df_out[i,2]] = df_out[i,3];
end;

function Temp_Storage_forklift_plot(Enter_IT,Enter_OT,Leave_IT,Leave_OT,Load_Tran6,Inf_Entr0,Inf_Leav0,Start_Load0,End_Load0)

    inbound_seq = Dict();
    outbound_seq = Dict();
    inbound_outbound = Dict();
    for i in 1:n_strip
        inbound_seq[i] = Array{Int32}(undef,0);
    end;
    for i in 1:n_stack
        outbound_seq[i] = Array{Int32}(undef,0);
    end;

    dummy_Leav_IT = copy(Leave_IT);

    for i in 1:n_strip
        while sum(dummy_Leav_IT[:,i]) > 0
            sel_truck = findmax(dummy_Leav_IT[:,i])[2];
            if dummy_Leav_IT[sel_truck,i] > 0
                push!(inbound_seq[i],sel_truck);
                dummy_Leav_IT[sel_truck,i] = 0;
            end;
        end;
        inbound_seq[i] = reverse(inbound_seq[i]);
    end;



    for i in 1:n_outbound
        sel_door = findmax(Leave_OT[i,:])[2];
        if Leave_OT[i,sel_door] > 0
            push!(outbound_seq[sel_door],i);
        end;
    end;

    for i in 1:n_inbound
        inbound_outbound[i] = Array{Int32}(undef,0);
        for j in 1:n_outbound
            if sum(Load_Tran6[i,j,:]) > 0
                push!(inbound_outbound[i],j);
            end;
        end;
    end;
    
    makespan = Int(ceil(maximum([sum(Leave_OT[j, :]) for j = 1:n_outbound])))+1;
    Storage_Used = zeros(Int64,n_strip,makespan);
    Storage_Add_Used = zeros(Int64,n_strip,makespan);
    Storage_Sub_Used = zeros(Int64,n_strip,makespan);
    Inbound_Entry = [maximum(Enter_IT[i,:]) for i in 1:n_inbound];
    Inbound_Leave = [maximum(Leave_IT[i,:]) for i in 1:n_inbound];
    Outbound_Entry = [maximum(Enter_OT[i,:]) for i in 1:n_outbound];
    Outbound_Leave = [maximum(Leave_OT[i,:]) for i in 1:n_outbound];
    cum_supply = zeros(Int64,n_strip);
    cum_remove = zeros(Int64,n_strip);
    for i in 1:n_strip
        seq_length = length(inbound_seq[i]);
        cnt = 1;
        for IT in inbound_seq[i]
            start_time = Enter_IT[IT,i];
            leave_time = Leave_IT[IT,i];
            t_load = Enter_IT[IT,i];
            for p in 1:n_products
                if Supply[IT,p] > 0
                    for k in 1:Supply[IT,p]
                        t_load += UnLoadTime[p];
                        cum_supply[i] += 1;
                        Storage_Add_Used[i,Int(ceil(t_load))] += 1 ;
                    end;
                end;
            end; 
            for OT in inbound_outbound[IT]
                t_out = Start_Load0[IT,OT];
                t = t_out;
                for p in 1:n_products
                    if Load_Tran6[IT,OT,p] > 0
                        for k in 1:Load_Tran6[IT,OT,p]
                            t += LoadTime[p];
                            cum_remove[i] += 1;
                            Storage_Sub_Used[i,Int(ceil(t))] += 1;
                        end;
                    end;
                end;
            end;
            #=
            if cnt < seq_length
                next_IT = inbound_seq[i][cnt+1];
                t_next = Entr_IT[next_IT,i];
                for t_place in Int(ceil(t_load)):Int(floor(t_next))
                    Storage_Add_Used[i,t_place] = Storage_Add_Used[i,Int(floor(t_load))];
                end;
            else
                for t_place in Int(ceil(t_load)):makespan
                    Storage_Add_Used[i,t_place] = Storage_Add_Used[i,Int(floor(t_load))];
                end;
            end;
            cnt += 1;
            =#
        end;
    end;


    for i in 1:n_strip
        prev_used = 0;
        for t in 1:makespan
            if t != 1
                prev_used = Storage_Used[i,t-1];
            end;
            Storage_Used[i,t] = prev_used+Storage_Add_Used[i,t]-Storage_Sub_Used[i,t];
            
            if Storage_Used[i,t] < 0
                Storage_Used[i,t] = 0;
            end;
            
        end;
    end;
    
    plot([1:makespan], Storage_Used[1,:], title = "Storage for Strip Doors", label = "Strip Door-1");
    
    prev_sum = Storage_Used[1,:];
    for i in 2:n_strip
        cum_sum = prev_sum + Storage_Used[i,:];
        plot!([1:makespan], cum_sum, label = "Strip Door-$(i)")
        prev_sum = cum_sum;
    end;
    println("$(sum(cum_supply[i] for i in 1:n_strip))");
    println("$(sum(cum_remove[i] for i in 1:n_strip))");
    xlabel!("Time Elapsed");
    ylabel!("Temporary Storage in Product units");
    savefig("Temp Storage Utilized-Separate Doors with forklift-Exact");
    
end;


include("code_v1.8_Forklift_Separate.jl");

sol_time = @elapsed begin
optimize!(CDSSP);
end;

Out_Entr = Dict();
Out_Leav = Dict();
Inf_Entr = zeros(Float64,n_inbound,n_forklifts); # Latest times available for forklifts
Inf_Leav = zeros(Float64,n_inbound,n_forklifts);
for i in 1:n_outbound
    for j in 1:n_forklifts
        Out_Entr[(i,j)] = Array{Float64}(undef,0);
        Out_Leav[(i,j)] = Array{Float64}(undef,0);
    end;
end;

for i in 1:n_inbound
    for j in 1:n_outbound
        Start_Load[i,j] = value(Sl[i,j]);
        End_Load[i,j] = value(Ll[i,j]);
    end;
end;



resultfile = open("Solution_file_Forklift_Sequential_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "w");

println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_strip, n_stack, n_inbound, n_outbound, n_forklifts, n_products, C_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, New_Entry_time, New_Leave_time, Assigned_forklift, forklift_start_time, forklift_leave_time, Tardiness, Objective value");

for i in 1:n_inbound
    forklift_assign = sum(k for k = 1:n_forklifts if value(UnFork[i, k]) > 0);
    Inf_Entr[i,forklift_assign] = value(St[i]);
    Inf_Leav[i,forklift_assign] = value(St[i])+sum(Supply[i,p]*UnLoadTime[p] for p = 1:n_products);
    println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_forklifts), $(n_products), $(C_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(value(Entry1Inbound[i])), $(value(Leave1Inbound[i])), $(forklift_assign), $(value(St[i])), $(value(St[i])+sum(Supply[i,p]*UnLoadTime[p] for p = 1:n_products)), NA, $(objective_value(CDSSP)) ");
end;

for i in 1:n_outbound
    for j in 1:n_inbound
        for g in 1:n_forklifts
            if value(LFork[j,i,g]) > 0.1
                push!(Out_Entr[(i,g)],value(Sl[j,i]));
                push!(Out_Leav[(i,g)],value(Ll[j,i]));
                println(resultfile, "OT$(i), Outbound trailer, NA, $(DeadLines[i]), $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_forklifts), $(n_products), $(C_time), NA, $(sum(Demand[i,p] for p = 1:n_products)), $(Penalty[i]), $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(value(Entry1Outbound[i])), $(value(Leave1Outbound[i])), $(g), $(value(Sl[j,i])), $(value(Ll[j,i])), $(value(Tardiness1[j])), $(objective_value(CDSSP))");
            end;
        end;
    end;
end;
close(resultfile);

Temp_Storage_forklift_plot(Entry_IT,Entry_OT,Leave_IT,Leave_OT,Load_Tran,Inf_Entr,Inf_Leav,Start_Load,End_Load);
Job_Storage = Array{Float64}(undef,0);
Job_Delay = Array{Float64}(undef,0);
for i in 1:n_inbound
    for j in 1:n_outbound
        if sum(Load_Tran[i,j,:]) > 0
            Job_arrival_time = sum(Entry_IT[i,:])+sum(Supply[i, :].*UnLoadTime)/2;
            OT_arrival_time = sum(Entry_OT[j,:]);
            Job_Storage_single = max(0,OT_arrival_time-Job_arrival_time);
            push!(Job_Storage,Job_Storage_single);
        end;
    end;
end;

req_num = 0;
for i in 1:length(Job_Storage)
    global req_num;
    if Job_Storage[i] == 0
        req_num += 1;
    end;
end;
println("Avg Job Storage time $(mean(Job_Storage))");
println("Maximum Job Storage time $(maximum(Job_Storage))");
println("Number of Jobs with no wait time $(req_num)");
println("Proportion of jobs with zero wait time $(req_num/length(Job_Storage))");
                      #--------------------------------------------------------------------------#
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