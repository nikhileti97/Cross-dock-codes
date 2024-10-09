include("data_v1.2.jl");

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


rnd_tally = 1;
Door_Availability = zeros(Float64,n_strip);
copy_Arrival = copy(Arrival)
latest_availability = 0;
for i in 1:n_inbound
    global rnd_tally,copy_Arrival,Door_Availability,latest_availability;
    
    sel_truck = findmin(copy_Arrival)[2];
    if rnd_tally == 1
        latest_availability = maximum(Door_Availability);
        #println("$(Door_Availability)");
    end;
    Door_Availability[rnd_tally] = max(Arrival[sel_truck],latest_availability)+sum(Supply[sel_truck, :].*UnLoadTime)+Change_time;
    rnd_tally += 1;
    if rnd_tally > n_strip 
        rnd_tally = 1;
    end;
    copy_Arrival[sel_truck] = 1e5;
end;

M1 = findmax(Door_Availability)[1]-1*Change_time;
max_out = Int(ceil(n_outbound/n_stack));
Tot_Demand = [sum(Demand[j,p] for p in 1:n_products) for j in 1:n_outbound];
Tot_Demand = sort(Tot_Demand, rev=true);
M2 = M1 + sum(Tot_Demand[1:max_out]*LoadTime[1])+(max_out-2)*Change_time;
M2 = min(M2,HorizEnd);
println("M1 - $(M1); M2 - $(M2)");

M1 = HorizEnd;
M2 = HorizEnd;

println("Method Used - Exact - Convex Hull");
println("Inbound - $(n_inbound), Outbound - $(n_outbound), Strip - $(n_strip) Stack - $(n_stack)");

function Temp_Storage_Exact_Separate_plot(Enter_IT,Enter_OT,Leave_IT,Leave_OT,Load_Tran6)

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
                if Enter_IT[IT,i]+sum(Supply[IT,:].*UnLoadTime)/2 > maximum(Enter_OT[OT,:])
                    t_out = Enter_IT[IT,i]+sum(Supply[IT,:].*UnLoadTime)/2;
                else
                    t_out = maximum(Enter_OT[OT,:]);
                end;
                for p in 1:n_products
                    if Load_Tran6[IT,OT,p] > 0
                        for k in 1:Load_Tran6[IT,OT,p]
                            t_out += LoadTime[p];
                            cum_remove[i] += 1;
                            Storage_Sub_Used[i,Int(ceil(t_out))] += 1;
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
    xlabel!("Time Elapsed");
    ylabel!("Temporary Storage in Product units");
    savefig("Temp Storage Utilized-Exact-Convex-Hull-Separate Doors");
    
end;

Enter_IT = zeros(Float64, n_inbound, n_strip);
Leave_IT = zeros(Float64, n_inbound, n_strip);    
Enter_OT = zeros(Float64, n_outbound, n_stack);
Leave_OT = zeros(Float64, n_outbound, n_stack);
Load_Tran6 = zeros(Int32, n_inbound, n_outbound, n_products);


include("code_Separate_Convex_Hull.jl");


sol_time = @elapsed begin
optimize!(CDSP);
end;

if (termination_status(CDSP) == MOI.OPTIMAL) || ((termination_status(CDSP) == MOI.TIME_LIMIT) && has_values(CDSP))

    #include("schedule_gantt_exact_Convex_Hull.jl");

    resultfile = open("Solution_file_exact_Convex_Hull.csv", "w");
    println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_strip, n_stack, n_inbound, n_outbound, n_products, Change_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, Assigned_DD, DD_entry_time, DD_leave_time, Tardiness, Objective value, Termination Status , Solution_time");
    for i in 1:n_inbound
        door_assign = sum(k for k = 1:n_strip if value(InbAssign[i, k]) > 0.5);
        println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), Strip Door $(door_assign), $(value(EntryInbound[i])), $(value(LeaveInbound[i])), NA, $(objective_value(CDSP)), $(termination_status(CDSP)), $(sol_time)");
        Enter_IT[i,door_assign] = value(EntryInbound[i]);
        Leave_IT[i,door_assign] = value(EntryInbound[i])+sum(Supply[i,:].*UnLoadTime);
    end;
    for j in 1:n_outbound
        door_assign = sum(g for g = 1:n_stack if value(OutAssign[j, g]) > 0.5);
        println(resultfile, "OT$(j), Outbound trailer, $(0), $(DeadLines[j]), $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), NA, $(sum(Demand[j,p] for p = 1:n_products)), $(Penalty[j]), $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), Stack Door $(door_assign), $(value(EntryOutbound[j])), $(value(LeaveOutbound[j])), $(value(Tardiness[j])), $(objective_value(CDSP)), $(termination_status(CDSP)), $(sol_time)");
        Enter_OT[j,door_assign] = value(EntryOutbound[j]);
        Leave_OT[j,door_assign] = value(LeaveOutbound[j]);
    end;
    close(resultfile);

    resultfile = open("Flows_exact_Convex_Hull.csv", "w");
    println(resultfile, "Inbound Trailer, Outbound Trailer, Product, Quantity");
    for i in 1:n_inbound
        str_dd = sum(k for k = 1:n_strip if value(InbAssign[i, k]) == 1);
        for j in 1:n_outbound
            sta_dd = sum(g for g = 1:n_stack if value(OutAssign[j, g]) == 1);
            for p in 1:n_products
                if value(LoadTrans[i, j, p]) != 0
                    println(resultfile, "IT$(i), OT$(j), P$(p), $(value(LoadTrans[i, j, p]))");
                    Load_Tran6[i,j,p] = value(LoadTrans[i,j,p]);
                end;
            end;
        end;
    end;
    close(resultfile);
    Temp_Storage_Exact_Separate_plot(Enter_IT,Enter_OT,Leave_IT,Leave_OT,Load_Tran6);
    Job_Storage = Array{Float64}(undef,0);
    Job_Delay = Array{Float64}(undef,0);
    for i in 1:n_inbound
        for j in 1:n_outbound
            if sum(Load_Tran6[i,j,:]) > 0
                Job_arrival_time = sum(Enter_IT[i,:])+sum(Supply[i, :].*UnLoadTime)/2;
                OT_arrival_time = sum(Enter_OT[j,:]);
                Job_Storage_single = max(0,OT_arrival_time-Job_arrival_time);
                push!(Job_Storage,Job_Storage_single);
            end;
        end;
    end;
    local req_num = 0
    for i in 1:length(Job_Storage)
        if Job_Storage[i] == 0
            req_num += 1;
        end;
    end;
    println("Avg Job Storage time $(mean(Job_Storage))");
    println("Maximum Job Storage time $(maximum(Job_Storage))");
    println("Number of Jobs with no wait time $(req_num)");
    println("Proportion of jobs with zero wait time $(req_num/length(Job_Storage))");

else

    println("No Feasible solution found in 24 hours time limit")

end;
