#v1.3: updated to remove crossovers and upgrade SA to PBSA

import Pkg;
Pkg.add("Plots");
Pkg.add("Random");
Pkg.add("Distributions");

using Plots;
#using VegaLite, DataFrames;

n_max = 16;

#obj_list = zeros(Float64, n_max-30+1);
#Ar_list = zeros(Float64,n_max-30+1);

include("data_v1.2.jl");

include("Heuristic_v1.3.jl");


resultfile1 = open("Supply_Demand_Forklift_New.csv", "w");
cols = vcat(["Trailer_num", ", Trailer_type", ", Supply/Demand/Loadtime"], [", P$(p)" for p in 1:n_products]);
cols1 = join(cols)
println(resultfile1, cols1);
ltime = vcat(["NA", ", NA", ", Loadtime"], [", $(LoadTime[p])" for p in 1:n_products]);
ltime1 = join(ltime);
println(resultfile1, ltime1);
for i in 1:n_inbound
    entry = vcat(["IT$(i)", ", Inbound trailer", ", Supply"], [", $(Supply[i, p])" for p in 1:n_products]);
    entry1 = join(entry)
    println(resultfile1, entry1);
end;
for j in 1:n_outbound
    entry2 = vcat(["OT$(j)", ", Outbound trailer", ", Demand"], [", $(Demand[j, p])" for p in 1:n_products]);
    entry3 = join(entry2)
    println(resultfile1, entry3);
end;
close(resultfile1);
#=
include("Swaps_v1.4.jl");

population = 20;

mut_iterations = 250;

mut_temp = 600;

mut_temp_grad = 0.90 #[0.85,0.96] Ke-Lin Du et. al, Search and Optimization by metaheuristic, Springer 2016
=#

#n_doors = 12;
n = n_doors;


#include("Heuristic_v1.3.jl");

Inb_wt1 = 0.01;

Inb_wt2 = 0.99; 

Out_wt3 = 0.5;

Out_wt4 = 0.5;

n_strip = Int(round(mean(UnLoadTime)/(mean(UnLoadTime)+mean(LoadTime))*n_doors));
n_stack = n_doors - n_strip;
TransTime = zeros(Float64,n_strip,n_stack);

#Temp_storage
sol_time = @elapsed begin
Entr_IT1, Leav_IT1, Entr_OT1, Leav_OT1, Load_Tran1, In_out= Heuristic_CDS(Inb_wt1, Inb_wt2, Out_wt3, Out_wt4,n_strip,n_stack,TransTime);
end;

#include("schedule_gantt_heuristic.jl");

# stime = @benchmark Heuristic_CDS(Inb_wt1, Inb_wt2, Out_wt3, Out_wt4);
#
# sol_time = mean(stime.times)*10e-10;

#obj_value = maximum([sum(Leav_OT[j, :]) for j = 1:n_outbound]);

obj_value = 0;

Tardiness1 = zeros(Float64, n_outbound);
for j in 1:n_outbound
    global obj_value;
    Tardiness1[j] = max(0, sum(Leav_OT1[j, g] for g = 1:n_stack) - DeadLines[j]);
    if Tardiness1[j] <= Breakpt[1]
        obj_value += Multiplier[1]*Penalty[j]*Tardiness1[j];
    elseif Tardiness1[j] <= Breakpt[2]
        obj_value += Multiplier[2]*Penalty[j]*Tardiness1[j] + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]);
    else
        obj_value += Multiplier[3]*Penalty[j]*Tardiness1[j] +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]);
    end
end;
function Temp_Storage_forklift_plot(Enter_IT,Enter_OT,Leave_IT,Leave_OT,Load_Tran6,Inf_Entr0,Inf_Leav0,Start_Load0)

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
    
    xlabel!("Time Elapsed");
    ylabel!("Temporary Storage in Product units");
    savefig("Temp Storage Utilized-Separate Doors with forklift-Inexact");
    
end;

#-----------------------------------(PBSA)-----------------------------------------------------------------------------#

include("Integrated_Swaps.jl");

init_IT, init_OT = inbound_rand_seq(Leav_IT1,Leav_OT1);

populations = [20];
#populations = [15];
mut_iterations_list = [1000];
#mut_iterations_list = [250];
mut_temp = obj_value/20;
final_best = 0;
best_final_IT = zeros(Int32,n_inbound);
best_final_OT = zeros(Int32,n_outbound);

mut_temp_grad = 0.90
sol_time2 = @elapsed begin
n_sample = 10;
Sol_Sample = zeros(Float64,length(populations),length(mut_iterations_list),n_sample);
i1 = 0;   
for population in populations
    global i1,final_best,best_final_IT,best_final_OT;
    #local resultfile;
    i1 = i1+1;
    i2 = 0;
    for mut_iterations in mut_iterations_list
        i2 = i2+1
        for iter in 1:n_sample
            sol_time3 = @elapsed begin
            final_IT,final_OT = PBSA_sequence(init_IT, init_OT, mut_iterations, mut_temp, mut_temp_grad, population);
            Sol_Sample[i1,i2,iter] = calc_fitness(final_IT,final_OT);
            end;
            if i1 == 1 && i2 == 1 && iter ==1
                resultfile = open("Integrated_Tuning_Swaps_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "w");
                println(resultfile, "$(population), $(mut_iterations), $(iter), $(Sol_Sample[i1,i2,iter]), $(sol_time3)");
                println();
                close(resultfile);
                final_best = Sol_Sample[i1,i2,iter];
                best_final_IT = final_IT;
                best_final_OT = final_OT;
            else
                resultfile = open("Integrated_Tuning_Swaps_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "a");
                println(resultfile, "$(population), $(mut_iterations), $(iter), $(Sol_Sample[i1,i2,iter]), $(sol_time3)");
                println();
                close(resultfile);
                if Sol_Sample[i1,i2,iter] < final_best
                    best_final_IT = final_IT;
                    best_final_OT = final_OT;
                    final_best = Sol_Sample[i1,i2,iter];
                end;
            end;
            
            Entr_IT,Leav_IT,Entr_OT,Leav_OT,Load_Tran,Inf_Entr,Inf_Leav,Outf_Entr,Outf_Leav,Out_Entr,Out_Leav,Job_seq,Start_Load1 = Integrated_Decoder(best_final_IT,best_final_OT);
            
            Tardiness = zeros(Float64, n_outbound);
            for j in 1:n_outbound
                Tardiness[j] = max(0, sum(Leav_OT[j, g] for g = 1:n_stack) - DeadLines[j]);
            end;
            resultfile = open("Solution_file_Integrated_Swaps_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "w");
            println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_stack, n_strip, n_inbound, n_outbound, n_products, Change_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, Assigned_DD, DD_entry_time, DD_leave_time, Tardiness, Objective_value, Solution_time(s)");
            for i in 1:n_inbound
                door_assign = sum(k for k = 1:n_strip if Leav_IT[i, k] > 0);
                println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_stack), $(n_strip), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), Strip Door $(door_assign), $(Entr_IT[i,door_assign]), $(Leav_IT[i,door_assign]), NA, $(final_best), $(sol_time3)");
            end;
            for j in 1:n_outbound
                door_assign = sum(g for g = 1:n_stack if Leav_OT[j, g] > 0);
                println(resultfile, "OT$(j), Outbound trailer, $(0), $(DeadLines[j]), $(n_stack), $(n_strip), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), NA, $(sum(Demand[j,p] for p = 1:n_products)), $(Penalty[j]), $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), Stack Door $(door_assign), $(Entr_OT[j, door_assign]), $(Leav_OT[j, door_assign]), $(Tardiness[j]), $(final_best), $(sol_time3)");
            end;
            close(resultfile);
            
            resultfile = open("Solution_file_Forklift_Initial_Swap_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "w");
            println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_strip, n_stack, n_inbound, n_outbound, n_forklifts, n_products, C_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, Assigned_forklift, forklift_start_time, forklift_leave_time, Tardiness, Objective value");
            for i in 1:n_inbound
                forklift_assign = sum(k for k = 1:n_forklifts if Inf_Leav[i, k] > 1);
                println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_forklifts), $(n_products), $(C_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(forklift_assign), $(Inf_Entr[i,forklift_assign]), $(Inf_Leav[i,forklift_assign]), NA, $(final_best) ");
            end;

            for j in 1:n_outbound
                for g in 1:n_forklifts
                    cnt = 1
                    while cnt <= length(Out_Leav[(j,g)])
                        println(resultfile, "OT$(j), Outbound trailer, $(0), $(DeadLines[j]), $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_forklifts), $(n_products), $(C_time), NA, $(sum(Demand[j,p] for p = 1:n_products)), $(Penalty[j]), $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(g), $(Out_Entr[(j,g)][cnt]), $(Out_Leav[(j,g)][cnt]), $(Tardiness[j]), $(final_best)");
                        cnt = cnt+1;
                    end;
                end;
            end;
            close(resultfile);
            
            
            if i1 == length(populations) && i2 == length(mut_iterations_list) && iter == n_sample 
                Job_Storage = Array{Float64}(undef,0);
                Job_Delay = Array{Float64}(undef,0);
                Temp_Storage_forklift_plot(Entr_IT,Entr_OT,Leav_IT,Leav_OT,Load_Tran,Inf_Entr,Inf_Leav,Start_Load1);
                for i in 1:n_inbound
                    for j in 1:n_outbound
                        if sum(Load_Tran[i,j,:]) > 0
                            Job_arrival_time = sum(Entr_IT[i,:])+sum(Supply[i, :].*UnLoadTime)/2;
                            OT_arrival_time = sum(Entr_OT[j,:]);
                            Job_Storage_single = max(0,OT_arrival_time-Job_arrival_time);
                            push!(Job_Storage,Job_Storage_single);
                        end;
                    end;
                end;
                local req_num = 0;
                for i in 1:length(Job_Storage)
                    if Job_Storage[i] == 0
                        req_num += 1;
                    end;
                end;
                
                println("Avg Job Storage time $(mean(Job_Storage))");
                println("Maximum Job Storage time $(maximum(Job_Storage))");
                println("Number of Jobs with no wait time $(req_num)");
                println("Proportion of jobs with zero wait time $(req_num/length(Job_Storage))");
            end;
            
            sol_time4 = @elapsed begin
                                                
            num_iterations = 1000;num_pop = 20; T_max = final_best/20; T_grad = 0.9;
         
            Entr_IT3,Leav_IT3,Entr_OT3,Leav_OT3,better_obj,final_seq,prelim_list = Job_SA(Entr_IT,Leav_IT,Entr_OT,Leav_OT,Load_Tran,Job_seq,final_best,num_iterations,T_max,num_pop,T_grad);
                                        
            Inf_Entr, Inf_Leav, Outf_Entr, Outf_Leav, Out_Entr, Out_Leav, Entr_IT3, Leav_IT3, Entr_OT3, Leav_OT3, final1_seq = Job_Decoder(Entr_IT3,Leav_IT3,Entr_OT3,Leav_OT3,Load_Tran,final_seq);
            final_obj = better_obj;
            #final_obj = final_objective(Entr_IT3, Leav_IT3, Entr_OT3, Leav_OT3);
            end;
            if i1 == 1 && i2 == 1 && iter == 1
                resultfile = open("Forklift_Swaps_Objective_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "w");
                println(resultfile, "$(population), $(mut_iterations), $(iter), $(final_obj), $(sol_time4)");
                println();
                close(resultfile);
            else
                resultfile = open("Forklift_Swaps_Objective_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "a");
                println(resultfile, "$(population), $(mut_iterations), $(iter), $(final_obj), $(sol_time4)");
                println();
                close(resultfile);
            end;
                                        
            [Tardiness[j] = max(0, sum(Leav_OT3[j, g] for g = 1:n_stack) - DeadLines[j]) for j in 1:n_outbound];
            resultfile = open("Solution_file_Forklift_Swap_$(n_inbound)_$(n_outbound)_$(n_strip)_$(n_stack)_$(n_forklifts).csv", "w");
            println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_strip, n_stack, n_inbound, n_outbound, n_forklifts, n_products, C_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, Assigned_forklift, forklift_start_time, forklift_leave_time, Tardiness, Objective value");
            for i in 1:n_inbound
                forklift_assign = sum(k for k = 1:n_forklifts if Inf_Leav[i, k] > 1);
                println(resultfile, "IT$(i), Inbound trailer, $(Arrival[i]), NA, $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_forklifts), $(n_products), $(C_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(forklift_assign), $(Inf_Entr[i,forklift_assign]), $(Inf_Leav[i,forklift_assign]), NA, $(final_obj) ");
            end;

            for j in 1:n_outbound
                for g in 1:n_forklifts
                    cnt = 1
                    while cnt <= length(Out_Leav[(j,g)])
                        println(resultfile, "OT$(j), Outbound trailer, $(0), $(DeadLines[j]), $(n_strip), $(n_stack), $(n_inbound), $(n_outbound), $(n_forklifts), $(n_products), $(C_time), NA, $(sum(Demand[j,p] for p = 1:n_products)), $(Penalty[j]), $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(g), $(Out_Entr[(j,g)][cnt]), $(Out_Leav[(j,g)][cnt]), $(Tardiness[j]), $(final_obj)");
                        cnt = cnt+1;
                    end;
                end;
            end;
            close(resultfile);
            
            
        end;
    end;
end;
end;

