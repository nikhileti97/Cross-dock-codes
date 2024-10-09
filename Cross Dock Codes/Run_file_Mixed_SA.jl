#v1.3: updated to remove crossovers and upgrade SA to PBSA

import Pkg;
Pkg.add("Plots");
Pkg.add("Random");
Pkg.add("Distributions");
#Pkg.add("VegaLite");
Pkg.add("DataFrames");

using Plots;
using DataFrames;

n_max = 16;

#obj_list = zeros(Float64, n_max-30+1);
#Ar_list = zeros(Float64,n_max-30+1);

include("data_v1.2.jl");

include("Heuristic_v1.3.jl");

resultfile1 = open("Supply_Demand_Mixed_New.csv", "w");
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

n = n_doors;
#=
#while n <= n_max
    
if n % 2 == 0
    n_strip = Int(n/2);
    n_stack = Int(n/2);
else
    n_strip = div(n,2)+1
    n_stack = div(n,2)
end; 

Strip = Dict();
Stack = Dict();

#Ar_list[n-24] = Area = 15*7*15*(n_stack-6+1)

for i in 1:n_strip
    if i <= 6
        Strip[i] = (15*i,15*max(1,n_strip-6));
    else
        Strip[i] = (0,15*(i-6))
    end;
end;


for i in 1:n_stack
    if i <= 6
        Stack[i] = (15*(i-1),0);
    else
        Stack[i] = (90,15*(i-7));
    end;
end;


Distance_matrix = Array{Float64, 2}(undef, n_strip, n_stack);

for i in 1:n_strip, j = 1:n_stack
    Distance_matrix[i,j] = abs(Strip[i][1]-Stack[j][1]) + abs(Strip[i][2]-Stack[j][2])
end;

TransTime = (Distance_matrix)./(forklift_speed*60);
=#
#include("Heuristic_v1.3.jl");

Inb_wt1 = 0.5;

Inb_wt2 = 0.5; 

Out_wt3 = 0.5;

Out_wt4 = 0.5;

n_strip = Int(round(mean(UnLoadTime)/(mean(UnLoadTime)+mean(LoadTime))*n_doors));
n_stack = n_doors - n_strip;
TransTime = zeros(Float64,n_strip,n_stack);

#Temp_storage
#sol_time = @elapsed begin
Entr_IT, Leav_IT, Entr_OT, Leav_OT, Load_Tran, In_out= Heuristic_CDS(Inb_wt1, Inb_wt2, Out_wt3, Out_wt4,n_strip,n_stack,TransTime);
#end;

#include("schedule_gantt_heuristic.jl");

# stime = @benchmark Heuristic_CDS(Inb_wt1, Inb_wt2, Out_wt3, Out_wt4);
#
# sol_time = mean(stime.times)*10e-10;


#obj_value = maximum([sum(Leav_OT[j, :]) for j = 1:n_outbound]);
obj_value = 0;

Tardiness = zeros(Float64, n_outbound);
for j in 1:n_outbound
    global obj_value;
    Tardiness[j] = max(0, sum(Leav_OT[j, g] for g = 1:n_stack) - DeadLines[j]);
    if Tardiness[j] <= Breakpt[1]
        obj_value += Multiplier[1]*Penalty[j]*Tardiness[j];
    elseif Tardiness[j] <= Breakpt[2]
        obj_value += Multiplier[2]*Penalty[j]*Tardiness[j] + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]);
    else
        obj_value += Multiplier[3]*Penalty[j]*Tardiness[j] +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]);
    end
end;

cat_Entr_IT = zeros(Float64,n_inbound,n_stack)
Entr_IT = hcat(Entr_IT,cat_Entr_IT);
Leav_IT = hcat(Entr_IT,cat_Entr_IT);
cat_Entr_OT = zeros(Float64,n_outbound,n_strip);
Entr_OT = hcat(cat_Entr_OT,Entr_OT);
Leav_OT = hcat(cat_Entr_OT,Leav_OT);
#=
Strip = Dict();
iter1 = 1;
a = 1;
while a <= n_doors
    global a;
    global iter1;
    if a <= 6
        Strip[a] = (15*a,15*max(1,Int((n_doors-12)/2)));
        a = a+1;
    elseif (a >= 7) && (a <= 12)
        Strip[a] = (15*(a-7),0);
        a = a+1;
    else
        Strip[a] = (0,15*(iter1));
        Strip[a+1] = (90,15*(iter1-1));
        iter1 = iter1+1;
        a = a+2;
    end;
end;

Distance_matrix = Array{Float64, 2}(undef, n_doors, n_doors);

for i in 1:n_doors, j = 1:n_doors
    if i!=j
        Distance_matrix[i,j] = abs(Strip[i][1]-Strip[j][1]) + abs(Strip[i][2]-Strip[j][2]);
    else
        Distance_matrix[i,j] = 0;
    end;
end;

TransTime = (Distance_matrix)./(forklift_speed*60);
=#

TransTime = zeros(Float64,n_doors,n_doors);


#-----------------------------------(PBSA)-----------------------------------------------------------------------------#

function Temp_Storage_Mixed_plot(Enter_IT,Enter_OT,Leave_IT,Leave_OT,Load_Tran6)

    inbound_seq = Dict();
    outbound_seq = Dict();
    inbound_outbound = Dict();
    for i in 1:n_doors
        inbound_seq[i] = Array{Int32}(undef,0);
    end;
    for i in 1:n_doors
        outbound_seq[i] = Array{Int32}(undef,0);
    end;

    dummy_Leav_IT = copy(Leave_IT);

    for i in 1:n_doors
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
    Storage_Used = zeros(Int64,n_doors,makespan);
    Storage_Add_Used = zeros(Int64,n_doors,makespan);
    Storage_Sub_Used = zeros(Int64,n_doors,makespan);
    Inbound_Entry = [maximum(Enter_IT[i,:]) for i in 1:n_inbound];
    Inbound_Leave = [maximum(Leave_IT[i,:]) for i in 1:n_inbound];
    Outbound_Entry = [maximum(Enter_OT[i,:]) for i in 1:n_outbound];
    Outbound_Leave = [maximum(Leave_OT[i,:]) for i in 1:n_outbound];
    cum_Supply = zeros(Int64,n_doors);
    cum_remove = zeros(Int64,n_doors);
    for i in 1:n_doors
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
                        cum_Supply[i] += 1;
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


    for i in 1:n_doors
        prev_used = 0;
        for t in 1:makespan
            if t != 1
                prev_used = Storage_Used[i,t-1];
            end;
            Storage_Used[i,t] = prev_used+Storage_Add_Used[i,t]-Storage_Sub_Used[i,t];
            #=
            if Storage_Used[i,t] < 0
                Storage_Used[i,t] = 0;
            end;
            =#
        end;
    end;
    
    plot([1:makespan], Storage_Used[1,:], title = "Storage for Doors", label = "Door-1");
    prev_sum = Storage_Used[1,:];
    for i in 2:n_doors
        cum_sum = prev_sum + Storage_Used[i,:];
        plot!([1:makespan], cum_sum, label = "Door-$(i)")
        prev_sum = cum_sum;
    end;
    xlabel!("Time Elapsed");
    ylabel!("Temporary Storage Utilized in product units");
    savefig("Temp Storage Utilized-Mixed Doors_$(n_inbound)_$(n_outbound)_$(n_doors)");
    
end;

println("Inbound - $(n_inbound), Outbound - $(n_outbound), Doors - $(n_doors)");
include("Swaps_Mixed.jl");

OT_chrom, IT_chrom = sch_to_chrom1(Leav_IT,Leav_OT);
req_gate_seq = chrom_to_comb(OT_chrom,IT_chrom,Leav_OT,Leav_IT);
up_obj = fitness_fxn1(req_gate_seq);
populations = [20];
#populations = [15];
mut_iterations_list = [1000];
#mut_iterations_list = [250];
mut_temp = 4*up_obj/20;

final_best_seq = Dict();
final_best = 0;

mut_temp_grad = 0.90
sol_time2 = @elapsed begin
n_sample = 50;
Sol_Sample = zeros(Float64,length(populations),length(mut_iterations_list),n_sample);
i1 = 0;   
for population in populations
    global i1,final_best_seq,final_best;
    #local iter2;
    i1 = i1+1;
    i2 = 0;
    for mut_iterations in mut_iterations_list
        i2 = i2+1
        for iter in 1:n_sample
            sol_time = @elapsed begin
            PBSA_sol= PBSA1(req_gate_seq, up_obj, mut_iterations, mut_temp, mut_temp_grad, population);
            end;
            Sol_Sample[i1,i2,iter] = fitness_fxn1(PBSA_sol);
            if i1 == 1 && i2 == 1 && iter ==1
                resultfile = open("Mixed_Tuning_Swaps_New_$(n_inbound)_$(n_outbound)_$(n_doors).csv", "w");
                println(resultfile, "Inbound Trailers, Outbound Trailers, Doors,population, iterations, iter, Objective, Solution time");
                println(resultfile, "$(n_inbound), $(n_outbound), $(n_doors), $(population), $(mut_iterations), $(iter), $(Sol_Sample[i1,i2,iter]), $(sol_time)");
                println();
                close(resultfile);
                final_best = Sol_Sample[i1,i2,iter];
                final_best_seq = PBSA_sol;
            else
                resultfile = open("Mixed_Tuning_Swaps_New_$(n_inbound)_$(n_outbound)_$(n_doors).csv", "a");
                println(resultfile, "$(n_inbound), $(n_outbound), $(n_doors), $(population), $(mut_iterations), $(iter), $(Sol_Sample[i1,i2,iter]), $(sol_time)");
                println();
                close(resultfile);
                if Sol_Sample[i1,i2,iter] < final_best
                    final_best = Sol_Sample[i1,i2,iter];
                    final_best_seq = PBSA_sol;
                end;    
            end;
            Leave_IT,Entry_IT,Leave_OT,Entry_OT,Load_Dist,leftover1,left_Supply1,left_Demand1 = seq_to_sch(final_best_seq);
            obj_value2 = fitness_fxn1(final_best_seq);
            resultfile = open("Solution_file_Mix_Swap_New_$(n_inbound)_$(n_outbound)_$(n_doors).csv", "w");
            println(resultfile, "Trailer_num, Trailer_type, Arrival_time, Deadline, n_doors, n_inbound, n_outbound, n_products, Change_time, Supply, Demand, Penalty, BP1, BP2, M1, M2, M3, Assigned_DD, DD_entry_time, DD_leave_time, Tardiness, Objective_value, Solution_time(s)");
            for i in 1:n_inbound
                door_assign = sum(k for k = 1:n_doors if Leave_IT[i, k] > 0);
                println(resultfile, "$(i), Inbound trailer, $(Arrival[i]), NA, $(n_doors), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), $(sum(Supply[i,p] for p = 1:n_products)), NA, NA, $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(door_assign), $(Entry_IT[i,door_assign]), $(Leave_IT[i,door_assign]), NA, $(obj_value2), $(sol_time)");
            end;
            for j in 1:n_outbound
                Tardiness[j] = max(0, sum(Leave_OT[j, g] for g = 1:n_stack) - DeadLines[j]);
                door_assign = sum(g for g = 1:n_doors if Leave_OT[j, g] > 0);
                println(resultfile, "$(j), Outbound trailer, $(0), $(DeadLines[j]), $(n_doors), $(n_inbound), $(n_outbound), $(n_products), $(Change_time), NA, $(sum(Demand[j,p] for p = 1:n_products)), $(Penalty[j]), $(Breakpt[1]), $(Breakpt[2]), $(Multiplier[1]), $(Multiplier[2]), $(Multiplier[3]), $(door_assign), $(Entry_OT[j, door_assign]), $(Leave_OT[j, door_assign]), $(Tardiness[j]), $(obj_value2), $(sol_time)");
            end;
            close(resultfile);
            Temp_Storage_Mixed_plot(Entry_IT,Entry_OT,Leave_IT,Leave_OT,Load_Dist);
                
            resultfile = open("Flows_heuristic_Mix_Swap_New_$(n_inbound)_$(n_outbound)_$(n_doors).csv", "w");
            println(resultfile, "Inbound Trailer, Outbound Trailer, Product, Quantity, Trans Ship Time");
            for i in 1:n_inbound
                str_dd = sum(k for k = 1:n_doors if Leave_IT[i, k] > 0);
                for j in 1:n_outbound
                    sta_dd = sum(g for g = 1:n_doors if Leave_OT[j, g] > 0);
                    for p in 1:n_products
                        if Load_Dist[i, j, p] != 0
                            println(resultfile, "IT$(i), OT$(j), P$(p), $(Load_Dist[i, j, p]), $(TransTime[str_dd, sta_dd])");
                        end;
                    end;
                end;
            end;
            close(resultfile);
                
            
        end;
    end;
end;
end;

