#updated Supply Demand generation function to make sure supply = demand

import Pkg;
Pkg.add("JuMP");
Pkg.add("Random");
Pkg.add("Distributions");
Pkg.add("Gurobi");
#Pkg.add("VegaLite");
Pkg.add("DataFrames");
Pkg.add("CSV");
using Random,Distributions,Plots;

rng = MersenneTwister(1234);

# Assuming I-shaped cross-dock 90ft wide and 12 ft offset between dock-doors
#(3,5)
n_inbound = 20; # number of inbound trailers 'i'

n_outbound = 10; # number of outbound trailers 'j'

n_forklifts = 6;

C_time = 2; # forklift changeover time

n_doors = 6;



#=
if n % 2 == 0
    n_strip = Int(n/2);
    n_stack = Int(n/2);
else
    n_strip = div(n,2)+1
    n_stack = div(n,2)
end        
=#

n_products = 5; # number of products 'p'

HorizEnd = 48*60; # end of planning horizon

forklift_speed = 4; # ft/s or 7km/h

#(240-150) , (300-220)
Arrival = Array{Float64,1}(undef, n_inbound);  # arrival time of inbound trailer 'i'
Arrival = round.(rand(rng, Float64, n_inbound).*(900 - 150) .+ 150, digits = 0);
#Arrival = round.(rand(rng,Normal(195,20), n_inbound));
#Arrival = round.(rand(rng,TriangularDist(150, 240, 190), n_inbound));
#Arrival = round.(rand(rng,Exponential(200), n_inbound));
#Arrival = [203,219,201,227];
DeadLines = Array{Float64,1}(undef, n_outbound); #outbound trailers leave by this time
DeadLines = round.(rand(rng, Float64, n_outbound).*(1400 - 400) .+ 400, digits = 0); 

#DeadLines = round.(rand(rng,Normal(260,20), n_outbound));
#DeadLines = round.(rand(rng,TriangularDist(220, 300, 260), n_outbound));
#DeadLines = round.(rand(rng,Exponential(280), n_outbound));
Change_time = 30; #minimum time between two successive trailer dockings

LoadTime = Array{Float64,1}(undef, n_products);  # (W) unloading/loading time of product 'p'
LoadTime = round.(rand(rng, Float64, n_products).*(3.0 - 3.0) .+ 3.0, digits = 1);

UnLoadTime = Array{Float64,1}(undef, n_products);  # (W) unloading/loading time of product 'p'
UnLoadTime = round.(rand(rng, Float64, n_products).*(1.5 - 1.5) .+ 1.5, digits = 1);



#=
Distance_matrix = Array{Float64, 2}(undef, n_strip, n_stack);

TransTime = (Distance_matrix)./forklift_speed;
=#

Penalty = Array{Float64,1}(undef, n_outbound);  # (B) penalty for outbound trailer 'j' for leaving before loading all required products
Penalty = round.(rand(rng, Float64, n_outbound).*(2 - 2) .+ 3, digits = 0);

Breakpt = [120, 240]; # a

Multiplier = [1, 2, 4]; # lambda


function RandomSupDem_Gen(UnifRange_low, UnifRange_high, Matrix_dense)

    Supply = zeros(Int32, n_inbound, n_products);

    for i in 1:n_inbound
        allocation_1 = Array{Bool, n_products};
        prod_alloc = 0;
        while prod_alloc == 0
            allocation_1 = rand(rng, Bool, n_products);
            prod_alloc = sum(allocation_1);
        end;
        for p in 1:n_products
            if allocation_1[p] == 1
                Supply[i,p] = round(rand(rng, Float64)*(UnifRange_high - UnifRange_low) + UnifRange_low, digits = 0);
            else
                continue;
            end;
        end;
    end;

    for p in 1:n_products
        allocation_2 = Array{Bool, n_inbound};
        tra_alloc = sum(Supply[:, p])
        if tra_alloc == 0
            while tra_alloc == 0
                allocation_2 =  rand(rng, Bool, n_inbound);
                tra_alloc = sum(allocation_2);
            end;
            for i in 1:n_inbound
                if allocation_2[i] == 1
                    Supply[i,p] = round(rand(rng, Float64)*(UnifRange_high - UnifRange_low) + UnifRange_low, digits = 0);
                else
                    continue;
                end;
            end;
        else
            continue;
        end;
    end;

    allocation_3 = [if Supply[i, p] > 0 1 else 0 end for i = 1:n_inbound, p = 1:n_products];
    tot_alloc = sum(allocation_3)/length(allocation_3);
    if tot_alloc < Matrix_dense
        while tot_alloc < Matrix_dense
            unallocated = [(i,p) for i = 1:n_inbound, p = 1:n_products if allocation_3[i,p] == 0];
            count_unalloc = length(unallocated);
            selected_mem = rand(rng, 1:count_unalloc);
            Supply[unallocated[selected_mem][1], unallocated[selected_mem][2]] = round(rand(rng, Float64)*(UnifRange_high - UnifRange_low) + UnifRange_low, digits = 0);
            allocation_3 = [if Supply[i, p] > 0 1 else 0 end for i = 1:n_inbound, p = 1:n_products];
            tot_alloc = sum(allocation_3)/length(allocation_3);
        end;
    end;

    Demand = zeros(Int32, n_outbound, n_products);

    for j in 1:n_outbound
        allocation_4 = Array{Bool, n_products};
        p_alloc = 0;
        while p_alloc == 0
            allocation_4 = rand(rng, Bool, n_products);
            p_alloc = sum(allocation_4);
        end;
        for p in 1:n_products
            if allocation_4[p] == 1
                Demand[j, p] = round(rand(rng, Float64)*(UnifRange_high - UnifRange_low) + UnifRange_low, digits = 0);
            else
                continue;
            end;
        end;
    end;

    for p in 1:n_products
        allocation_5 = Array{Bool, n_outbound};
        t_alloc = sum(Demand[:, p])
        if t_alloc == 0
            while t_alloc == 0
                allocation_5 =  rand(rng, Bool, n_outbound);
                t_alloc = sum(allocation_5);
            end;
            for j in 1:n_outbound
                if allocation_5[j] == 1
                    Demand[j, p] = round(rand(rng, Float64)*(UnifRange_high - UnifRange_low) + UnifRange_low, digits = 0);
                else
                    continue;
                end;
            end;
        else
            continue;
        end;
    end;

    allocation_6 = [if Demand[j, p] > 0 1 else 0 end for j = 1:n_outbound, p = 1:n_products];
    tota_alloc = sum(allocation_6)/length(allocation_6);
    if tota_alloc < Matrix_dense
        while tota_alloc < Matrix_dense
            unallocated_d = [(j, p) for j = 1:n_outbound, p = 1:n_products if allocation_6[j, p] == 0];
            count_unalloc_d = length(unallocated_d);
            selected_mem_d = rand(rng, 1:count_unalloc_d);
            Supply[unallocated_d[selected_mem_d][1], unallocated_d[selected_mem_d][2]] = round(rand(rng, Float64)*(UnifRange_high - UnifRange_low) + UnifRange_low, digits = 0);
            allocation_6 = [if Demand[j, p] > 0 1 else 0 end for j = 1:n_outbound, p = 1:n_products];
            tota_alloc = sum(allocation_6)/length(allocation_6);
        end;
    end;

    for p in 1:n_products
        demand_imb = (sum(Demand[j, p] for j =1:n_outbound) - sum(Supply[i, p] for i =1:n_inbound))/sum(Demand[j, p] for j =1:n_outbound);
        Demand[:, p] = floor.((1 - demand_imb).*Demand[:, p]);
        gap = sum(Demand[:, p]) - sum(Supply[:, p]);
        allocated_dd = [j for j = 1:n_outbound if Demand[j, p] != 0];
        count_all  = length(allocated_dd);
        sel_mem = rand(rng, 1:count_all);
        while (Demand[allocated_dd[sel_mem], p] - gap) == 0
            sel_mem = rand(rng, 1:count_all);
        end;
        Demand[allocated_dd[sel_mem], p] = Demand[allocated_dd[sel_mem], p] - gap;
    end;

    return Supply, Demand;
end
#10,25
Supply, Demand = RandomSupDem_Gen(5, 20, 0.25);

Demand_max = Array{Int32,1}(undef, n_outbound); # (D_max_j) maximum demand
Demand_max = [sum(Demand[j,p] for p=1:n_products) for j = 1:n_outbound];

resultfile = open("Supply_Demand_exact_$(n_inbound)_$(n_outbound).csv", "w");
cols = vcat(["Trailer_num", ", Trailer_type", ",Arrival_time", ",Deadline", ", Supply/Demand/Loadtime"], [", P$(p)" for p in 1:n_products]);
cols1 = join(cols)
println(resultfile, cols1);
ltime = vcat(["NA", ", NA", ", min", ",min", ", Loadtime"], [", $(LoadTime[p])" for p in 1:n_products]);
ltime1 = join(ltime);
println(resultfile, ltime1);
for i in 1:n_inbound
    entry = vcat(["IT$(i)", ", Inbound trailer", ",$(Arrival[i])", ",NA", ", Supply"], [", $(Supply[i, p])" for p in 1:n_products]);
    entry1 = join(entry)
    println(resultfile, entry1);
end;
for j in 1:n_outbound
    entry2 = vcat(["OT$(j)", ", Outbound trailer", ", NA", ",$(DeadLines[j])", ", Demand"], [", $(Demand[j, p])" for p in 1:n_products]);
    entry3 = join(entry2)
    println(resultfile, entry3);
end;
close(resultfile);