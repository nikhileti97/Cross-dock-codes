#v1.4: removed crossover function and upgraded mutation to PBSA
function Hybrid_meta(Entr_OT, Leav_OT, Entr_IT, Leav_IT, obj_value, mut_iterations, mut_temp, mut_temp_grad, population)

    function sch_to_chrom(Leave_IT, Leave_OT)

        OT_time = copy(Leave_OT);

        IT_time = copy(Leave_IT);

        OT_chrom = [(0,0) for x = 1:n_outbound];

        IT_chrom = [(0,0) for x = 1:n_inbound];

        for door in 1:n_stack
            # calculate number of trailers for each door 
            num_trail = 0;
            for OT in 1:n_outbound
                if Leave_OT[OT, door] > 0
                    num_trail += 1;
                end
            end
            # find positions of trailers at door from decreasing order. 
            if num_trail > 0
                pos = num_trail;
                 while maximum(OT_time[:, door]) > 0
                    # find latest truck at the door
                     max_OT = findmax(OT_time[:, door])[2];
                     # tuples of OT_chrom only mentioned for outbound trucks - (door utilised,position of truck)
                     OT_chrom[max_OT] = (door, pos);
                     pos -= 1;
                     OT_time[max_OT, door] = 0;
                end;
            end;

        end;

        for door in 1:n_strip

            num_trail = 0;
            for IT in 1:n_inbound
                if Leave_IT[IT, door] > 0
                    num_trail += 1;
                end
            end

            if num_trail > 0
                pos = num_trail;
                 while maximum(IT_time[:, door]) > 0
                     max_IT = findmax(IT_time[:, door])[2];
                     IT_chrom[max_IT] = (door, pos);
                     pos -= 1;
                     IT_time[max_IT, door] = 0;
                end;
            end;

        end;

        return OT_chrom, IT_chrom;

    end;

    #mutation functions
    function swap_trailers(Chrom_ST)
        # pick any two parts or gene of a chromosome and swap it
        Chrom = copy(Chrom_ST);
        gene1 = rand(1:length(Chrom));
        gene2 = rand(1:length(Chrom));

        while gene2 == gene1
            gene2 = rand(1:length(Chrom));
        end;

        temp_gene = Chrom[gene1];
        Chrom[gene1] = Chrom[gene2];
        Chrom[gene2] = temp_gene;

        return Chrom;

    end;

    function insert_trailer(Chrom_SO, type)
        
        Chrom = copy(Chrom_SO);
        gene_to_swap = rand(1:length(Chrom));
        gene_gate = Chrom[gene_to_swap][1];
        gene_pos = Chrom[gene_to_swap][2];

        if type == "IT"
            num_trailer = n_inbound;
            n_doors = n_strip;
        else
            num_trailer = n_outbound;
            n_doors = n_stack;
        end;
        # obtain number of doors, num_trailer = n_inbound/n_outbound
        veh_per_gate = zeros(Int8, n_doors);

        for trailer in 1:num_trailer
            gate = Chrom[trailer][1];
            veh_per_gate[gate] += 1;
        end;
        # number of trailers in each gate
        avail_gates = [x for x in 1:n_doors if x != gene_gate];
        # sel_gate is not gene_gate
        sel_gate = avail_gates[rand(1:length(avail_gates))];
        # extra vehicle on sel_gate
        avail_pos = veh_per_gate[sel_gate] + 1;

        sel_pos = rand(1:avail_pos);
        # one outbound truck(gene_to_swap) moved to another gate and another position
        Chrom[gene_to_swap] = (sel_gate, sel_pos);
        # change positions for trailers other than the selected gene_to_swap trailer
        for trailer in 1:length(Chrom)
            # at the selcted gate or inserted gate, shift all subsequent trailers backward at sel_gate by 1
            if (trailer != gene_to_swap) && (Chrom[trailer][1] == sel_gate) && (Chrom[trailer][2] >= sel_pos)
                Chrom[trailer] = (sel_gate, Chrom[trailer][2] + 1);
            #at the gate from which trailer is removed, shift all subsequent trailers forward at gene_gate by 1
            elseif (trailer != gene_to_swap) && (Chrom[trailer][1] == gene_gate) && (Chrom[trailer][2] >= gene_pos)
                Chrom[trailer] = (gene_gate, Chrom[trailer][2] - 1);
            end;
        end;

        return Chrom;

    end;

    function inb_chrom_to_sch(IT_chrom)

        veh_per_gate = zeros(Int8, n_strip);
        # no. of trailers each gate
        for IT in 1:n_inbound
            gate = IT_chrom[IT][1];
            veh_per_gate[gate] += 1;
                                    end;
        # array of trucks each gate
        gate_seq = Dict();
        for gate in 1:n_strip
            gate_seq[gate] = zeros(Int8, veh_per_gate[gate]);
        end;
        # provide sequence of inbound trucks each gate
        for IT in 1:n_inbound
            gate = IT_chrom[IT][1];
            pos = IT_chrom[IT][2];
            gate_seq[gate][pos] = IT;
        end;

        Entry_time = zeros(Float64, n_inbound, n_strip);
        Leave_time = zeros(Float64, n_inbound, n_strip);
        
        for gate in keys(gate_seq)
            #outer loop for each gate 
            pos = 1;
            for IT in gate_seq[gate]
                if pos == 1
                    # entry time and leave time of first inbound truck for each gate
                    Entry_time[IT, gate] = Arrival[IT];
                    Leave_time[IT, gate] = Entry_time[IT, gate] + sum(Supply[IT, :].*UnLoadTime);
                else
                    # entry and leave time of later inbound trucks ; pos = 2, pos-1 = 1; next_iteration pos = 3, pos-1 = 2
                    prev_IT = gate_seq[gate][pos-1];
                    # entry time is maximum of (arrival of current inbound truck, leave time of previous outbound inbound truck) 
                    Entry_time[IT, gate] = max(Arrival[IT], Leave_time[prev_IT, gate] + Change_time);
                    Leave_time[IT, gate] = Entry_time[IT, gate] + sum(Supply[IT, :].*UnLoadTime);
                end;
                pos += 1;
            end;
        end;

        return Entry_time, Leave_time

    end;

    function out_chrom_to_sch(OT_chrom, IT_chrom, Entry_time, Leave_time)

        Position_dict = Dict();
        # sequence of outbound trucks in each position
        for OT in 1:n_outbound
            # dictionary of position and set of trucks at that position
            pos = OT_chrom[OT][2];
            if pos in keys(Position_dict)
                push!(Position_dict[pos], OT)
            else
                Position_dict[pos] = [OT];
            end;

        end;

        Entry2_OT = zeros(Float64, n_outbound, n_stack);

        Leave2_OT = zeros(Float64, n_outbound, n_stack);

        Demand_iter = copy(Demand);

        Supply_iter = copy(Supply);

        Load_Tran3 = zeros(Int32, n_inbound, n_outbound, n_products);

        for pos in sort(collect(keys(Position_dict)));
            # for each position
                                    
            Unsch_OT = Position_dict[pos];

            while length(Unsch_OT) > 0

                AIT = Dict();
                Departure = Dict();
                Dep_deadl = Dict();

                for OT in Unsch_OT

                    AIT[OT] = [];

                    IT_list = [IT for IT in 1:n_inbound];

                    Demand_OT = copy(Demand_iter[OT, :]);

                    #create AIT for each OT at position "pos"
                    while sum(Demand_OT) > 0

                        Prod_avail_time = Dict();

                        for IT in IT_list
                            #Supply_iter is updated at the end of the current loop 
                            if sum([min(Demand_OT[p], Supply_iter[IT,p]) for p in 1:n_products]) > 0

                                Entry2_IT = sum(Entry_time[IT, :]);

                                Prod_avail_time[IT] = Entry2_IT + sum(Supply[IT, :].*UnLoadTime)/2;

                            end;

                        end;
                        # find list of IT or AIT[OT] with earliest availability for each OT.
                        next_IT = findmin(Prod_avail_time)[2];

                        push!(AIT[OT], next_IT);

                        Demand_OT = Demand_OT - Supply_iter[next_IT,:];
                        Demand_OT[Demand_OT .< 0] .= 0;

                        deleteat!(IT_list, findall(x -> x == next_IT, IT_list));

                    end;

                    #calculate the criteria values for OT selection strategy
                    OT_DD = OT_chrom[OT][1];

                    Leave_AIT = Dict();
                    Dem_copy = copy(Demand[OT, :]);
                    for IT in AIT[OT]
                        IT_DD = IT_chrom[IT][1];
                        Load_trans = [min(Dem_copy[p], Supply_iter[IT, p]) for p in 1:n_products];
                        Leave_AIT[IT] = sum(Entry_time[IT, :]) + sum(Supply[IT, :].*UnLoadTime)/2 + TransTime[IT_DD, OT_DD] + sum(Load_trans.*LoadTime);
                        Dem_copy = Dem_copy - Load_trans;
                    end;

                    if pos == 1
                        temp_Entry_OT = 0;
                    else
                        # for pos = 2,3,4 so on, leave_OT increases gradually with each iteration
                        temp_Entry_OT = maximum(Leave2_OT[:, OT_DD]) + Change_time;
                    end

                    OT_min_load = sum(Demand[OT, :].*LoadTime);
                    Departure[OT] = max(temp_Entry_OT + OT_min_load, findmax(Leave_AIT)[1]);

                    if Departure[OT] > DeadLines[OT]
                        Dep_deadl[OT] = Departure[OT] - DeadLines[OT];
                    else
                        Dep_deadl[OT] = 0;
                    end;

                end;

                #select an OT based on OT selection strategy
                # OT selected from set of OT's with same position
                if findmax(Dep_deadl)[1] > 0
                    next_OT = findmax(Dep_deadl)[2];
                else
                    next_OT = findmin(Departure)[2];
                end;

                OT_DD = OT_chrom[next_OT][1];

                #update the unscheduled OTs at position "pos"
                deleteat!(Unsch_OT, findall(x -> x == next_OT, Unsch_OT));

                #calculate load transferred between OT and ITs in AITs, update supplies
                for IT in AIT[next_OT]
                    Load_trans = [min(Supply_iter[IT, p], Demand_iter[next_OT, p]) for p in 1:n_products];
                    Load_Tran3[IT, next_OT, :] = Load_trans;
                    # Update Supply and demand of trucks for AIT calculation of next outbound trucks
                    Supply_iter[IT, :] = Supply_iter[IT, :] - Load_trans;
                    Demand_iter[next_OT, :] = Demand_iter[next_OT, :] - Load_trans;
                end;

                #update OT Leave time
                Leave2_OT[next_OT, OT_DD] = Departure[next_OT];

                Entry2_OT[next_OT, OT_DD] = Leave2_OT[next_OT, OT_DD] - sum(Demand[next_OT, :].*LoadTime);
                # loop for scheduling OTs sharing same position
            end;
            # loop for traversing across positions
        end;

        return Entry2_OT, Leave2_OT, Load_Tran3;

    end;

    function fitness_fxn(OT_chrom, IT_chrom)

        Entry_IT, Leave_IT = inb_chrom_to_sch(IT_chrom);

        Leave_OT = out_chrom_to_sch(OT_chrom, IT_chrom, Entry_IT, Leave_IT)[2]; #we only need leave time

        #object_value = maximum([sum(Leave_OT[j, :]) for j = 1:n_outbound]);
        object_value = 0;

        Tardiness3 = zeros(Float64, n_outbound);
        Earliness3 = zeros(Float64, n_outbound);
        for j in 1:n_outbound
            Tardiness3[j] = max(0, sum(Leave_OT[j, g] for g = 1:n_stack) - DeadLines[j]);
            if Tardiness3[j] <= Breakpt[1]
                object_value += Multiplier[1]*Penalty[j]*Tardiness3[j];
            elseif Tardiness3[j] <= Breakpt[2]
                object_value += Multiplier[2]*Penalty[j]*Tardiness3[j] + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]);
            else
                object_value += Multiplier[3]*Penalty[j]*Tardiness3[j] +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]);
            end
        end;

        return object_value;

    end;

    function PBSA(Chrom_OT, Chrom_IT, obj_value, max_iterations, mut_temp, mut_temp_grad, population)

        incumbent_sol = copy(obj_value);
        best_OT_sol = copy(Chrom_OT); #starting solution from heuristic
        best_IT_sol = copy(Chrom_IT);
        iter = 1;
        temp_curr = mut_temp;
        sol_list = Array{Float64,1}(undef, max_iterations);

        while iter <= max_iterations

            Neighbor_pop = [([],[]) for i = 1:population];
            fitness_pop = Array{Float64,1}(undef, population);

            for sol in 1:population

                neighbor_options_1 = ["Swap", "Do Nothing"];
                neighbor_options_2 = ["Swap", "Insert", "Do Nothing"];
                # select any options from neighbor if n_stack = 1 or n_stack > 1
                if n_stack == 1
                    selection_OT = neighbor_options_1[rand(1:2)];
                else
                    selection_OT = neighbor_options_2[rand(1:3)];
                end

                if n_strip == 1
                    selection_IT = neighbor_options_1[rand(1:2)];
                else
                    selection_IT = neighbor_options_2[rand(1:3)];
                end
                # try to not get "Do Nothing" option for both inbound and outbound trucks
                while (selection_OT == "Do Nothing") &&  (selection_IT == "Do Nothing")
                    if n_stack == 1
                        selection_OT = neighbor_options_1[rand(1:2)];
                    else
                        selection_OT = neighbor_options_2[rand(1:3)];
                    end

                    if n_strip == 1
                        selection_IT = neighbor_options_1[rand(1:2)];
                    else
                        selection_IT = neighbor_options_2[rand(1:3)];
                    end
                end;
                #=
                selection_IT = "Do Nothing";
                selection_OT = "Do Nothing";
                =#
                # gene mutation from given chromosome(inbound and outbound)
                if selection_OT == "Swap"
                    Ne_chrom_OT = swap_trailers(best_OT_sol);
                elseif selection_OT == "Insert"
                    Ne_chrom_OT = insert_trailer(best_OT_sol, "OT");
                else
                    Ne_chrom_OT = best_OT_sol;
                end;

                if selection_IT == "Swap"
                    Ne_chrom_IT = swap_trailers(best_IT_sol);
                elseif selection_IT == "Insert"
                    Ne_chrom_IT = insert_trailer(best_IT_sol, "IT");
                else
                    Ne_chrom_IT = best_IT_sol;
                end;
                # get a set of neighbor solutuions based on the given best solution
                Neighbor_pop[sol] = Ne_chrom_OT, Ne_chrom_IT;
                fitness_pop[sol] = fitness_fxn(Ne_chrom_OT, Ne_chrom_IT)

            end
            # find best solution among neighborhood
            current_sol = findmin(fitness_pop)[1] #value of min element
            sol_best_fitness = findmin(fitness_pop)[2] #index of min element
            # new chrome is current best chrome obtained from population
            New_chrom_OT, New_chrom_IT = Neighbor_pop[sol_best_fitness];
            # current best solution - previous best solution
            delta_sol = current_sol - incumbent_sol;
            #prob_accpt = 1;
            prob_accpt = exp(-delta_sol/temp_curr);

            if delta_sol < 0
                best_OT_sol = New_chrom_OT;
                best_IT_sol = New_chrom_IT;
                incumbent_sol = current_sol;
            # allow chance of accepting worse solutions with probability. Higher temperature, higher chance of accepting worse solutions
            elseif prob_accpt >= rand(1)[1]; 
                best_OT_sol = New_chrom_OT;
                best_IT_sol = New_chrom_IT;
                incumbent_sol = current_sol;
            end;
            # add the current solution to the list of solutions, incumbent solution is the value of the chromosome
            sol_list[iter] = incumbent_sol;
            temp_curr = mut_temp_grad*temp_curr;
            iter += 1;

        end;
        # sol_list[iter+1] is derived from sol_list[iter]
        plot([1:max_iterations], sol_list, title = "PBSA Convergence", label = "Incumbent_Sol");
        xlabel!("Iteration");
        ylabel!("Objective Value");
        savefig("PBSA_Converg.png");
        # returning one best chromosome after running (multiple iterations-250)
        return best_OT_sol, best_IT_sol;

    end;
    # 1) convert sequence of genes with given leave times of inbound and outbound trucks
    Chrom_OT, Chrom_IT =  sch_to_chrom(Leav_IT, Leav_OT);

    #println("************************Population Based Simulated Annealing************************")

    #println("Heuristic Solution: ", obj_value);
    # 2) best solution/final solution after number of iterations on initial heuristic solution
    PBSA_sol_OT, PBSA_sol_IT = PBSA(Chrom_OT, Chrom_IT, obj_value, mut_iterations, mut_temp, mut_temp_grad, population);
    # 3) final optimal objective value  
    PBSA_obj = fitness_fxn(PBSA_sol_OT, PBSA_sol_IT)

    #println("Solution after $(mut_iterations) PBSA iterations: ", PBSA_obj);
    # gene sequence to entry and leave times
    Entry_IT, Leave_IT = inb_chrom_to_sch(PBSA_sol_IT);
    Entry_OT, Leave_OT, Load_Tran4 = out_chrom_to_sch(PBSA_sol_OT, PBSA_sol_IT, Entry_IT, Leave_IT);

    return Entry_OT, Leave_OT, Entry_IT, Leave_IT, Load_Tran4, PBSA_obj;

end;
