#v1.4: removed crossover function and upgraded mutation to PBSA

    function sch_to_chrom1(Leave_IT, Leave_OT)

        OT_time = copy(Leave_OT);

        IT_time = copy(Leave_IT);

        OT_chrom = [(0,0) for x = 1:n_outbound];

        IT_chrom = [(0,0) for x = 1:n_inbound];

        for door in 1:n_doors
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
                     max_OT = findmax(OT_time[:, door])[2];
                     # tuples of OT_chrom only mentioned for outbound trucks - (door utilised,position of truck)
                     OT_chrom[max_OT] = (door, pos);
                     pos -= 1;
                     OT_time[max_OT, door] = 0;
                end;
            end;

        end;

        for door in 1:n_doors

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
    #=
    #mutation functions
    function swap_trailers(Chrom_ST)
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
    =#
    function swap_trailers1(swap_sq)
        swap_gate_sq = deepcopy(swap_sq);    
        veh_per_gate = zeros(Int8,length(swap_gate_sq));
        num = 0;
        for gate in sort(collect(keys(swap_gate_sq)))
            veh_per_gate[gate]=length(swap_gate_sq[gate]);
            if length(swap_gate_sq[gate]) > 0
                num = num+1
            end;
        end;
        if num == 1
            ind = findall(x -> x > 0, veh_per_gate)[1]
            gene1 = ind;
            gene2 = ind;
        else
            gene1 = rand(1:length(swap_gate_sq));
            gene2 = rand(1:length(swap_gate_sq));
            while length(swap_gate_sq[gene1]) < 1
                gene1 = rand(1:length(swap_gate_sq));
            end;
            while length(swap_gate_sq[gene2]) < 1
                gene2 = rand(1:length(swap_gate_sq));
            end;
        end;


        gene11 = rand(1:length(swap_gate_sq[gene1]));
        gene22 = rand(1:length(swap_gate_sq[gene2]));
        #=
        if gene1 == gene2
            while gene22 == gene11
                gene22 = rand(1:length(swap_gate_sq[gene2]));
            end;
        end;
        =#
        temp_gene = swap_gate_sq[gene1][gene11];
        swap_gate_sq[gene1][gene11] = swap_gate_sq[gene2][gene22];
        swap_gate_sq[gene2][gene22] = temp_gene;

        return swap_gate_sq;
    end;
    
    function insert_trailer1(gate_sq)
        req_sq = deepcopy(gate_sq);
        gene1 = rand(1:length(req_sq));
        gene2 = rand(1:length(req_sq));

        # to pick a non-zero gate 
        while length(req_sq[gene1]) < 1
            gene1 = rand(1:length(req_sq));
        end;

        # pick a different gate
        while gene1 == gene2
            gene2 = rand(1:length(req_sq));
        end;
        # pick any truck from gate 1
        gene11 = rand(1:length(req_sq[gene1]));
        temp = req_sq[gene1][gene11];
        # delete truck from gate 1
        deleteat!(req_sq[gene1], findall(x -> x == temp, req_sq[gene1]));
        # push truck in gate 2
        push!(req_sq[gene2],temp);
        if length(req_sq[gene2]) > 1
            prev_seq=copy(req_sq[gene2])
            gene22 = rand(1:length(req_sq[gene2])-1)
            req_sq[gene2][gene22] = temp;
            i = gene22+1;
            while i <= length(req_sq[gene2])
                req_sq[gene2][i] = prev_seq[i-1];
                i = i+1;
            end;
        end;

        return req_sq;

    end;
    
    function chrom_to_comb(OT_chrom,IT_chrom,Leave_OT,Leave_IT)
    
        veh_per_gate_in = zeros(Int8, n_doors);
        veh_per_gate_out = zeros(Int8, n_doors);
        OT_time = copy(Leave_OT);
        IT_time = copy(Leave_IT);
        # no. of trailers for each gate
        for IT in 1:n_inbound
            gate = IT_chrom[IT][1];
            veh_per_gate_in[gate] += 1;
        end;

        for OT in 1:n_outbound
            gate = OT_chrom[OT][1];
            veh_per_gate_out[gate] += 1;
        end;

        gate_seq_in = Dict();
        gate_seq_out = Dict();
        for gate in 1:n_doors
            gate_seq_in[gate] = zeros(Int8, veh_per_gate_in[gate]);
            gate_seq_out[gate] = zeros(Int8, veh_per_gate_out[gate]);
        end;
        # arrange trailers in ascending order
        for gate in 1:n_doors
            num = veh_per_gate_in[gate];
            if num > 0
                while maximum(IT_time[:,gate]) > 0
                    max_IT = findmax(IT_time[:,gate])[2]
                    gate_seq_in[gate][num] = max_IT;
                    num = num-1;
                    IT_time[max_IT,gate] = 0;
                end;
            end;
        end;

        for gate in 1:n_doors
            num1 = veh_per_gate_out[gate];
            if num1 > 0
                while maximum(OT_time[:,gate]) > 0
                    max_OT = findmax(OT_time[:,gate])[2]
                    gate_seq_out[gate][num1] = max_OT;
                    num1 = num1-1;
                    OT_time[max_OT,gate] = 0;
                end;
            end;
        end;
        # Ordered sequence of trailers for each gate
        gate_seq = Dict();
        prev_gate_seq = Dict();
        for gate in 1:n_doors
            # get total number of vehicles
            #print("--------$(gate)-------\n")
          tot_veh = veh_per_gate_in[gate]+veh_per_gate_out[gate];
          gate_seq[gate] = [("In",0) for x = 1:tot_veh]
          num = 1;
            # initialize ITs to the comb. seq
          for IT in gate_seq_in[gate]
            gate_seq[gate][num] = ("In",IT)
            num = num+1;
          end;
            # create a duplicate for copying
          prev_gate_seq[gate] =  copy(gate_seq[gate]);
          num = 0;
          num3 = 0;
          for OT in gate_seq_out[gate]
              #print("$(OT)\n");
              if Leave_OT[OT,gate] < maximum(Leave_IT[:,gate])
                num3 = num3+1;
                for IT in gate_seq_in[gate]
                    if Leave_OT[OT,gate] < Leave_IT[IT,gate]
                      req_pos = findall(x -> x == ("In",IT) ,prev_gate_seq[gate])[1];
                      temp = ("In",IT);
                      gate_seq[gate][req_pos] = ("Out",OT);
                      num2 = req_pos+1;
                      while num2 <= tot_veh
                          gate_seq[gate][num2] = prev_gate_seq[gate][num2-1];
                          num2 = num2+1;
                      end;
                      break;
                    end;
                end;
              else
                gate_seq[gate][veh_per_gate_in[gate]+num3+1+num] = ("Out",OT);
                num = num+1;
              end;
              prev_gate_seq[gate] = copy(gate_seq[gate]);
              #print("--$(prev_gate_seq[gate])")
          end;
        end;

        return gate_seq

    end;

    function seq_to_sch(gate_seq)
        
        sch_seq = copy(gate_seq);
        Entry2_OT = zeros(Float64, n_outbound, n_doors);
        Leave2_OT = zeros(Float64, n_outbound, n_doors);
        Entry2_IT = zeros(Float64, n_inbound, n_doors);
        Leave2_IT = zeros(Float64, n_inbound, n_doors);
        Demand_iter = copy(Demand);
        Supply_iter = copy(Supply);
        Load_Tran3 = zeros(Int32, n_inbound, n_outbound, n_products);
        Unsch_OT = [j for j in 1:n_outbound];
        In_gate = [Array{Int8}(undef,0) for x = 1:length(sch_seq)];
        OT_iter = zeros(Int8, length(sch_seq))
        veh_per_gate_out = zeros(Int8, n_doors);
        for gate in sort(collect(keys(gate_seq)))
            num = 0;
            for obj in req_gate_seq[gate]
                if obj[1] == "Out"
                    num = num+1;
                end;
            end;
            veh_per_gate_out[gate] = num;
        end;
        max_iteration = maximum(veh_per_gate_out)
        #OT_iter = [Array{Int8}(undef,0) for x = 1:length(sch_seq)];
        Dem_req = Array{Int8}(undef,0)
        Sup_req = Array{Int8}(undef,0)
        prev_Sup_req = Array{Int8}(undef,0)

        # ------------------ outer Iteration unscheduled OT possibly starting from here -------------------------------#
        # ------------------ inner iteration unscheduled IT possibly starting from here--------------------------------#
        # get all inbound trucks before the outbound truck in other iterations
        iter=1;
        OT_iter_prev = copy(OT_iter);
        sup_flag = 0;
        rem_dem = 0;
        while length(Unsch_OT) > 0
            #print("\n----------------------------$(iter)------------------------------------------\n")
            #print("\nList of previous inbound trucks - $(prev_Sup_req)\n") 
            #track = Array{Int8}(undef,0);
            # initialize new unscheduled list of inbound trucks
            In_gate = [Array{Int8}(undef,0) for x = 1:length(sch_seq)];
            # get all inbound trucks before the first outbound truck
            if iter == 1
                for gate in sort(collect(keys(sch_seq)))
                    for obj in sch_seq[gate]
                        if obj[1] == "In"
                            push!(In_gate[gate],obj[2])
                        else
                            OT_iter[gate]=obj[2]
                            push!(Dem_req,obj[2])
                            break
                        end
                    end
                end
            else
                for i in 1:length(OT_iter)
                    if OT_iter[i] == 0
                        # Exception for index to tally the last outbound truck at gate
                        if OT_iter_prev[i] == 0
                            continue;
                        else
                            index = findall(x -> x == ("Out",OT_iter_prev[i]),sch_seq[i])[1]
                            for obj in sch_seq[i][index+1:end]
                                if obj[1] == "In"
                                    push!(In_gate[i],obj[2])
                                else
                                    OT_iter[i]=obj[2]
                                    push!(Dem_req,obj[2])
                                    break
                                end;
                            end;
                        end;
                    end;
                end;
            end;
            #print("\n$(In_gate)\n")
            OT_iter_prev = copy(OT_iter);
            #get entry and leave times of inbound trucks based on previous inbound and outbound trucks
            for pos in 1:length(In_gate)
                if iter == 1
                    if length(In_gate[pos]) > 0
                        num = 1;
                        while num <= length(In_gate[pos])
                            if num == 1;
                                IT = In_gate[pos][num]
                                Entry2_IT[IT,pos] = Arrival[IT];
                                Leave2_IT[IT,pos] = Arrival[IT]+sum(Supply[IT, :].*UnLoadTime);
                            else
                                prev_IT = In_gate[pos][num-1]
                                IT = In_gate[pos][num]
                                Entry2_IT[IT,pos] = max(Arrival[IT],Leave2_IT[prev_IT,pos]+Change_time);
                                Leave2_IT[IT,pos] = Entry2_IT[IT,pos]+sum(Supply[IT, :].*UnLoadTime);
                            end;
                            push!(Sup_req,IT);
                            num = num+1;
                        end;
                    end;
                else
                    if length(In_gate[pos]) > 0
                        num1 = 1;
                        while num1 <= length(In_gate[pos])
                            if num1 == 1;
                                IT = In_gate[pos][num1]
                                Entry2_IT[IT,pos] = max(Arrival[IT],maximum(Leave2_OT[:,pos])+Change_time);
                                Leave2_IT[IT,pos] = Entry2_IT[IT,pos]+sum(Supply[IT, :].*UnLoadTime);
                            else
                                prev_IT = In_gate[pos][num1-1]
                                IT = In_gate[pos][num1]
                                Entry2_IT[IT,pos] = max(Arrival[IT],Leave2_IT[prev_IT,pos]+Change_time);
                                Leave2_IT[IT,pos] = Entry2_IT[IT,pos]+sum(Supply[IT, :].*UnLoadTime);
                            end;
                            num1 = num1+1;
                            push!(Sup_req,IT);
                        end;
                    end;
                end;
            end;


            if length(Sup_req) == length(prev_Sup_req)
                sort_Sup_req = sort(Sup_req);
                sort_prev_req = sort(prev_Sup_req);
                if length(sort_Sup_req) == 0
                    sup_flag = 1;
                else
                    for ind in 1:length(sort_Sup_req)
                        sup_flag = 1;
                        if sort_Sup_req[ind] != sort_prev_req[ind]
                            sup_flag = 0;
                            break;
                        end;
                    end;
                end;
            end;

            for i in 1:length(OT_iter)
                if OT_iter_prev[i] != OT_iter[i]
                    sup_flag = 0;
                    break;
                end;
            end;

            #=
            #print("\n$(sup_flag)\n")
            #print("\n$(rem_dem)\n")
            =#
            # Dem_req is sequence of trucks
            max_iter = length(Dem_req)
            iterate = 0;

            #print("List of Scheduled Inbound Trucks - \n$(Sup_req)\n")
            #print("\nList of Considered Demand Trucks - $(OT_iter)\n")


            if sup_flag == 1 && rem_dem > 0
                break;
            end;

            sup_flag = 0;
            #print("\n List of Current Outbound Trucks - $(OT_iter)\n")
            #print("\n$(Dem_req)\n")
            while length(Dem_req) > 0
                if iterate > max_iter
                    break;
                end;
                # need a break statement
                AIT = Dict();
                Departure = Dict();
                Dep_deadl = Dict();
                flag = 0;
                for OT in Dem_req
                    AIT[OT] = [];
                    IT_list = [IT for IT in Sup_req];
                    Demand_iter_OT = copy(Demand_iter[OT,:]);

                    while sum(Demand_iter_OT) > 0

                        Prod_avail_time = Dict();
                        sub_flag = 0;

                        for IT in IT_list
                            #Supply_iter is updated at the end of the current loop 
                            if sum([min(Demand_iter_OT[p], Supply_iter[IT,p]) for p in 1:n_products]) > 0

                                Entry_IT = sum(Entry2_IT[IT, :]);

                                Prod_avail_time[IT] = Entry_IT + sum(Supply[IT, :].*UnLoadTime)/2;

                            end;
                            # no matched products between any available inbound trucks and the outbound truck
                            if sum([min(Demand_iter_OT[p], Supply_iter[IT,p]) for p in 1:n_products]) == 0
                                sub_flag = sub_flag +1;
                            end;

                        end;
                        # find list of IT or AIT[OT] with earliest availability for each OT. 
                        # no available transfer of goods and demand not satisfied -> escape the loop of meeting the demand
                        # no AIT
                        #=
                        if length(IT_list) == sub_flag & sum(Demand_iter_OT) > 0
                            flag = 1;
                            break;
                        end;
                        =#
                        if length(Prod_avail_time) == 0
                            flag = 1;
                            break;
                        end;

                        next_IT = findmin(Prod_avail_time)[2];

                        push!(AIT[OT], next_IT);

                        Demand_iter_OT = Demand_iter_OT - Supply_iter[next_IT,:];
                        Demand_iter_OT[Demand_iter_OT .< 0] .= 0;

                        deleteat!(IT_list, findall(x -> x == next_IT, IT_list));
                        # end of while loop to meet demand
                    end;

                    #calculate the criteria values for OT selection strategy
                    #if previous loop broke, continue for next truck
                    # Loop broke -> no departure, delay times
                    if flag == 1;
                        flag = 0;
                        continue;
                    else
                        #else calculate the departure time, delay of the outbound truck. 
                        OT_DD = findall(x -> x == OT,OT_iter)[1];

                        Leave_AIT = Dict();
                        Dem_copy = copy(Demand[OT, :]);
                        for IT in AIT[OT]
                            IT_DD = findmax(Entry2_IT[IT,:])[2];
                            Load_trans = [min(Dem_copy[p], Supply_iter[IT, p]) for p in 1:n_products];
                            Leave_AIT[IT] = sum(Entry2_IT[IT, :]) + sum(Supply[IT, :].*UnLoadTime)/2 + TransTime[IT_DD, OT_DD] + sum(Load_trans.*LoadTime);
                            Dem_copy = Dem_copy - Load_trans;
                        end;

                        # for pos = 2,3,4 so on, leave_OT increases gradually with each iteration
                        # if loop involving number of inbound trailers between
                        temp_Entry_OT = max(maximum(Leave2_IT[:, OT_DD]),maximum(Leave2_OT[:, OT_DD])) + Change_time;
                        OT_min_load = sum(Demand[OT, :].*LoadTime);
                        Departure[OT] = max(temp_Entry_OT + OT_min_load, findmax(Leave_AIT)[1]);

                        if Departure[OT] > DeadLines[OT]
                            Dep_deadl[OT] = Departure[OT] - DeadLines[OT];
                        else
                            Dep_deadl[OT] = 0;
                        end;
                    end;

                end;

                # Unable to get departure times of any outbound trucks due to unsatisfied demand
                # push those trucks into the next iteration.
                if length(Departure) == 0
                    #sup_flag = 1;
                    break;
                end;
                    #select an OT based on OT selection strategy
                    # OT selected from set of OT's with same position
                if findmax(Dep_deadl)[1] > 0
                    next_OT = findmax(Dep_deadl)[2];
                else
                    next_OT = findmin(Departure)[2];
                end;
                #print("$(next_OT)\n")
                # OT_DD find the door for the next scheduled truck
                OT_DD = findall(x -> x == next_OT,OT_iter)[1];
                OT_iter[OT_DD] = 0;
                deleteat!(Dem_req, findall(x -> x == next_OT, Dem_req));
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
                iterate = iterate+1;
            end;

            #=
            Sup_iter = [Array{Int8}(undef,0) for x = 1:length(Sup_req)];
            Dem_iter = [Array{Int8}(undef,0) for x = 1:length(OT_iter)];
            ind = 0;
            for IT in Sup_req
                ind = ind+1
                for p in 1:n_products
                    push!(Sup_iter[ind],Supply_iter[IT,p]);
                end;
            end;
            =#
            #prev_Sup_req = copy(Sup_req)
            #print("\nList of updated inbound trucks - $(prev_Sup_req)\n")
            IT_dummy = [IT for IT in Sup_req];
            for IT in IT_dummy
                if sum(Supply_iter[IT,:]) == 0
                    #print("\n$(IT)\n")
                    deleteat!(Sup_req, findall(x -> x == IT,Sup_req))
                end;
            end;
            iter = iter+1;

            Sup_sum = Array{Int32}(undef,0);
            for IT in Sup_req
                push!(Sup_sum,sum(Supply_iter[IT,:]))
            end;

            Dem_sum = Array{Int32}(undef,0);
            for OT in Dem_req
                push!(Dem_sum,sum(Demand_iter[OT,:]))
            end;
            rem_dem = sum(Dem_sum)

            #=
            print("\nLeft Over Supply of Supply Trucks - $(Sup_sum)\n")
            print("\nLeft Over Demand of Demand Trucks - $(Dem_sum)\n")

            print("\nList of Demand Trucks - $(Dem_req)\n")
            =#
            prev_Sup_req = copy(Sup_req)        
            #print("\nPrev - $(OT_iter_prev)\n")
            #print("\nCurrent - $(OT_iter)\n")
        end;

        return Leave2_IT,Entry2_IT,Leave2_OT,Entry2_OT,Load_Tran3,length(Unsch_OT),Supply_iter,Demand_iter


    end;
    

    function fitness_fxn1(gate_seq)
        copy_gate_seq = copy(gate_seq)
        Leave_IT,Entry_IT,Leave_OT,Entry_OT,Load_Tran4,leftover,left_Supply,left_Demand = seq_to_sch(copy_gate_seq);
        local_value = 0;
        if leftover > 0
            local_value = 1e8;
            #print("\nInfeasible Schedule\n")
        else
            Tardiness = zeros(Float64, n_outbound);
            #local_value = maximum([sum(Leave_OT[j, :]) for j = 1:n_outbound])
            for j in 1:n_outbound
                #global obj_value;
                Tardiness[j] = max(0, sum(Leave_OT[j, g] for g = 1:n_doors) - DeadLines[j]);
                if Tardiness[j] <= Breakpt[1]
                    local_value += Multiplier[1]*Penalty[j]*Tardiness[j];
                elseif Tardiness[j] <= Breakpt[2]
                    local_value += Multiplier[2]*Penalty[j]*Tardiness[j] + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]);
                else
                    local_value += Multiplier[3]*Penalty[j]*Tardiness[j] +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]);
                end
            end;
        end;

        return local_value
    end;
    
    function PBSA1(req_gate_seq, upd_obj, max_iterations, mut_temp, mut_temp_grad, population);
    
        best_seq = copy(req_gate_seq);
        Next_seq = copy(req_gate_seq);
        incumbent_sol = upd_obj;
        local_optimal = upd_obj;
        temp_curr = mut_temp;
        infeasible_iter = zeros(Int32, max_iterations);

        iter = 1
        while iter <= max_iterations

            sol = 1;
            while sol <= population
                neighbor_options_1 = ["Swap", "Do Nothing"];
                neighbor_options_2 = ["Swap", "Insert"];
                p1 = [1.0,0];
                p2 = [0.5,0.5];
                d1 = Categorical(p1);
                d2 = Categorical(p2);
                # select any options from neighbor if n_stack = 1 or n_stack > 1
                selection_OT = neighbor_options_2[rand(d2)];
                #=
                if n_doors == 1
                    selection_OT = neighbor_options_1[rand(d1)];
                else
                    selection_OT = neighbor_options_2[rand(d2)];
                end
                =#
                #=
                if n_doors == 1
                    selection_IT = neighbor_options_1[rand(1:2)];
                else
                    selection_IT = neighbor_options_2[rand(1:3)];
                end
                # try to not get "Do Nothing" option for both inbound and outbound trucks
                while (selection_OT == "Do Nothing") &&  (selection_IT == "Do Nothing")
                    if n_doors == 1
                        selection_OT = neighbor_options_1[rand(1:2)];
                    else
                        selection_OT = neighbor_options_2[rand(1:3)];
                    end

                    if n_doors == 1
                        selection_IT = neighbor_options_1[rand(1:2)];
                    else
                        selection_IT = neighbor_options_2[rand(1:3)];
                    end
                end;
                =#
                if selection_OT == "Swap"
                    Ne_seq_OT = swap_trailers1(best_seq);
                else
                    Ne_seq_OT = insert_trailer1(best_seq);
                end;
                #=
                if selection_IT == "Swap"
                    Ne_chrom_IT = swap_trailers(best_IT_sol);
                elseif selection_IT == "Insert"
                    Ne_chrom_IT = insert_trailer(best_IT_sol, "IT");
                else
                    Ne_chrom_IT = best_IT_sol;
                end;
                =#
                # get a set of neighbor solutuions based on the given best solution
                #print("\n$(sol)\n")

                #Neighbor_pop[sol] = Ne_seq_OT;
                curr_soln = fitness_fxn1(Ne_seq_OT);
                if curr_soln < local_optimal
                    Next_seq = Ne_seq_OT;
                    local_optimal = curr_soln;
                end;    

                if curr_soln == 1e8
                    infeasible_iter[iter] += 1;
                end;
                sol += 1;

            end;

            delta_sol = local_optimal - incumbent_sol

            if delta_sol > 0 && delta_sol < 10*temp_curr
                prob_accpt = exp(-delta_sol/temp_curr);
            else
                prob_accpt = 0;
            end;


            if delta_sol < 0 || prob_accpt >= rand(1)[1] 
                best_seq = Next_seq;
                incumbent_sol = local_optimal;
            end;

            temp_curr = mut_temp_grad*temp_curr;
            iter = iter+1;

        end;
    
        return best_seq
    end;



