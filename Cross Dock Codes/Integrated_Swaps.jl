function inbound_rand_seq(Leave_IT1,Leave_OT1)
    
    #IT_time = copy(Leave_IT1);
    IT_arrive = zeros(Float64,n_inbound);
    OT_arrive = zeros(Float64,n_outbound);
    key_array = zeros(Float64,n_inbound);
    key_out_array = zeros(Float64,n_outbound);
    raw_in_seq = zeros(Int32,n_inbound);
    raw_out_seq = zeros(Int32,n_outbound);
    for i in 1:n_inbound
        IT_arrive[i] = maximum(Leave_IT1[i,:]);
    end;
    for j in 1:n_outbound
        OT_arrive[j] = maximum(Leave_OT1[j,:]);
    end;
    place_array = rand(rng,Float64,n_inbound);
    place_out = rand(rng,Float64,n_outbound);
    
    cnt = n_inbound;

    #sort_array = sort(key_array);
    
    cnt = 1;
    while cnt <= n_inbound
        # find the last inbound truck,find the maximum of array
        IT_sel = findmax(IT_arrive)[2];
        raw_in_seq[n_inbound-cnt+1] = IT_sel;
        key_array[IT_sel] = findmax(place_array)[1];
        IT_arrive[IT_sel] = 0;
        place_array[findmax(place_array)[2]] = 0;
        cnt += 1;
    end;
    
    cnt = 1;
    while cnt <= n_outbound
        # find the last inbound truck,find the maximum of array
        OT_sel = findmax(OT_arrive)[2];
        raw_out_seq[n_outbound-cnt+1] = OT_sel;
        key_out_array[OT_sel] = findmax(place_out)[1];
        OT_arrive[OT_sel] = 0;
        place_out[findmax(place_out)[2]] = 0;
        cnt += 1;
    end;
    
    
    return raw_in_seq,raw_out_seq,key_array,key_out_array;
    
end;
       
function swap(req_seq)
    gene1 = rand(1:length(req_seq));
    gene2 = rand(1:length(req_seq));
    temp_seq = copy(req_seq);
    while gene1 == gene2
        gene2 = rand(1:length(req_seq));
    end;
    
    temp_truck = temp_seq[gene1];
    temp_seq[gene1] = temp_seq[gene2];
    temp_seq[gene2] = temp_truck;
    
    return temp_seq;
end;

function insert(req_seq)
        
    temp_seq = copy(req_seq);
    gene1 = rand(1:length(req_seq));
    gene2 = rand(1:length(req_seq));
    index = gene1;
    while index != gene2;
        if gene1 < gene2
            temp_truck = temp_seq[index];
            temp_seq[index] = temp_seq[index+1];
            temp_seq[index+1] = temp_truck;
            index += 1;
        else
            temp_truck = temp_seq[index];
            temp_seq[index] = temp_seq[index-1];
            temp_seq[index-1] = temp_truck;
            index -= 1;
        end;
    end;
    
    return temp_seq;
end;


function Integrated_Decoder(inbound_array,outbound_array)
    
    #in_priori = zeros(Int32,n_inbound);
    #out_priori = zeros(Int32,n_outbound);
    in_priori = copy(inbound_array);
    out_priori = copy(outbound_array);
    load_factor = 0.01;
    
    In_Entr = Dict(); # for possible consideration of multiple instances of same forklift to same truck
    In_Leav = Dict();
    Out_Entr = Dict();
    Out_Leav = Dict();
    AIT_time = Dict();
    Inf_Entr = zeros(Float64,n_inbound,n_forklifts); # Latest times available for forklifts
    Inf_Leav = zeros(Float64,n_inbound,n_forklifts);
    Outf_Entr = zeros(Float64,n_outbound,n_forklifts);
    Outf_Leav = zeros(Float64,n_outbound,n_forklifts);
    Avail_fork = zeros(Float64, n_forklifts);
    Stor_Load = zeros(Int32,n_inbound,n_products);
    Start_Load = zeros(Float64,n_inbound,n_outbound);
    rally = zeros(Int32,n_outbound);
    for i in 1:n_outbound
        for j in 1:n_forklifts
            Out_Entr[(i,j)] = Array{Float64}(undef,0);
            Out_Leav[(i,j)] = Array{Float64}(undef,0);
        end;
    end;
    
    for j in 1:n_outbound
        AIT_time[j] = Array{Float64}(undef,0);
    end;
    #=
    cnt = n_inbound;
    while cnt > 0
        IT_sel = findmax(in_seq)[2];
        in_priori[cnt] = IT_sel;
        in_seq[IT_sel] = 0;
        cnt -= 1;
    end;
    cnt = n_outbound
    while cnt > 0
        OT_sel = findmax(out_seq)[2];
        out_priori[cnt] = OT_sel;
        out_seq[OT_sel] = 0;
        cnt -= 1;
    end;
    =#
    Load_Trans = zeros(Int32, n_inbound, n_outbound, n_products);
    Entry_IT = zeros(Float64, n_inbound, n_strip);
    Leave_IT = zeros(Float64, n_inbound, n_strip);
    Entry_OT = zeros(Float64, n_outbound, n_stack);
    Leave_OT = zeros(Float64, n_outbound, n_stack);
    Earl_Avail_Strip = zeros(Float64, n_strip);
    Earl_Avail_Stack = zeros(Float64, n_stack);
    Stack_Door = zeros(Int32, n_outbound);
    Demand_iter = copy(Demand);
    Supply_iter = copy(Supply);
    # use in_priori
    cnt = 1;
    integrated_seq = Dict();
    Used_Doors = [];
    Free_Doors = [i for i in 1:n_stack];
    
    function Stack_Door_finder(Earl_Avail_Stack,Avail_fork,forklift_ind)
        Idle_Stack = Dict();
        req_door_compare = max(Avail_fork[forklift_ind],minimum(Earl_Avail_Stack))
        for k in 1:n_stack
            if req_door_compare >= Earl_Avail_Stack[k]
                Idle_Stack[k] =  req_door_compare - Earl_Avail_Stack[k];
            end;
        end;
        
        return findmin(Idle_Stack)[2];
    end;
    
    function Sub_Stack_Door_finder(Earl_Avail_Stack,Avail_fork,free_door,forklift_ind)
        Idle_Stack = Dict();
        req_door_compare = max(Avail_fork[forklift_ind],minimum(Earl_Avail_Stack[free_door]))
        for k in free_door
            if req_door_compare >= Earl_Avail_Stack[k]
                Idle_Stack[k] =  req_door_compare - Earl_Avail_Stack[k];
            end;
        end;
        
        return findmin(Idle_Stack)[2];
    end;
    
    function temp_place_forklift(req_time,Avail_fork)
        temp_forklift_avail = Dict();
        for i in 1:n_forklifts
            if req_time >= Avail_fork[i]
                temp_forklift_avail[i] = req_time-Avail_fork[i];
            end;
        end;
        
        return findmin(temp_forklift_avail)[2];
    
    end;

    
    
    
    function Forklift_finder(Earliest_Avail,door_index,Avail_fork)
        Idle_forklift = Dict();
        req_compare = max(Earliest_Avail[door_index],minimum(Avail_fork))
        for k in 1:n_forklifts
            if req_compare >= Avail_fork[k]
                Idle_forklift[k] = req_compare - Avail_fork[k];
            end;
        end;
        
        return findmin(Avail_fork)[2];
        
    end;
    # initialize set of scheduled and unscheduled inbound trucks
    Sch_IT = [];
    Unsch_IT = copy(in_priori);
    #----------------------------------------Start Iteration----------------------------------------------------#
    seq_index = 1;
    #n_outbound
    while cnt <= n_outbound
        Earl_Avail_DOT = 0;
        OT_sel = out_priori[cnt];
        #=
        println("---------------------------------Iteration $(cnt)-------------------------");
        println("$(cnt) -------- $(OT_sel)");
        println("Before Demand fulfillment $(sum(Demand_iter[OT_sel,:]))");
        println("0th $(Stack_Door)")
        println("Free Doors $(Free_Doors)")
        println("0th instance $(Stack_Door[OT_sel])")
        =#
        # fill up outbound truck from temporary storage
        if sum(Demand_iter[OT_sel,:]) > 0 
            Prod_Avail_time = Dict();
            # get goods from all possible scheduled trucks
            for IT in Sch_IT
                if sum([min(Stor_Load[IT, x], Demand_iter[OT_sel,x]) for x = 1:n_products]) > 0 #if scheduled truck has demanded products
                    Prod_Avail_time[IT] = maximum(Entry_IT[IT, :]);
                end;
            end;
            # first clear products from storage area.
            Storage_IT = [];
            if (length(Prod_Avail_time)) > 0
                while length(Prod_Avail_time) > 0
                    Sch_IT_sel = findmin(Prod_Avail_time)[2];
                    #integrated_seq[seq_index] = (Sch_IT_sel,OT_sel);
                    #seq_index += 1;
                    delete!(Prod_Avail_time, Sch_IT_sel);
                    # determine transport of stored products to OT_sel;
                    Load_Trans[Sch_IT_sel,OT_sel,:] = [min(Demand_iter[OT_sel,p],Stor_Load[Sch_IT_sel,p]) for p = 1:n_products];
                    # update demand of OT_sel;
                    if sum(Load_Trans[Sch_IT_sel,OT_sel,:]) == 0
                        continue;
                    end;
                    Demand_iter[OT_sel,:] = Demand_iter[OT_sel,:] - Stor_Load[Sch_IT_sel, :];
                    Demand_iter[Demand_iter .< 0] .= 0;
                    # update storage amount;
                    Stor_Load[Sch_IT_sel, :] = Stor_Load[Sch_IT_sel, :] - Load_Trans[Sch_IT_sel,OT_sel,:];
                    Stor_Load[Stor_Load .< 0] .= 0;
                    push!(Storage_IT, Sch_IT_sel);
                    if sum(Demand_iter[OT_sel,:]) == 0
                        break;
                    end;
                end;
            end;
            #=
            println("The required Storage IT is $(Storage_IT)")
            println("After Demand fulfillment from storage $(sum(Demand_iter[OT_sel,:]))");
            =#
            # select AIT from remaining unscheduled trucks;
            # if there is still demand left to be fulfilled
            AIT = [];
            if sum(Demand_iter[OT_sel,:]) > 0
                for IT_sel in Unsch_IT
                    if sum([min(Demand_iter[OT_sel,p], Supply[IT_sel,p]) for p in 1:n_products]) > 0
                        Stor_Load[IT_sel, :] = Supply[IT_sel, :] - Demand_iter[OT_sel,:]; # Update storage,demand after fulfilling from temporary storage
                        Stor_Load[Stor_Load .< 0] .= 0;
                        Load_Trans[IT_sel, OT_sel, :] = [min(Supply[IT_sel, p], Demand_iter[OT_sel,p]) for p in 1:n_products];
                        Demand_iter[OT_sel,:] = Demand_iter[OT_sel,:] - Supply[IT_sel, :];   # demand left after supply from inbound of AIT
                        Demand_iter[Demand_iter .< 0] .= 0;
                        push!(AIT,IT_sel);
                    end;
                    # fill AIT until demand is fulfilled
                    if sum(Demand_iter[OT_sel,:]) == 0
                        break;
                    end;
                end;
            end;
            #=    
            println("Final Demand fulfillment $(sum(Demand_iter[OT_sel,:]))");
            println("The required AIT is $(AIT)")
            =#
            n_count = Int(ceil(length(AIT)/n_strip));
            Sub_AIT = Dict();
            # divide AIT into groups of n_strip;
            for ind in 1:n_count
                Sub_AIT[ind] = AIT[(ind-1)*n_strip+1:min(ind*n_strip,length(AIT))];
            end;
            entry_time = Array{Float64}(undef,0);
            if length(AIT) > 0
                for ind in 1:n_count
                    IT_num = 1;
                    place_hold = Array{Int8}(undef,0);
                    # assign strip doors and forklifts to inbound trucks;
                    while IT_num <= length(Sub_AIT[ind]) 
                        # assign door and forklift to each inbound truck one at a time
                        IT_sel = Sub_AIT[ind][IT_num];
                        #println("---IT selected Outer ---------$(IT_sel), $(IT_num)")
                        Strip_Door = findmin(Earl_Avail_Strip)[2];
                        req_Strip_forklift = Forklift_finder(Earl_Avail_Strip,Strip_Door,Avail_fork);
                        Entry_IT[IT_sel,Strip_Door] = max(Arrival[IT_sel], Earl_Avail_Strip[Strip_Door], Avail_fork[req_Strip_forklift]);
                        Inf_Entr[IT_sel,req_Strip_forklift] = Entry_IT[IT_sel,Strip_Door];
                        Leave_IT[IT_sel,Strip_Door] = Entry_IT[IT_sel,Strip_Door] + sum(Supply[IT_sel, :].*UnLoadTime);
                        Inf_Leav[IT_sel,req_Strip_forklift] = Leave_IT[IT_sel,Strip_Door];
                        Avail_fork[req_Strip_forklift] = Inf_Leav[IT_sel,req_Strip_forklift]+C_time;
                        Earl_Avail_Strip[Strip_Door] = Leave_IT[IT_sel,Strip_Door]+Change_time;
                        integrated_seq[seq_index] = (IT_sel,0);
                        seq_index += 1;
                        deleteat!(Unsch_IT,findall(x -> x == IT_sel,Unsch_IT));
                        push!(Sch_IT,IT_sel);
                        # track arrival of the next truck
                        if IT_num < length(Sub_AIT[ind])
                            next_IT_sel = Sub_AIT[ind][IT_num+1];
                        else
                            if ind < n_count
                                next_IT_sel = Sub_AIT[ind+1][1];
                            else
                                IT_num += 1;
                                continue;
                            end;
                        end;
                        Strip_Door = findmin(Earl_Avail_Strip)[2];
                        # a placeholder forklift for the next inbound truck
                        req_next_forklift = Forklift_finder(Earl_Avail_Strip,Strip_Door,Avail_fork);
                        next_arrival = max(Arrival[next_IT_sel], Earl_Avail_Strip[Strip_Door], Avail_fork[req_next_forklift]);
                        #println("---next arrival ---------$(next_IT_sel),$(next_arrival)")
                        #curr_load = maximum(Inf_Entr[IT_sel,:])+sum(Supply[IT_sel,:].*UnLoadTime/2);
                        if ind == 1
                            if length(Storage_IT) > 0
                                stor_cnt = 1;
                                while length(Storage_IT) > 0
                                    sel_IT = Storage_IT[stor_cnt];
                                    stor_load = max(maximum(Inf_Entr[sel_IT,:])+sum(Supply[sel_IT,:].*UnLoadTime/2),minimum(Avail_fork));
                                    #=
                                    temp_forklift_avail = Dict();
                                    for i in 1:n_forklifts
                                        if stor_load >= Avail_fork[i]
                                            temp_forklift_avail[i] = stor_load-Avail_fork[i];
                                        end;
                                    end;
                                    =#
                                    # find the earliest forklift for the task
                                    place_forklift = temp_place_forklift(stor_load,Avail_fork);
                                    if Stack_Door[OT_sel] == 0
                                        stor_load = max(stor_load,Avail_fork[place_forklift]);
                                    else
                                        stor_load = max(stor_load,Earl_Avail_Stack[Stack_Door[OT_sel]],Avail_fork[place_forklift]);
                                    end;
                                    stor_load = max(stor_load,Avail_fork[place_forklift]);
                                    if stor_load < next_arrival 

                                        if Stack_Door[OT_sel] == 0            
                                            Stack_Door[OT_sel] = Sub_Stack_Door_finder(Earl_Avail_Stack,Avail_fork,Free_Doors,place_forklift);
                                            stor_load = max(stor_load,Earl_Avail_Stack[Stack_Door[OT_sel]]);
                                            if stor_load > next_arrival
                                                Stack_Door[OT_sel] = 0;
                                                break;
                                            else
                                                push!(Used_Doors, Stack_Door[OT_sel]);
                                                deleteat!(Free_Doors, findall(x -> x == Stack_Door[OT_sel], Free_Doors));
                                            end;
                                        end;

                                        # find forklift after getting the door
                                        #place_forklift = Forklift_finder(Earl_Avail_Stack,Stack_Door,Avail_fork);
                                        Outf_Entr[OT_sel,place_forklift] = max(stor_load,Earl_Avail_Stack[Stack_Door[OT_sel]],Avail_fork[place_forklift]);
                                        Start_Load[sel_IT,OT_sel] = Outf_Entr[OT_sel,place_forklift];
                                        #push!(entry_time,max(stor_load,Earl_Avail_Stack[Stack_Door[OT_sel]],Avail_fork[place_forklift]));
                                        Outf_Leav[OT_sel,place_forklift] = Outf_Entr[OT_sel,place_forklift] +sum(Load_Trans[sel_IT,OT_sel,:].*LoadTime)+TransTime[findmax(Entry_IT[sel_IT,:])[2],Stack_Door[OT_sel]];
                                        push!(Out_Entr[(OT_sel,place_forklift)],Outf_Entr[OT_sel,place_forklift]);
                                        push!(Out_Leav[(OT_sel,place_forklift)],Outf_Leav[OT_sel,place_forklift]);
                                        push!(AIT_time[OT_sel],Outf_Leav[OT_sel,place_forklift]);
                                        # modification do not allow current door in present iteration

                                        # 1st instance of door allocation before laoding barrage 
                                        push!(Used_Doors,Stack_Door[OT_sel]);
                                        deleteat!(Free_Doors, findall(x -> x == Stack_Door[OT_sel], Free_Doors));
                                        #=
                                        if Stack_Door[OT_sel] != 0
                                            Earl_Avail_Stack[Stack_Door[OT_sel]] = 1e6;
                                        end;
                                        =#    
                                        Avail_fork[place_forklift] = findmax(Out_Leav[OT_sel, place_forklift])[1]+C_time;
                                        integrated_seq[seq_index] = (sel_IT,OT_sel);
                                        seq_index += 1;
                                        deleteat!(Storage_IT,findall(x -> x == sel_IT,Storage_IT));
                                        # now find if we can start loading from the new inbound truck in the given time
                                    else
                                        # proceed to unloading of next truck
                                        break;
                                    end;
                                end;
                            end;
                        end;
                        #println("1st instance $(Stack_Door[OT_sel])");
                        # start from the beginning of the subset
                        AIT_IT_cnt = 1;
                        # purpose of rem_tally - to check number of loading assignments before next truck arrival
                        rem_tally = 0;
                        # loading jobs from current Sub_AIT[ind] into current OT
                        #AIT_IT_cnt <= length(Sub_AIT[ind]);
                        while length(Sub_AIT[ind]) > 0;
                            AIT_IT_sel = Sub_AIT[ind][AIT_IT_cnt];
                            temp_forklift_avail = Dict();
                            #println("Set -------------$(Sub_AIT[ind])------------,$(rem_tally) ");
                            #println("Select truck -----$(AIT_IT_sel)--------$(sum(Inf_Entr[AIT_IT_sel,:]))");

                            if sum(Inf_Entr[AIT_IT_sel,:]) != 0
                                curr_load = maximum(Inf_Entr[AIT_IT_sel,:])+sum(Supply[AIT_IT_sel,:].*UnLoadTime/2);    
                            else
                                if rem_tally == 0
                                    # if not possible with first truck, then not possible with subsequent trucks
                                    IT_num += 1;
                                    break;
                                else
                                    break;
                                end;
                            end;

                            for i in 1:n_forklifts
                                if max(curr_load,minimum(Avail_fork)) >= Avail_fork[i]
                                    temp_forklift_avail[i] = curr_load-Avail_fork[i];
                                end;
                            end;
                        # find the forklift nearest to product availability- no need to account idle time with docking outbond truck
                            place_forklift = findmin(temp_forklift_avail)[2];
                            if Stack_Door[OT_sel] == 0
                                curr_load = max(curr_load,Avail_fork[place_forklift]);
                            else
                                curr_load = max(curr_load,Avail_fork[place_forklift],Earl_Avail_Stack[Stack_Door[OT_sel]]);
                            end;
                            if curr_load < next_arrival
                                # condition block if about to assign door for the first time.
                                if Stack_Door[OT_sel] == 0
                                    Stack_Door[OT_sel] = Sub_Stack_Door_finder(Earl_Avail_Stack,Avail_fork,Free_Doors,place_forklift);
                                    curr_load = max(curr_load,Earl_Avail_Stack[Stack_Door[OT_sel]]);
                                    if curr_load > next_arrival
                                        Stack_Door[OT_sel] = 0;
                                    else
                                        push!(Used_Doors, Stack_Door[OT_sel]);
                                        deleteat!(Free_Doors, findall(x -> x == Stack_Door[OT_sel], Free_Doors));
                                    end;
                                end;
                                    #=
                                    if curr_load > next_arrival

                                        # rally only counts unloading from previous iterations, doesn't consider current storage loading
                                        # loop stops removing stack door assignment from previous iteration
                                        if rally[OT_sel] == 0
                                            Stack_Door[OT_sel] = 0;
                                        end;
                                    end;
                                    =#

                                if curr_load > next_arrival                
                                    if rem_tally == 0
                                        IT_num += 1;
                                    end;
                                    break;
                                end;

                                req_AIT_forklift = place_forklift;
                    #req_AIT_forklift = Forklift_finder(Earl_Avail_Stack,Stack_Door,Avail_fork);
                                Outf_Entr[OT_sel,req_AIT_forklift] = max(curr_load,Avail_fork[req_AIT_forklift]);
                                Start_Load[AIT_IT_sel,OT_sel] = Outf_Entr[OT_sel,req_AIT_forklift];
                                #push!(entry_time,max(curr_load,Avail_fork[req_AIT_forklift]));
                                Outf_Leav[OT_sel,req_AIT_forklift] = Outf_Entr[OT_sel,req_AIT_forklift] +sum(Load_Trans[AIT_IT_sel,OT_sel,:].*LoadTime)+TransTime[findmax(Entry_IT[AIT_IT_sel,:])[2],Stack_Door[OT_sel]];
                                push!(Out_Entr[(OT_sel,req_AIT_forklift)],Outf_Entr[OT_sel,req_AIT_forklift]);
                                push!(Out_Leav[(OT_sel,req_AIT_forklift)],Outf_Leav[OT_sel,req_AIT_forklift]);
                                ##println("New OT-sel $(OT_sel), $(req_AIT_forklift) $(Outf_Leav[OT_sel,req_AIT_forklift])");
                                ##println("AIT_time $(AIT_time[OT_sel])");
                                push!(AIT_time[OT_sel],Outf_Leav[OT_sel,req_AIT_forklift]);
                                Avail_fork[req_AIT_forklift] = findmax(Out_Leav[OT_sel, req_AIT_forklift])[1]+C_time;
                                integrated_seq[seq_index] = (AIT_IT_sel,OT_sel);

                                # required door modification 

                                # 2nd instance of door allocation before laoding barrage - under AIT 
                                #=
                                push!(Used_Doors,Stack_Door[OT_sel]);
                                deleteat!(Free_Doors, findall(x -> x == Stack_Door[OT_sel]));
                                =#
                                #=
                                if Stack_Door[OT_sel] != 0
                                    Earl_Avail_Stack[Stack_Door[OT_sel]] = 1e6;
                                end;
                                =#
                                seq_index += 1;
                                push!(place_hold,AIT_IT_sel);
                                deleteat!(Sub_AIT[ind],findall(x -> x == AIT_IT_sel,Sub_AIT[ind]));

                                #IT_num += 1;
                                rem_tally += 1;

                                if rem_tally > 1
                                    IT_num = max(IT_num-1,1);
                                    AIT_IT_cnt = max(AIT_IT_cnt-1,1);

                                end;

                            else
                                if rem_tally == 0
                                    IT_num += 1;
                                end;
                                break;
                            end;
                            #AIT_IT_cnt += 1;
                        end;
                        #println("2nd instance $(Stack_Door[OT_sel])")
                        dum_cnt = cnt+1;
                        cut_flag = 0;
                        # pick next outbound truck
                        #println(" Next Outbound trucks assignment ")
                        # place hold indicates availability of current AIT trucks for future trucks
                        # you can start docking outbound trucks iff you can transfer some products from place_hold to the current outbound truck
                        if length(place_hold) > 0 
                            while dum_cnt <= n_outbound
                                
                                dum_OT_sel = out_priori[dum_cnt];
                                #println("Selected dum_cnt $(dum_cnt) - Future OT - $(dum_OT_sel)");
                                if length(Free_Doors) <  1 
                                    break;
                                end;

                                if rem(dum_cnt-1,n_stack) == 0
                                    break;
                                end;
                                # rally variable to check docking from first loading job 
                                # pick previously unloaded inbound trucks in the same subset
                                for dum_sel_IT in Sch_IT
                                    #println("Future OT selected $(dum_OT_sel)");
                                    if sum([min(Demand_iter[dum_OT_sel,p],Stor_Load[dum_sel_IT,p]) for p in 1:n_products]) > load_factor*sum(Stor_Load[dum_sel_IT,:])
                                        new_dock_load = max(maximum(Inf_Entr[dum_sel_IT,:])+sum(Supply[dum_sel_IT,:].*UnLoadTime/2),minimum(Avail_fork));
                                        # forklift availbility with respect to product availability
                                        temp_forklift = temp_place_forklift(new_dock_load,Avail_fork);                                      
                                        # dock(temporary) OT only when dock is unassigned and meets minimum load requirement
                                        if Stack_Door[dum_OT_sel] == 0
                                            new_dock_load = max(new_dock_load,Avail_fork[temp_forklift]);
                                        else
                                            new_dock_load = max(new_dock_load,Avail_fork[temp_forklift],Earl_Avail_Stack[Stack_Door[dum_OT_sel]]);
                                        end;
                                        # find stack door for the rest of the loop
                                        if Stack_Door[dum_OT_sel] == 0         
                                            Stack_Door[dum_OT_sel] = Sub_Stack_Door_finder(Earl_Avail_Stack,Avail_fork,Free_Doors,temp_forklift);
                                            #Entry_OT[dum_OT_sel,Stack_Door[dum_OT_sel]] = Earl_Avail_Stack[Stack_Door[dum_OT_sel]];
                                            new_dock_load = max(Earl_Avail_Stack[Stack_Door[dum_OT_sel]],Avail_fork[temp_forklift],new_dock_load);
                                            if new_dock_load > next_arrival
                                                Stack_Door[dum_OT_sel] = 0;
                                                break;
                                            else
                                                push!(Used_Doors,Stack_Door[dum_OT_sel]);
                                                deleteat!(Free_Doors, findall(x -> x == Stack_Door[dum_OT_sel], Free_Doors));
                                            end;

                                        end;


                                        #println("OT - $(dum_OT_sel) Stack Door - $(Stack_Door[dum_OT_sel])");
                                        # get latest temp loading time from available docking time
                                        #new_dock_load = max(Earl_Avail_Stack[Stack_Door[dum_OT_sel]],Avail_fork[temp_forklift],new_dock_load);
                                        if new_dock_load < next_arrival
                                            # if first loading job not possible, remove docking assignment.
                                            # else keep entry time and docking door
                                            #=
                                            if rally[dum_OT_sel] == 0
                                                Entry_OT[dum_OT_sel,Stack_Door[dum_OT_sel]] = 0;
                                                Stack_Door[dum_OT_sel] = 0;                                  
                                            end;
                                            # to escape from outer loop
                                            # if first scheduled IT not good enough, discontinue the entire outer loop

                                            #cut_flag = 1;
                                            break;
                                            =#
                                            # assign forklift to docked new OT and to check if the next truck is assigned door prior to iteration
                                            #println("OT - $(dum_OT_sel) Stack Door - $(Stack_Door[dum_OT_sel])");
                                            Outf_Entr[dum_OT_sel,temp_forklift] = new_dock_load;
                                            Start_Load[dum_sel_IT,dum_OT_sel] = Outf_Entr[dum_OT_sel,temp_forklift];
                                            push!(Out_Entr[(dum_OT_sel,temp_forklift)],Outf_Entr[dum_OT_sel,temp_forklift]);
                                            # get product transfer to new outbound truck
                                            Load_Trans[dum_sel_IT,dum_OT_sel,:] = [min(Demand_iter[dum_OT_sel,p],Stor_Load[dum_sel_IT,p]) for p = 1:n_products];
                                            # update demand of next truck
                                            Demand_iter[dum_OT_sel,:] = Demand_iter[dum_OT_sel,:] - Stor_Load[dum_sel_IT, :];   # demand left after supply from inbound of AIT
                                            Demand_iter[Demand_iter .< 0] .= 0;
                                            # update storage of a scheduled inbound truck in the current subset
                                            Stor_Load[dum_sel_IT, :] = Stor_Load[dum_sel_IT, :] - Load_Trans[dum_sel_IT,dum_OT_sel,:]; # Update storage,demand after fulfilling from temporary storage
                                            Stor_Load[Stor_Load .< 0] .= 0;
                                            Outf_Leav[dum_OT_sel,temp_forklift] = Outf_Entr[dum_OT_sel,temp_forklift] +sum(Load_Trans[dum_sel_IT,dum_OT_sel,:].*LoadTime)+TransTime[findmax(Entry_IT[dum_sel_IT,:])[2],Stack_Door[dum_OT_sel]];
                                            push!(Out_Leav[(dum_OT_sel,temp_forklift)],Outf_Leav[dum_OT_sel,temp_forklift]);
                                            push!(AIT_time[dum_OT_sel],Outf_Leav[dum_OT_sel,temp_forklift]);
                                            Avail_fork[temp_forklift] = Outf_Leav[dum_OT_sel,temp_forklift]+C_time;
                                            integrated_seq[seq_index] = (dum_sel_IT,dum_OT_sel);

                                            #println("$(seq_index) - ($(dum_sel_IT),$(dum_sel_OT))")
                                            seq_index += 1;
                                            rally[dum_OT_sel] += 1;
                                            #=
                                            if rally[dum_OT_sel] == 1
                                                push!(Used_Doors,Stack_Door[OT_sel]);
                                                deleteat!(Free_Doors, findall(x -> x == Stack_Door[OT_sel]));
                                            end;
                                            =#
                                            #=
                                            if Stack_Door[dum_OT_sel] != 0
                                                Earl_Avail_Stack[dum_OT_sel] = 1e6;
                                            end;
                                            =#
                                        end;
                                    else
                                        continue;
                                        #cut_flag = 1;
                                    end;
                                    #println("Cut flag from inner iteration $(cut_flag)");
                                end;
                                dum_cnt += 1;
                                #println("Advanced leave time - $(AIT_time[dum_OT_sel])");
                                #println("OT - $(dum_OT_sel) Stack Door - $(Stack_Door[dum_OT_sel])");
                                #=
                                if rem(dum_cnt+1,n_stack) == 0
                                    #println("Paul McShane")
                                    break;
                                else
                                    dum_cnt += 1;
                                end;
                                =#
                            end;
                        end;
                            # to tally for the last truck 
                        #IT_num += 1;
                        

                    end;
                    

                    ##println("2nd instance $(Stack_Door)")
                    # at first assignment    
                    if ind == 1
                        # get stack door
                        if Stack_Door[OT_sel] == 0
                            Stack_Door[OT_sel] = Sub_Stack_Door_finder(Earl_Avail_Stack,Avail_fork,Free_Doors,findmin(Avail_fork)[2]);
                            push!(Used_Doors,Stack_Door[OT_sel]);
                            deleteat!(Free_Doors, findall(x -> x == Stack_Door[OT_sel], Free_Doors));
                        end;
                        Earl_Avail_DOT = Earl_Avail_Stack[Stack_Door[OT_sel]];
                        if length(vcat(Storage_IT,Sub_AIT[ind])) == 0
                            continue;
                        end;
                        for IT_sel in vcat(Storage_IT,Sub_AIT[ind])
                            req_AIT_forklift = Forklift_finder(Earl_Avail_Stack,Stack_Door[OT_sel],Avail_fork);
                            Outf_Entr[OT_sel,req_AIT_forklift] = max(maximum(Inf_Entr[IT_sel,:])+sum(Supply[IT_sel,:].*UnLoadTime/2),Earl_Avail_Stack[Stack_Door[OT_sel]],Avail_fork[req_AIT_forklift]);
                            Start_Load[IT_sel,OT_sel] = Outf_Entr[OT_sel,req_AIT_forklift];
                            Outf_Leav[OT_sel,req_AIT_forklift] = Outf_Entr[OT_sel,req_AIT_forklift] +sum(Load_Trans[IT_sel,OT_sel,:].*LoadTime)+TransTime[findmax(Entry_IT[IT_sel,:])[2],Stack_Door[OT_sel]];
                            push!(Out_Entr[(OT_sel,req_AIT_forklift)],Outf_Entr[OT_sel,req_AIT_forklift]);
                            push!(Out_Leav[(OT_sel,req_AIT_forklift)],Outf_Leav[OT_sel,req_AIT_forklift]);
                            push!(AIT_time[OT_sel],Outf_Leav[OT_sel,req_AIT_forklift]);
                            integrated_seq[seq_index] = (IT_sel,OT_sel);
                            seq_index += 1;
                            Avail_fork[req_AIT_forklift] = findmax(Out_Leav[OT_sel, req_AIT_forklift])[1]+C_time;
                            #iter_cnt += 1;
                        end;
                    else
                        if length(Sub_AIT[ind]) == 0
                            continue;
                        end;
                        for IT_sel in Sub_AIT[ind]
                            req_AIT_forklift = Forklift_finder(Earl_Avail_Stack,Stack_Door[OT_sel],Avail_fork);
                            Outf_Entr[OT_sel,req_AIT_forklift] = max(maximum(Inf_Entr[IT_sel,:])+sum(Supply[IT_sel,:].*UnLoadTime/2),Earl_Avail_Stack[Stack_Door[OT_sel]],Avail_fork[req_AIT_forklift]);
                            Start_Load[IT_sel,OT_sel] = Outf_Entr[OT_sel,req_AIT_forklift];
                            Outf_Leav[OT_sel,req_AIT_forklift] = Outf_Entr[OT_sel,req_AIT_forklift] + +sum(Load_Trans[IT_sel,OT_sel,:].*LoadTime)+TransTime[findmax(Entry_IT[IT_sel,:])[2],Stack_Door[OT_sel]];
                            push!(Out_Entr[(OT_sel,req_AIT_forklift)],Outf_Entr[OT_sel,req_AIT_forklift]);
                            push!(Out_Leav[(OT_sel,req_AIT_forklift)],Outf_Leav[OT_sel,req_AIT_forklift]);
                            push!(AIT_time[OT_sel],Outf_Leav[OT_sel,req_AIT_forklift]);
                            integrated_seq[seq_index] = (IT_sel,OT_sel);
                            seq_index += 1;
                            Avail_fork[req_AIT_forklift] = findmax(Out_Leav[OT_sel, req_AIT_forklift])[1]+C_time;
                        end;
                    end;
                    #println("3rd instnace $(Stack_Door[OT_sel])");
                end;
            else
                if Stack_Door[OT_sel] == 0
                    Stack_Door[OT_sel] = Sub_Stack_Door_finder(Earl_Avail_Stack,Avail_fork,Free_Doors,findmin(Avail_fork)[2]);
                    push!(Used_Doors,Stack_Door[OT_sel]);
                    deleteat!(Free_Doors, findall(x -> x == Stack_Door[OT_sel], Free_Doors));
                end;
                Earl_Avail_DOT = Earl_Avail_Stack[Stack_Door[OT_sel]];
                for IT_sel in Storage_IT
                    req_AIT_forklift = Forklift_finder(Earl_Avail_Stack,Stack_Door[OT_sel],Avail_fork);
                    Outf_Entr[OT_sel,req_AIT_forklift] = max(maximum(Inf_Entr[IT_sel,:])+sum(Supply[IT_sel,:].*UnLoadTime/2),Earl_Avail_Stack[Stack_Door[OT_sel]],Avail_fork[req_AIT_forklift]);
                    Start_Load[IT_sel,OT_sel] = Outf_Entr[OT_sel,req_AIT_forklift];
                    #push!(entry_time,Earl_Avail_Stack[Stack_Door],Avail_fork[req_AIT_forklift]);
                    Outf_Leav[OT_sel,req_AIT_forklift] = Outf_Entr[OT_sel,req_AIT_forklift] + +sum(Load_Trans[IT_sel,OT_sel,:].*LoadTime)+TransTime[findmax(Entry_IT[IT_sel,:])[2],Stack_Door[OT_sel]];
                    push!(Out_Entr[(OT_sel,req_AIT_forklift)],Outf_Entr[OT_sel,req_AIT_forklift]);
                    push!(Out_Leav[(OT_sel,req_AIT_forklift)],Outf_Leav[OT_sel,req_AIT_forklift]);
                    push!(AIT_time[OT_sel],Outf_Leav[OT_sel,req_AIT_forklift]);
                    integrated_seq[seq_index] = (IT_sel,OT_sel);
                    seq_index += 1;
                    Avail_fork[req_AIT_forklift] = findmax(Out_Leav[OT_sel, req_AIT_forklift])[1]+C_time;
                end;
            end;
        end;
        #Entry_OT[OT_sel,Stack_Door] = minimum(entry_time);
        #Leave_OT[OT_sel, Stack_Door] = maximum(AIT_time);
        if Entry_OT[OT_sel,Stack_Door[OT_sel]] == 0
            Entry_OT[OT_sel,Stack_Door[OT_sel]] = Earl_Avail_Stack[Stack_Door[OT_sel]];
        end;
        Leave_OT[OT_sel, Stack_Door[OT_sel]] = max(Entry_OT[OT_sel,Stack_Door[OT_sel]]+sum(Demand[OT_sel,:].*LoadTime),maximum(AIT_time[OT_sel]));
        Earl_Avail_Stack[Stack_Door[OT_sel]] = Leave_OT[OT_sel, Stack_Door[OT_sel]]+Change_time;
        cnt += 1;
        push!(Free_Doors, Stack_Door[OT_sel]);
        deleteat!(Used_Doors, findall(x -> x == Stack_Door[OT_sel], Used_Doors));
        #println("Final Door Allocation $(Stack_Door)")
        ##println("Selected Stack Door $(Stack_Door)");
        ##println("$(Earl_Avail_Stack)");
        ##println("Left over trucks is $(Unsch_IT)");
    end;
                
    return Entry_IT,Leave_IT,Entry_OT,Leave_OT,Load_Trans,Inf_Entr,Inf_Leav,Outf_Entr,Outf_Leav,Out_Entr,Out_Leav,integrated_seq,Start_Load;
        
end;

function calc_fitness(in_seq,out_seq)

    Entr_IT5,Leav_IT5,Entr_OT5,Leav_OT5,Load_Tran = Integrated_Decoder(in_seq,out_seq);
    #print(out_seq)
    #obj_value = maximum([sum(Leav_OT[j, :]) for j = 1:n_outbound]);
    obj_value = 0;

    Tardiness = zeros(Float64, n_outbound);
    for j in 1:n_outbound
        #global obj_value;
        Tardiness[j] = max(0, sum(Leav_OT5[j, g] for g = 1:n_stack) - DeadLines[j]);
        if Tardiness[j] <= Breakpt[1]
            obj_value += Multiplier[1]*Penalty[j]*Tardiness[j];
        elseif Tardiness[j] <= Breakpt[2]
            obj_value += Multiplier[2]*Penalty[j]*Tardiness[j] + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]);
        else
            obj_value += Multiplier[3]*Penalty[j]*Tardiness[j] +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]);
        end
    end;
    
    return obj_value;
    
end;

function PBSA_sequence(in_IT, out_OT, max_iterations, mut_temp, mut_temp_grad, population);
    
    best_IT = copy(in_IT);
    best_OT = copy(out_OT);
    curr_IT = copy(in_IT);
    curr_OT = copy(out_OT);
    curr_soln = calc_fitness(in_IT,out_OT);
    best_soln = calc_fitness(in_IT,out_OT);
    temp_curr = mut_temp;
    iter = 1;
    while iter <= max_iterations
        
        
        sol = 1;
        while sol <= population
            neighbor_options_1 = ["Swap", "Do Nothing"];
            neighbor_options_2 = ["Swap", "Insert", "Do Nothing"];
            p1 = [0.9,0.1];
            p2 = [0.45,0.45,0.1];
            d1 = Categorical(p1);
            d2 = Categorical(p2);
            # select any options from neighbor if n_stack = 1 or n_stack > 1
            if n_stack == 1
                selection_OT = neighbor_options_1[rand(d1)];
            else
                selection_OT = neighbor_options_2[rand(d2)];
            end
            
            if n_strip == 1
                selection_IT = neighbor_options_1[rand(d1)];
            else
                selection_IT = neighbor_options_2[rand(d2)];
            end
            
            while (selection_OT == "Do Nothing") &&  (selection_IT == "Do Nothing")
                if n_stack == 1
                    selection_OT = neighbor_options_1[rand(d1)];
                else
                    selection_OT = neighbor_options_2[rand(d2)];
                end

                if n_strip == 1
                    selection_IT = neighbor_options_1[rand(d1)];
                else
                    selection_IT = neighbor_options_2[rand(d2)];
                end
            end;
          
            if selection_OT == "Swap"
                Next_OT = swap(curr_OT);
            elseif selection_OT == "Insert"
                Next_OT = insert(curr_OT);
            else
                Next_OT = curr_OT;
            end;
            
            if selection_IT == "Swap"
                Next_IT = swap(curr_IT);
            elseif selection_IT == "Insert"
                Next_IT = insert(curr_IT);
            else
                Next_IT = curr_IT;
            end;
            # get objective for the current element in population
            next_soln = calc_fitness(Next_IT,Next_OT);
            
            if next_soln < best_soln
                best_soln = next_soln;
                best_IT = Next_IT;
                best_OT = Next_OT;
            end;
            
            delta_sol = next_soln - curr_soln
            
            if delta_sol > 0 && delta_sol < 10*temp_curr
                prob_accpt = exp(-delta_sol/temp_curr);
            else
                prob_accpt = 0;
            end;
            
            if delta_sol < 0 || prob_accpt >= rand(1)[1]
                curr_IT = Next_IT;
                curr_OT = Next_OT;
                curr_soln = next_soln;
            end;
                
            sol += 1;

        end;
        # update best global solution
        temp_curr = mut_temp_grad*temp_curr;
        iter = iter+1;

    end;
    #=
    plot([1:max_iterations], infeasible_iter, title = "Number of Infeasible Solutions", label = "Supplementary");
    xlabel!("Iteration");
    ylabel!("No. of infeasible solutions");
    savefig("Infeasible Convergence");
    =#
    return best_IT,best_OT
end;


function swap_Job(curr_Job,ind1,ind2)
    swap_Job = deepcopy(curr_Job);
    temp = swap_Job[ind1];
    swap_Job[ind1] = swap_Job[ind2];
    swap_Job[ind2] = temp;
    
    return swap_Job;
    
end;

function insert_Job(curr_Job,ind1,ind2)
    insert_Job = deepcopy(curr_Job);
    index = ind1;
    while index != ind2;
        if ind1 < ind2
            insert_Job = swap_Job(insert_Job,index,index+1);
            index += 1;
        else
            insert_Job = swap_Job(insert_Job,index,index-1);
            index -= 1;
        end;
    end;
    
    return insert_Job;
    
end;

function Job_Decoder(Entr_IT0, Leav_IT0, Entr_OT0, Leav_OT0, Load_Tran0, init_sequence)
    Unsch_IT = [i for i = 1:n_inbound];

    Sch_IT = [];

    Unsch_OT = [j for j = 1:n_outbound];

    Sch_OT = [];

    Entr1_IT = copy(Entr_IT0);

    Leav1_IT = copy(Leav_IT0);

    Entr1_OT = copy(Entr_OT0);

    Leav1_OT = copy(Leav_OT0);

    In_Entr = Dict(); # for possible consideration of multiple instances of same forklift to same truck

    In_Leav = Dict();

    Out_Entr = Dict();

    Out_Leav = Dict();

    Inf_Entr = zeros(Float64,n_inbound,n_forklifts); # Latest times available for forklifts

    Inf_Leav = zeros(Float64,n_inbound,n_forklifts);

    Outf_Entr = zeros(Float64,n_outbound,n_forklifts);

    Outf_Leav = zeros(Float64,n_outbound,n_forklifts);

    Earliest_Avail_DT = zeros(Float64, n_forklifts);
    
    Earliest_out_dock = zeros(Float64, n_outbound);
    
    num_loads = zeros(Int32, n_outbound);
    
    Job_sequence = deepcopy(init_sequence);
    
    new_sequence = Dict();
    
    Entry_time = Dict();
    
    for j in 1:n_outbound
        Entry_time[j] = Array{Float64}(undef,0);
    end;
    
    for i in 1:n_outbound
        for j in 1:n_forklifts
            Out_Entr[(i,j)] = Array{Float64}(undef,0);
            Out_Leav[(i,j)] = Array{Float64}(undef,0);
        end;
    end;
    
    function swap_inside(ind1,ind2)
        temp = Job_sequence[ind1];
        Job_sequence[ind1] = Job_sequence[ind2];
        Job_sequence[ind2] = temp;
    end;
    
    function outdoor_seq1(Leav1_OT)

        OT_time = copy(Leav1_OT);

        outdoor_pos = Dict();

        for door in 1:n_stack
            # calculate number of trailers for each door
            num_trail = 0;
            for OT in 1:n_outbound
                if OT_time[OT, door] > 0
                    num_trail += 1;
                end
            end
            # find positions of trailers at door from decreasing order.
            if num_trail > 0
                outdoor_pos[door] = zeros(Int32,num_trail);
                pos = num_trail;
                 while maximum(OT_time[:, door]) > 0
                     max_OT = findmax(OT_time[:, door])[2];
                     # tuples of OT_chrom only mentioned for outbound trucks - (door utilised,position of truck)
                     outdoor_pos[door][pos] = max_OT;
                     pos -= 1;
                     OT_time[max_OT, door] = 0;
                end;
            end;

        end;

        return outdoor_pos;

    end;
    
    function indoor_seq1(Leav1_IT)

        IT_time = copy(Leav1_IT);

        indoor_pos = Dict();

        for door in 1:n_strip
            # calculate number of trailers for each door
            num_trail = 0;
            for IT in 1:n_inbound
                if IT_time[IT, door] > 0
                    num_trail += 1;
                end
            end
            # find positions of trailers at door from decreasing order.
            if num_trail > 0
                indoor_pos[door] = zeros(Int32,num_trail);
                pos = num_trail;
                 while maximum(IT_time[:, door]) > 0
                     max_IT = findmax(IT_time[:, door])[2];
                     # tuples of OT_chrom only mentioned for outbound trucks - (door utilised,position of truck)
                     indoor_pos[door][pos] = max_IT;
                     pos -= 1;
                     IT_time[max_IT, door] = 0;
                end;
            end;

        end;

        return indoor_pos;

    end;
    #=
    for match in sort(collect(keys(init_sequence)))
        next_OT = init_sequence[match][2];
        if next_OT == 0
            continue;
        else
            num_loads[next_OT] += 1;
        end;
    end;
    =#
    for i in 1:n_inbound
        for j in 1:n_outbound
            if sum(Load_Tran0[i,j,:]) > 0
                num_loads[j] += 1;
            end;
        end;
    end;
    #println(num_loads);
 
    
    init_loads = copy(num_loads);
    Earliest_Avail_DT = zeros(Float64, n_forklifts);
    req_outdoor_seq = outdoor_seq1(Leav1_OT);
    req_indoor_seq = indoor_seq1(Leav1_IT);
    #=
    for match in sort(collect(keys(req_outdoor_seq)))
        println("$(match) - $(req_outdoor_seq[match])")
    end;
    =#
    #print(req_outdoor_seq)
    tally = [];
    search_index = 1;
    travel_index = 1;
    temp_leave = zeros(Float64, n_outbound);
    #=
    check_list = [];
    for ind in sort(collect(keys(init_sequence)))
        push!(check_list,init_sequence[ind]);
    end;
    =#
    flag = 0; 
    while travel_index <= length(Job_sequence)
        flag = 0;

        #Earliest_Avail_DT = DOT(Inf_Leav,Outf_Leav);
        next_IT = Job_sequence[travel_index][1];
        next_OT = Job_sequence[travel_index][2];
        if next_OT == 0
            req_door_in = findmax(Leav1_IT[next_IT,:])[2];
            req_ind_in = findall(x -> x == next_IT,req_indoor_seq[req_door_in])[1];
            
            for ind2 in 1:req_ind_in
                prev_IT = req_indoor_seq[req_door_in][ind2];
                if sum(Inf_Leav[prev_IT,:]) == 0
                    req_index = findall(x -> x == (prev_IT,0),Job_sequence)[1];
                    while req_index > travel_index
                        swap_inside(req_index,req_index-1);
                        req_index -= 1;
                    end;
                    break;
                end;
            end;
            next_IT = Job_sequence[travel_index][1];
            next_OT = Job_sequence[travel_index][2];
            
        else
            if sum(Inf_Leav[next_IT,:]) == 0
                req_door_in = findmax(Leav1_IT[next_IT,:])[2];
                req_ind_in = findall(x -> x == next_IT,req_indoor_seq[req_door_in])[1];
                
                for ind2 in 1:req_ind_in
                    prev_IT = req_indoor_seq[req_door_in][ind2];
                    if sum(Inf_Leav[prev_IT,:]) == 0
                        req_index = findall(x -> x == (prev_IT,0),Job_sequence)[1];
                        while req_index > travel_index
                            swap_inside(req_index,req_index-1);
                            req_index -= 1;
                        end;
                        break;
                    end;
                end;
                next_IT = Job_sequence[travel_index][1];
                next_OT = Job_sequence[travel_index][2];
                
            else
                req_door_out = findmax(Leav1_OT[next_OT,:])[2];
                # find position of truck at the door
                req_ind = findall(x -> x == next_OT,req_outdoor_seq[req_door_out])[1];
                prev_OT = 0;
                prev_loads = 0;
                # previous OT doesn't exist for first outbound truck
                # try finding the sum of previous loading jobs;
                # create another function to calculate sum of leftover loading jobs;
                if req_ind != 1
                    prev_OT = req_outdoor_seq[req_door_out][req_ind-1];
                    prev_loads = num_loads[prev_OT];
                end;
                
                if prev_OT != 0 && prev_loads != 0 
                    # placeholder index to search for other problematic elements
                    place_ind = travel_index+1;
                    for ind in place_ind:length(Job_sequence)
                        # OT in the next element
                        place_IT = Job_sequence[ind][1];
                        place_OT = Job_sequence[ind][2];
                        search_index = ind;
                        # loading job before unloading job problematic, or came across an unloading job
                        
                        if place_OT != 0
                            if sum(Inf_Leav[place_IT,:]) == 0
                                continue;
                            end;
                        end;
                        
                        req_place_ind = 1;
                        # look at the search index over here
                        # find door of the OT in the next element
                        # if unloading job break out of the loop
                        if place_OT != 0
                            place_door_out = findmax(Leav1_OT[place_OT,:])[2];
                            # find position of the OT at the door
                            req_place_ind = findall(x -> x == place_OT, req_outdoor_seq[place_door_out])[1];
                    
                        else
                            flag = 1;
                            break;
                        
                        end;
                        place_prev_OT = 0;
                        # if not the first outbound truck, find the previous truck, else exit the loop
                        if req_place_ind != 1
                            place_prev_OT = req_outdoor_seq[place_door_out][req_place_ind-1];
                        
                        else
                            flag = 1;
                            break;
                        
                        end;
                        # if first truck, exit or if no loads left exit, else continue to find next sequential problematic elements
                        if place_prev_OT == 0 && sum(Inf_Leav[place_IT,:]) != 0
                            flag = 1;
                            break;
                        elseif num_loads[place_prev_OT] == 0 && sum(Inf_Leav[place_IT,:]) != 0
                            flag = 1;
                            break;
                        else
                            continue;
                        end;
                    end;
                    # swap out consecutive elements, to send the first non-problematic element early

                    while search_index > travel_index
                        swap_inside(search_index,search_index-1);
                        search_index -= 1;
                    end;
                    next_IT = Job_sequence[travel_index][1];
                    next_OT = Job_sequence[travel_index][2];                        
                end;
            end;
        end;
        
        if next_OT == 0
            req_door_in = findmax(Leav1_IT[next_IT,:])[2];
            req_ind_in = findall(x -> x == next_IT,req_indoor_seq[req_door_in])[1];
            for ind2 in 1:req_ind_in
                prev_IT = req_indoor_seq[req_door_in][ind2];
                if sum(Inf_Leav[prev_IT,:]) == 0
                    req_index = findall(x -> x == (prev_IT,0),Job_sequence)[1];
                    while req_index > travel_index
                        swap_inside(req_index,req_index-1);
                        req_index -= 1;
                    end;
                    break;
                end;
            end;
        end;
        
        next_IT = Job_sequence[travel_index][1];
        next_OT = Job_sequence[travel_index][2];
        
        #deleteat!(check_list,findall(x -> x == (next_IT,next_OT),check_list));
        # find trucks for the new element
        # do not delete selected job for next iteration
        #deleteat!(priori_val, findall(x -> x == findmin(priori_val)[1], priori_val));
        push!(tally,(next_IT,next_OT));
        #println("$(tally)");
        #println("$(check_list)");
        
        if next_OT == 0
            
            Idle_in = Dict();
            req_door_in = findmax(Leav1_IT[next_IT,:])[2];
            req_ind_in = findall(x -> x == next_IT,req_indoor_seq[req_door_in])[1];
            Earliest_dock = 0;
            if req_ind_in > 1
                prev_IT = req_indoor_seq[req_door_in][req_ind_in-1];
                Earliest_dock = max(Leav1_IT[prev_IT,req_door_in]+Change_time,Arrival[next_IT]);
            else
                Earliest_dock = min(Arrival[next_IT],Entr1_IT[next_IT,req_door_in]);
            end;
            
            Dock_Time = max(Earliest_dock,minimum(Earliest_Avail_DT)); 
            for k in 1:n_forklifts
                if Dock_Time >= Earliest_Avail_DT[k]
                    Idle_in[k] = Dock_Time - Earliest_Avail_DT[k]
                end;
            end;
            
            req_forklift = findmin(Idle_in)[2];    
            Inf_Entr[next_IT,req_forklift] = max(Earliest_dock,Earliest_Avail_DT[req_forklift]);
            Inf_Leav[next_IT,req_forklift] = Inf_Entr[next_IT,req_forklift]+sum(Supply[next_IT,:].*UnLoadTime);
            temp_in_dep = Inf_Leav[next_IT,req_forklift];

            Earliest_Avail_DT[req_forklift] = Inf_Leav[next_IT,req_forklift]+C_time;
            # to consider savings in time from better scheduling
            time_lag_in = maximum(Inf_Leav[next_IT,:])-maximum(Leav1_IT[next_IT,:]);
            for ind in req_ind_in:length(req_indoor_seq[req_door_in])
                curr_IT = req_indoor_seq[req_door_in][ind];
                Entr1_IT[curr_IT,req_door_in] = Entr1_IT[curr_IT,req_door_in]+time_lag_in;
                Leav1_IT[curr_IT,req_door_in] = Leav1_IT[curr_IT,req_door_in]+time_lag_in;
            end;
                    
        else
            # write a similar function to sch_to_chrom for getting a list of sequences at each door;
            # check if the previous trucks are filled out
            # if not shift the job rightwards by a single step
            # shift this code chunk to the front
                
            Idle_out = Dict();
            req_door_out = findmax(Leav1_OT[next_OT,:])[2];
            req_ind_out = findall(x -> x == next_OT,req_outdoor_seq[req_door_out])[1];
            if num_loads[next_OT] == init_loads[next_OT]
                if req_ind_out > 1
                    prev_OT = req_outdoor_seq[req_door_out][req_ind_out-1];
                    Earliest_out_dock[next_OT] = Leav1_OT[prev_OT,req_door_out]+Change_time;
                    #=
                    if Leav1_OT[prev_OT,req_door_out]+Change_time < Entr1_OT[next_OT,req_door_out] 
                        Earliest_out_dock[next_OT] = Leav1_OT[prev_OT,req_door_out]+Change_time;
                    else
                        Earliest_out_dock[next_OT] = Entr1_OT[next_OT,req_door_out];
                    end;
                    =#
                else
                    Earliest_out_dock[next_OT] = Entr1_OT[next_OT,req_door_out];
                end;
                Entr1_OT[next_OT,req_door_out] = Earliest_out_dock[next_OT];
            end;
            
            Out_Time = max(Earliest_out_dock[next_OT],minimum(Earliest_Avail_DT));
            for k in 1:n_forklifts
                if Out_Time >= Earliest_Avail_DT[k]
                    Idle_out[k] = Out_Time - Earliest_Avail_DT[k];
                end;
            end;
            
            req_forklift = findmin(Idle_out)[2];
            Outf_Entr[next_OT,req_forklift] = max(maximum(Inf_Entr[next_IT,:])+sum(Supply[next_IT,:].*UnLoadTime/2),Entr1_OT[next_OT,req_door_out],Earliest_Avail_DT[req_forklift]);
            Outf_Leav[next_OT,req_forklift] = Outf_Entr[next_OT,req_forklift]+sum(Load_Tran0[next_IT,next_OT,:].*LoadTime);
            push!(Out_Entr[(next_OT,req_forklift)],Outf_Entr[next_OT,req_forklift]);
            push!(Out_Leav[(next_OT,req_forklift)],Outf_Leav[next_OT,req_forklift]);
            push!(Entry_time[next_OT],Outf_Leav[next_OT,req_forklift]);
            time_out_dep = Outf_Leav[next_OT,req_forklift];
            Earliest_Avail_DT[req_forklift] = Outf_Leav[next_OT,req_forklift]+C_time;
            num_loads[next_OT] = num_loads[next_OT]-1;
            temp_dur = sum(Demand[next_OT,:].*LoadTime);
            if num_loads[next_OT] == 0
                time_lag_out = max(Earliest_out_dock[next_OT]+temp_dur,maximum(Entry_time[next_OT]))-maximum(Leav1_OT[next_OT,:]);
                Leav1_OT[next_OT,req_door_out] = max(Earliest_out_dock[next_OT]+temp_dur,maximum(Entry_time[next_OT]));
                
                if req_ind_out < length(req_outdoor_seq[req_door_out])
                    for ind in req_ind_out+1:length(req_outdoor_seq[req_door_out])
                        curr_OT = req_outdoor_seq[req_door_out][ind];
                        Entr1_OT[curr_OT,req_door_out] = Entr1_OT[curr_OT,req_door_out]+time_lag_out;
                        Leav1_OT[curr_OT,req_door_out] = Leav1_OT[curr_OT,req_door_out]+time_lag_out;
                    end;
                end;
                
            end;
        end;
        #=
        if next_OT != 0
            println("OT - $(next_OT); IT - $(next_IT) =>$(num_loads[next_OT]); $(sum(Load_Tran0[next_IT,next_OT,:]))");
        end;
        =#
        travel_index += 1;
        search_index += 1;

    end;
    

    return Inf_Entr, Inf_Leav, Outf_Entr, Outf_Leav, Out_Entr, Out_Leav, Entr1_IT, Leav1_IT, Entr1_OT, Leav1_OT, Job_sequence;
    
end;

function final_objective(Entr_IT5, Leav_IT5, Entr_OT5, Leav_OT5)
    
    #obj3_value = maximum([sum(Leav_OT5[j, :]) for j = 1:n_outbound]);
    obj3_value = 0;
    Tardiness = zeros(Float64, n_outbound);
    for j in 1:n_outbound
        #global obj_value;
        Tardiness[j] = max(0, sum(Leav_OT5[j, g] for g = 1:n_stack) - DeadLines[j]);
        if Tardiness[j] <= Breakpt[1]
            obj3_value += Multiplier[1]*Penalty[j]*Tardiness[j];
        elseif Tardiness[j] <= Breakpt[2]
            obj3_value += Multiplier[2]*Penalty[j]*Tardiness[j] + Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]);
        else
            obj3_value += Multiplier[3]*Penalty[j]*Tardiness[j] +  Breakpt[1]*Penalty[j]*(Multiplier[1]-Multiplier[2]) +  Breakpt[2]*Penalty[j]*(Multiplier[2]-Multiplier[3]);
        end
    end;
    
    return obj3_value;
end;

function Job_SA(Entr_IT,Leav_IT,Entr_OT,Leav_OT,Load_Tran,init_seq,init_obj,num_iterations,T_curr,num_pop,T_grad)

    curr_obj = init_obj;
    best_obj = init_obj;
    curr_seq = deepcopy(init_seq);
    curr1_seq = deepcopy(init_seq);
    best_seq = deepcopy(init_seq);
    count = 1;
    prelim_list = [];
    while count <= num_iterations
        neighbor_iter = 1;
        while neighbor_iter <= num_pop
            ind1 = rand(1:length(init_seq));
            ind2 = rand(1:length(init_seq));
            while ind1 == ind2
                ind2 = rand(1:length(init_seq));
            end;
            choice_op = rand(1:2)
            if choice_op == 1
                curr1_seq = swap_Job(curr_seq,ind1,ind2);
            else
                curr1_seq = insert_Job(curr_seq,ind1,ind2);
            end;
            
            Inf_Entr2, Inf_Leav2, Outf_Entr2, Outf_Leav2, Out_Entr2, Out_Leav2,Entr_IT4,Leav_IT4,Entr_OT4,Leav_OT4 = Job_Decoder(Entr_IT,Leav_IT,Entr_OT,Leav_OT,Load_Tran,curr1_seq);
            curr1_obj = final_objective(Entr_IT4,Leav_IT4,Entr_OT4,Leav_OT4);
            delta_sol = curr1_obj - curr_obj;
            if delta_sol > 0 && delta_sol < 10*T_curr
                prob_accept = exp(-delta_sol/T_curr);
            else
                prob_accept = 0;
            end;
            if delta_sol < 0 || prob_accept >= rand(1)[1]
                curr_seq = curr1_seq;
                curr_obj = curr1_obj;
            end;
            if curr1_obj < best_obj
                best_seq = curr1_seq;
                best_obj = curr1_obj;
            end;
            neighbor_iter += 1;
        end;
        push!(prelim_list,best_obj);
        T_curr = T_curr*T_grad;
        count += 1;
    end;
    Inf_Entr3, Inf_Leav3, Outf_Entr3, Outf_Leav3, Out_Entr3, Out_Leav3,Entr_IT4,Leav_IT4,Entr_OT4,Leav_OT4 = Job_Decoder(Entr_IT,Leav_IT,Entr_OT,Leav_OT,Load_Tran,best_seq);
    best_obj = final_objective(Entr_IT4,Leav_IT4,Entr_OT4,Leav_OT4);
    
    return Entr_IT4,Leav_IT4,Entr_OT4,Leav_OT4,best_obj,best_seq,prelim_list;
end;
    
    
    






