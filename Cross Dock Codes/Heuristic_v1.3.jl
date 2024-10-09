function Heuristic_CDS(Inb_β1, Inb_β2, Out_β3, Out_β4, n_strip, n_stack, TransTime)

    Unsch_IT = [i for i = 1:n_inbound];

    Sch_IT = [];

    Unsch_OT = [j for j = 1:n_outbound];

    Sch_OT = [];

    Entr_IT = zeros(Float64, n_inbound, n_strip);

    Leav_IT = zeros(Float64, n_inbound, n_strip);

    Entr_OT = zeros(Float64, n_outbound, n_stack);

    Leav_OT = zeros(Float64, n_outbound, n_stack);

    Temp_storage = zeros(Int32, n_inbound, n_products);

    Load_Tran = zeros(Int32, n_inbound, n_outbound, n_products);
    
    c = 1;

    function MCDM1(Obj1, Obj2, β1, β2)

        best_sol1 = findmin(Obj1)[1];
        best_sol2 = findmin(Obj2)[1];

        Dev1 = Dict();
        Dev2 = Dict();
        Weight_Dev = Dict();

        for i in keys(Obj1)
            if best_sol1 == 0
                Dev1[i] = Obj1[i];
            else
                Dev1[i] = (Obj1[i] - best_sol1)/best_sol1;
            end
        end;

        for i in keys(Obj1)
            if best_sol2 == 0
                Dev2[i] = Obj2[i];
            else
                Dev2[i] = (Obj2[i] - best_sol2)/best_sol2;
            end
        end;

        for i in keys(Obj1)
            Weight_Dev[i] = β1*Dev1[i] + β2*Dev2[i];
        end;

        sel_Trailer = findmin(Weight_Dev)[2];

        return sel_Trailer;

    end;
    function MCDM2(Obj1, Obj2, Obj3, β3, β4, b5)

        best_sol1 = findmin(Obj1)[1];
        best_sol2 = findmin(Obj2)[1];
        best_sol3 = findmin(Obj3)[1];

        Dev1 = Dict();
        Dev2 = Dict();
        Dev3 = Dict();
        Weight_Dev = Dict();

        for i in keys(Obj1)
            if best_sol1 == 0
                Dev1[i] = Obj1[i];
            else
                Dev1[i] = (Obj1[i] - best_sol1)/best_sol1;
            end
        end;

        for i in keys(Obj1)
            if best_sol2 == 0
                Dev2[i] = Obj2[i];
            else
                Dev2[i] = (Obj2[i] - best_sol2)/best_sol2;
            end
        end;
        
        for i in keys(Obj1)
            if best_sol3 == 0
                Dev3[i] = Obj3[i];
            else
                Dev3[i] = (Obj3[i] - best_sol3)/best_sol3;
            end
        end;

        for i in keys(Obj1)
            Weight_Dev[i] = β3*Dev1[i] + β4*Dev2[i] + b5*Dev3[i] ;
        end;

        sel_Trailer = findmin(Weight_Dev)[2];

        return sel_Trailer;

    end;

    function Inb_Strategy1(OT_IS1, OT_list_IS1, IT_list_IS, Dem_IS, Stor_IS1, Supply, LoadTime, Leav_IT, Arrival, Change_time, n_strip, n_stack, n_inbound)

        IT_list_IS1 = copy(IT_list_IS);
        Dem_IS1 = copy(Dem_IS);

        AIT_OT = [];
        Stor_OT = zeros(Int32, n_inbound, n_products);

        Dem_IS1[OT_IS1, :] = Dem_IS1[OT_IS1, :] - sum(Stor_IS1[IT, :] for IT in 1:n_inbound); #storage will be used first to satisfy demand if possible
        Dem_IS1[Dem_IS1 .< 0] .= 0; #change negative values to zero

        Earliest_Avail_DT = zeros(Float64, n_strip); #available DT for ITs
        for k in 1:n_strip
            if maximum(Leav_IT[:,k]) == 0
                Earliest_Avail_DT[k] = maximum(Leav_IT[:,k]);
            else
                Earliest_Avail_DT[k] = maximum(Leav_IT[:,k]) + Change_time;
            end;
        end;

        while sum(Dem_IS1[OT_IS1, :]) > 0

            Exc_unl_time = Dict();
            Excess = Dict();
            Earliest_DT = Dict();
            Earliest_Prod_avail = Dict();

            for i in IT_list_IS1

                if sum([min(Dem_IS1[OT_IS1, p], Supply[i,p]) for p in 1:n_products]) > 0

                    Excess[i] = Supply[i,:] - Dem_IS1[OT_IS1, :]; #excess after satisfying demand of OT
                    Excess[i][Excess[i] .< 0] .= 0;

                    Other_OT = [x for x in 1:length(OT_list_IS1) if x != OT_IS1]; #calculating excess after next OTs [min(n_stack, Unsch_OT)] absorb the excess

                    if length(Other_OT) > 0
                        itera = min(n_stack, length(Other_OT));
                        cnt = 1;
                        while cnt <= itera
                            cnt2 = 1;
                            Exc = 0;
                            rem_OT = 0;
                            for j in Other_OT
                                if cnt2 == 1
                                    Exc = Excess[i] - Dem_IS1[j, :];
                                    Exc[Exc .< 0] .= 0;
                                    rem_OT = j;
                                else
                                    Exc2 = Excess[i] - Dem_IS1[j, :];
                                    Exc2[Exc2 .< 0] .= 0;
                                    if sum(Exc2.*LoadTime) < sum(Exc.*LoadTime)
                                        Exc = Exc2;
                                        rem_OT = j;
                                    end;
                                end;
                                cnt2 += 1;
                            end;
                            Excess[i] = Exc;
                            deleteat!(Other_OT,  findall(x -> x == rem_OT, Other_OT));
                            cnt += 1;
                        end;
                    end;

                    Exc_unl_time[i] = c*sum(Excess[i].*LoadTime); #unloading time for excess supply

                    Earliest_DT[i] = max(Arrival[i], minimum(Earliest_Avail_DT)); #earliest possible DT for the IT

                    Earliest_Prod_avail[i] = Earliest_DT[i] + sum(Supply[i, :].*UnLoadTime)/2;
                end;

            end;

            next_IT = MCDM1(Earliest_Prod_avail, Exc_unl_time, Inb_β1, Inb_β2);

            push!(AIT_OT, next_IT);

            Leave_next_IT = max(Arrival[next_IT], minimum(Earliest_Avail_DT)) + sum(Supply[next_IT, :].*UnLoadTime);
            Idle_time_DD = Dict();
            for k in 1:n_strip
                if Earliest_DT[next_IT] >= Earliest_Avail_DT[k]
                    Idle_time_DD[k] = Earliest_DT[next_IT] - Earliest_Avail_DT[k]
                end;
            end;
            Earliest_Avail_DT[findmin(Idle_time_DD)[2]] = Leave_next_IT + Change_time; #update the earliest available DT after scheduling an IT

            Dem_IS1[OT_IS1, :] = Dem_IS1[OT_IS1, :] - Supply[next_IT,:]; #update demand
            Dem_IS1[Dem_IS1 .< 0] .= 0;

            Stor_OT[next_IT, :] = Excess[next_IT];

            deleteat!(IT_list_IS1, findall(x -> x == next_IT, IT_list_IS1)); #remove the IT added to AIT

        end;

        return AIT_OT, Stor_OT;

    end;

    function Out_Strategy1(OT_list_OS1, AIT_OS1, Stor_OS1, Temp_Stor_OS1, Leav_OT, Leav_IT, n_stack, Demand, Arrival, Supply, TransTime, LoadTime, n_products, DeadLines, Sch_IT, Entr_IT)

        Earl_DT = minimum([maximum(Leav_OT[:, g]) for g in 1:n_stack]);
        if Earl_DT > 0
            Earl_DT = Earl_DT + Change_time;
        end;

        Departure = Dict();
        Dep_deadl = Dict();
        Unload_time = Dict();
        Storage_IT = Dict(); #ITs used for storage

        for OT in keys(AIT_OS1)

            Storage_IT[OT] = [];

            #unloading time for all the AITs of an OT
            if length(AIT_OS1[OT]) > 0
                Unload_time[OT] = sum(sum(Stor_OS1[OT, i, :].*LoadTime) for i in AIT_OS1[OT]);
            else
                Unload_time[OT] = 0;
            end;

            #Entry time of the OT
            Dep_OT = Earl_DT + sum(Demand[OT, :].*LoadTime);
            OT_DD = findmin([maximum(Leav_OT[:, g]) for g in 1:n_stack])[2]; #earliest available stack door

            #earliest available DD for AITs
            Earliest_Avail_DT = zeros(Float64, n_strip); #available DT for ITs
            for k in 1:n_strip
                if maximum(Leav_IT[:,k]) == 0
                    Earliest_Avail_DT[k] = maximum(Leav_IT[:,k]);
                else
                    Earliest_Avail_DT[k] = maximum(Leav_IT[:,k]) + Change_time;
                end;
            end

            #Earliest departures of OT relative to AITs
            Leave_AIT = Dict();

            if length(AIT_OS1[OT]) > 0
                for IT in AIT_OS1[OT]
                    IT_DD = findmin(Earliest_Avail_DT)[2]; #earliest available strip door
                    Leave_AIT[IT] = max(Arrival[IT], minimum(Earliest_Avail_DT)) + sum(Supply[IT, :].*UnLoadTime)/2 + TransTime[IT_DD, OT_DD] + sum(min(Supply[IT, p],Demand[OT, p]).*LoadTime[p] for p = 1:n_products);
                    Earliest_Avail_DT[IT_DD] = max(Arrival[IT], minimum(Earliest_Avail_DT))+sum(Supply[IT, :].*UnLoadTime)+ Change_time;
                end;
            else
                Leave_AIT["None"] = 0;
            end;

            # Create a list of all scheduled ITs used for storage: Earliest avaialble load is used first
            Prod_Avail_time = Dict();
            Demand_copy = copy(Demand[OT, :]);
            for IT in Sch_IT
                if sum([min(Temp_Stor_OS1[IT, x], Demand_copy[x]) for x = 1:n_products]) > 0 #if scheduled truck has demanded products
                    Prod_Avail_time[IT] = maximum(Entr_IT[IT, :]) + sum(Supply[IT,:].*UnLoadTime)/2 + TransTime[findmax(Entr_IT[IT, :])[2] ,OT_DD];
                end;
            end;

            #Earliest possible departure of the OT
            if (length(Prod_Avail_time)) > 0

                while (length(Prod_Avail_time) > 0)
                    Sch_IT_sel = findmin(Prod_Avail_time)[2];
                    delete!(Prod_Avail_time, Sch_IT_sel);
                    Demand_copy = Demand_copy - Temp_Stor_OS1[Sch_IT_sel, :]
                    Demand_copy[Demand_copy .< 0] .= 0;
                    push!(Storage_IT[OT], Sch_IT_sel);
                    if sum(Demand_copy) == 0
                        break;
                    end;
                end;
                #Calculate departure time of OT relative to the ITs used for storage
                Leave_Stor_IT = Dict();
                for IT in Storage_IT[OT]
                    Leave_Stor_IT[IT] = maximum(Entr_IT[IT, :]) + sum(Supply[IT,:].*UnLoadTime)/2 + TransTime[findmax(Entr_IT[IT, :])[2], OT_DD] + sum(min(Temp_Stor_OS1[IT, p], Demand[OT, p]).*LoadTime[p] for p=1:n_products);
                end;

                Departure[OT] = max(Dep_OT, findmax(Leave_AIT)[1], findmax(Leave_Stor_IT)[1]); #need to account for storage availability later
            else
                Departure[OT] = max(Dep_OT, findmax(Leave_AIT)[1]);
            end;

            if Departure[OT] > DeadLines[OT]
                Dep_deadl[OT] = Departure[OT] - DeadLines[OT];
            else
                Dep_deadl[OT] = 0;
            end;
        end;
        #=
        if iterate == 1
            next_OT = 3;
        elseif iterate == 2
            next_OT =  2;
        elseif iterate == 3
            next_OT = 4;
        else
        =#
        
            if findmax(Dep_deadl)[1] > 0
                next_OT = findmax(Dep_deadl)[2];
                #print(next_OT);
                #print(Dep_deadl);                                   
            else                      
            next_OT = MCDM1(Departure, Unload_time, Out_β3, Out_β4);
            end;
        #end;

        return next_OT, Storage_IT[next_OT];

    end;

    function Stage_1(OT_listS1, IT_listS1, Dem, Stor, Supply, LoadTime, Leav_IT, Arrival, Change_time, n_strip, n_stack, n_inbound)

        AIT_iter = Dict();
        Storage_S1 =  zeros(Int32, n_outbound, n_inbound, n_products);

        for j in 1:length(OT_listS1)
            OT_num = j; #row of demand corresponding to OT
            AIT_OT, Stor_OT = Inb_Strategy1(OT_num, OT_listS1, IT_listS1, Dem, Stor, Supply, LoadTime, Leav_IT, Arrival, Change_time, n_strip, n_stack, n_inbound); # function Inb_StrategyX to create associated inbound trailers and storage
            AIT_iter[OT_listS1[j]] = AIT_OT; #add AIT of the OT to the AIT dict
            Storage_S1[OT_listS1[j], :, :] = Stor_OT; # update storage for the AITs of an OT
        end;

        return AIT_iter, Storage_S1;

    end;

    function Stage_2(OT_listS2, AIT, Stor, Temp_Stor_S2, Leav_OT, Leav_IT, n_stack, Demand, Arrival, Supply, TransTime, LoadTime, n_products, DeadLines, Sch_IT, Entr_IT)

        OT_iter, Stor_IT_list = Out_Strategy1(OT_listS2, AIT, Stor, Temp_Stor_S2, Leav_OT, Leav_IT, n_stack, Demand, Arrival, Supply, TransTime, LoadTime, n_products, DeadLines, Sch_IT, Entr_IT);
                              
        return OT_iter, Stor_IT_list; #selectd OT based on OT Strategy

    end;

    iteration = 1;
    In_out = Dict();
    # Outer iteration, runs when any outbound trucks are unscheduled                          
    while length(Unsch_OT) > 0

        #println("---------------------------Heuristic Iteration: $(iteration)---------------------------")

        if length(Unsch_IT) > 0

            Demand_iter = copy(Demand[Unsch_OT, :]);

            AIT, Storage_iter = Stage_1(Unsch_OT, Unsch_IT, Demand_iter, Temp_storage, Supply, LoadTime, Leav_IT, Arrival, Change_time, n_strip, n_stack, n_inbound);
# provides AIT and Storage iter of each unsceduled OT 
        else

            AIT = Dict();
            for OT in Unsch_OT
                AIT[OT] = [];
            end;

            Storage_iter = 0;

        end;

        OT_sel, Stor_IT_Used = Stage_2(Unsch_OT, AIT, Storage_iter, Temp_storage, Leav_OT, Leav_IT, n_stack, Demand, Arrival, Supply, TransTime, LoadTime, n_products, DeadLines, Sch_IT, Entr_IT);
# provides the selected scheduled OT and storage IT's used .
        #=
        println("*Stage 1 => ")
       
        for OT in sort(collect(keys(AIT)))
            if length(AIT[OT]) > 0
                println("\tOT$(OT): ", ["IT$(x)" for x in AIT[OT]]);
            else
                println("\tOT$(OT): Empty Set");
            end;
        end;

        #println("*Stage 2 => ");

        println("\tOT_$(OT_sel) selected");
        =#
        
        #update OT entry time
        Earl_Avail_DOT = minimum([maximum(Leav_OT[:,g]) for g in 1:n_stack]); #earliest available docking time at stack door,maximum gives the latest time for an empty door
        Earl_Avail_Stack = argmin([maximum(Leav_OT[:,g]) for g in 1:n_stack]); #earliest available stack door

        temp_Arr_OT = 0;
        if Earl_Avail_DOT == 0
            temp_Arr_OT = Earl_Avail_DOT; #Entr_OT[OT_sel, Earl_Avail_Stack] = Earl_Avail_DOT;
        else
            temp_Arr_OT = Earl_Avail_DOT + Change_time; #Entr_OT[OT_sel, Earl_Avail_Stack] = Earl_Avail_DOT + Change_time;
        end;
        #update OT leave time 
        OT_min_load = sum(Demand[OT_sel, :].*LoadTime);

        #decrease storage caused by usage by the OT_sel
        if length(Stor_IT_Used) > 0
            #arr = Stor_IT_Used;
            len_iter = 1;
            Stor_tran_load = zeros(Float64, length(Stor_IT_Used));  # 
            dem_copy2 = copy(Demand[OT_sel, :]); 

            #println("\tStorage ITs Used: ", ["IT$(x)" for x in Stor_IT_Used]);
        # 
            for IT in Stor_IT_Used

                Stor_transfer = [min(Temp_storage[IT, x], dem_copy2[x]) for x in 1:n_products]; # minimum of demand of seleccted OT/temporary storage from Stor_IT_Used giving amount of transport from storage
                Load_Tran[IT, OT_sel, :] = Stor_transfer;  #transport between a stor_IT and selected OT
                Temp_storage[IT, :] = Temp_storage[IT, :] - Stor_transfer;  # temporary storage for each stor_IT reduced by transport from storage
                dem_copy2 = dem_copy2 - Stor_transfer;  # Demand for selected OT reduced by transport from Stor_IT to selected OT 
                dem_copy2[dem_copy2 .< 0] .= 0;

                Stor_tran_load[len_iter] = maximum(Entr_IT[IT, :]) + sum(Supply[IT, :].*UnLoadTime)/2 + TransTime[findmax(Entr_IT[IT, :])[2], Earl_Avail_Stack] + sum(Stor_transfer.*LoadTime);  #time for transport between stor_IT and OT , above expression of supply[] doubtful
                len_iter += 1;
            end;
        else
            #arr = [];
            #println("\tStorage ITs Used: Empty Set")
            Stor_tran_load = 0;
        end;

        #if AIT is non empty
        if length(AIT[OT_sel]) > 0

            #println("\tAIT: ", ["IT$(x)" for x in AIT[OT_sel]]);
            In_out[OT_sel] = AIT[OT_sel];
            # dem_copy provides demand left after utilizing temporary storage  
            if length(Stor_IT_Used) > 0
                dem_copy = copy(dem_copy2);
            else
                dem_copy = copy(Demand[OT_sel, :]);
            end;

            AIT_min_load =  zeros(Float64, length(AIT[OT_sel]));

            for IT_sel in AIT[OT_sel]

                #increase storage caused by AIT of the OT selected
                Temp_storage[IT_sel, :] = Supply[IT_sel, :] - dem_copy; # Update storage,demand after fulfilling from temporary storage
                Temp_storage[Temp_storage .< 0] .= 0;
                Load_Tran[IT_sel, OT_sel, :] = [min(Supply[IT_sel, p], dem_copy[p]) for p in 1:n_products];
                dem_copy = dem_copy - Supply[IT_sel, :];   # demand left after supply from inbound of AIT
                dem_copy[dem_copy .< 0] .= 0;
               
                #Earliest available docking time at each strip door                        
                Earl_DT_Strip = zeros(Float64, n_strip);
                for k in 1:n_strip
                    if maximum(Leav_IT[:, k]) > 0
                        Earl_DT_Strip[k] = maximum(Leav_IT[:, k]) + Change_time;
                    else
                        Earl_DT_Strip[k] = maximum(Leav_IT[:, k]); #or zero
                    end;
                end;

                Earl_Avail_DIT = minimum(Earl_DT_Strip); #earliest available docking time at strip door
                # Calculate idle time at each strip door
                Idle_time_Strip = Dict();
                for k in 1:n_strip
                    if max(Arrival[IT_sel], Earl_Avail_DIT) >= Earl_DT_Strip[k]
                        Idle_time_Strip[k] = max(Arrival[IT_sel], Earl_Avail_DIT) - Earl_DT_Strip[k]; # arrival of truck/earliest availability of any strip door - latest availability of strip door
                    end;
                end;

                Earl_Avail_Strip = findmin(Idle_time_Strip)[2]; #earliest available strip door [2] stands to extract corresponding index

                #IT Entry time
                Entr_IT[IT_sel, Earl_Avail_Strip] = max(Arrival[IT_sel], Earl_Avail_DIT); #latest of arrival of inbound truck/availability of the strip door 

                #IT Leave time
                Leav_IT[IT_sel, Earl_Avail_Strip] = Entr_IT[IT_sel, Earl_Avail_Strip] + sum(Supply[IT_sel, :].*UnLoadTime); #entry time+ total unloading time of the truck = leave time of IT

                push!(AIT_min_load, Entr_IT[IT_sel, Earl_Avail_Strip] + sum(Supply[IT_sel, :].*UnLoadTime)/2 + TransTime[Earl_Avail_Strip, Earl_Avail_Stack] +  sum(min(Supply[IT_sel, p], Demand[OT_sel, p]).*LoadTime[p] for p=1:n_products)); # AIT_min_load refers to entry+ unload+transport+loading time of product form each AIT.

                deleteat!(Unsch_IT, findall(x -> x == IT_sel, Unsch_IT)); # Update Unscheduled IT list

                push!(Sch_IT, IT_sel) # Update Scheduled IT list
            end;
# for each AIT in selected OT
            Leav_OT[OT_sel, Earl_Avail_Stack] = max(temp_Arr_OT + OT_min_load, maximum(AIT_min_load), maximum(Stor_tran_load)); #max(Entr_OT[OT_sel, Earl_Avail_Stack] + OT_min_load, maximum(AIT_min_load), maximum(Stor_tran_load));

        else

            #println("\tAIT: Empty Set")

            Leav_OT[OT_sel, Earl_Avail_Stack] = max(temp_Arr_OT + OT_min_load, maximum(Stor_tran_load));

        end;

        Entr_OT[OT_sel, Earl_Avail_Stack] = Leav_OT[OT_sel, Earl_Avail_Stack] - OT_min_load;

        deleteat!(Unsch_OT, findall(x -> x == OT_sel, Unsch_OT)); # Update Unscheduled OT list

        push!(Sch_OT, OT_sel);  # Update Scheduled OT list

    iteration += 1;
    end;

    return Entr_IT, Leav_IT, Entr_OT, Leav_OT, Load_Tran, In_out;
end;
