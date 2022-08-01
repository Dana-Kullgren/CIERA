#04/12/2022
#File with the functions to get the full formation history for an object

print("Script is now running...")

#Imports 
import numpy as np 
import pandas as pd
from time_conversion import read_units                  ##this is the line I from BH_mergers.py
# from ipynb.fs.full.time_conversion import read_units  ##this was the original line

print('Finished imports...')

#Functions
def check_semerge(path, OUTPUT_ID, t, ID):
    
    '''
    This function checks if there is an instance of a given ID in the semergedisrupt file
    
    PARAMETERS
    -----------------
    path : str
    OUTPUT_ID:  list,str of rv values
    t:  list, time of the interaction
    ID: list, int of number of interations per model
    
    
    
    RETURN
    -----------------
    result : remnant ID of the interaction
    df3 : data frame with the following information:
            BH_ID: flt,ID of the original BH
            Type: str, 'semergedisrupt', indicates that there was a merger
            Time: str, time of the merger
            Mass1: flt, mass of object 1
            Mass2: flt, mass of object 1
            Merger Mass: flt, mass of merger
            Mass_1_IDs: flt, ID of object 1
            Mass_2_IDs: flt, ID of object 2
            Merger_IDs: flt, ID of merger
            Mass1_type: flt, type of object 1
            Mass2_type: flt, type of object 2
            Merger_type: flt, type of merger
            Final_bh_mass: 0, used to combine with other dataframes
            
    '''
    
    f3 = open(path + 'initial.semergedisrupt.log','r')
    lines3= f3.readlines()
    columns= ["BH ID", "Type" ,"Time", "Mass1" , "Mass2" , "MergerMass" , "Mass1_IDs" , "Mass2_IDs" , "Merger_IDs" , "Mass1_type" , "Mass2_type" , "Merger_type" , "Final_bh_mass"]
    df3 = pd.DataFrame(columns=columns)
    result = -100
    
    for j in range(len(lines3)-1,0, -1):
        line3 = lines3[j]
        line3 = line3.split('\n')
        line3 = line3[0]
        data = line3.split(' ')
        type_int = data[1]
        if type_int != 'disruptboth':#disruptboth does not change ids, so we don't have to worry
            time = data[0].split('=')
            time = float(time[1])
            idr     = data[2].split('=')
            r_mass = idr[-1].split(')')
            r_mass = r_mass[0]
            idr     = idr[1].split('(')
            idr = float(idr[0])
            initial = data[3].split('=')
            id1     = initial[1].split('(')
            id1_mass = initial[2].split(')')
            id1_mass = float(id1_mass[0])
            id1     = float(id1[0])
            id2     =  initial[3].split('(')
            id2_mass = initial[4].split(')')
            id2_mass = float(id2_mass[0])
            id2     = float(id2[0])
            
            if t > time and (int(id1) == int(OUTPUT_ID)) or (int(id2) == int(OUTPUT_ID)):
                disrupt_time = str(time)
                disrupt_ID = float(idr)
                disrupt_mass = float(r_mass)
                disrupt_id1 = float(id1)
                disrupt_id2 = float(id2)
                disrupt_id1_type = data[6].split('=')
                disrupt_id1_type = float(disrupt_id1_type[1])
                disrupt_mass1 = float(id1_mass)
                disrupt_id2_type = data[7].split('=')
                disrupt_id2_type = float(disrupt_id2_type[1])
                disrupt_remant_type = data[5].split('=')
                disrupt_mass2 = float(id2_mass)
                disrupt_remant_type = float(disrupt_remant_type[1])
                run = 'semergedisrupt'
                fin_bh_mass = 0.0
                ID = float(ID)
                
                df2 = pd.DataFrame({"Type": [run], "BH ID":[ID], "Time" : disrupt_time, "Mass1" : disrupt_mass1, "Mass2" : disrupt_mass2, "MergerMass" : disrupt_mass, "Mass1_IDs" : disrupt_id1, "Mass2_IDs" : disrupt_id2, "Merger_IDs" : disrupt_ID, "Mass1_type" :disrupt_id1_type, "Mass2_type" : disrupt_id2_type , "Merger_type" : disrupt_remant_type, "Final_bh_mass" : fin_bh_mass})
                print(df2)
                df3 = df3.append(df2,ignore_index = True)[df2.columns.tolist()]

                if disrupt_mass1 > disrupt_mass2: #keep id of the most massive start
                    result = id1
                if disrupt_mass1 < disrupt_mass2:
                    result = id2


    return result, df3

def collision(OUTPUT_ID, stop_num,lines2,mass_bh, ID):
    
    '''
    This function checks if there is an instance of a given ID in the collision.log file
    
    PARAMETERS
    -----------------
    OUTPUT_ID: int, specific ID we are looking at
    stop_num: int, line number from which to begin the search 
    lines2: files
    mass_bh foat, mass of the BH
    ID: ID for the initial body 
    
    
    RETURN
    -----------------
    
    OUTPUT_ID: Last ID in the tree 
    ids_to_check: list of ids to check in the future to complete the formation tree 
    df : data frame with the following information:
            BH_ID: flt,ID of the original BH
            Type: str, collision type
            Time: str, time of the collision
            Mass1: flt, mass of object 1
            Mass2: flt, mass of object 1
            Merger Mass: flt, mass of collision
            Mass_1_IDs: flt, ID of object 1
            Mass_2_IDs: flt, ID of object 2
            Merger_IDs: flt, ID of merger
            Mass1_type: flt, type of object 1
            Mass2_type: flt, type of object 2
            Merger_type: flt, type of collision
            Final_bh_mass: mass of BH, used to combine with other dataframes
            
    '''
    
    #### Now read through output backwards from t=BH_formation_time to t=0
    
    columns= ["BH ID", "Type" ,"Time", "Mass1" , "Mass2" , "MergerMass" , "Mass1_IDs" , "Mass2_IDs" , "Merger_IDs" , "Mass1_type" , "Mass2_type" , "Merger_type" , "Final_bh_mass"]
    df = pd.DataFrame(columns=columns)
    t_array = []
    m_total_array = []
    m_input_array = []
    input_mass_array = []
    input_mass_array2 = []
    input_ID_array = []
    input_ID_array2 = []
    output_ID_array = []
    input_kw_array = []
    input_kw_array2 = []
    output_kw_array = []
    count = 0
    count_array = []
    b_array = []
    v_inf_array = []
    input_array = []
    collision_type = []
    bh_ID = []
    for j in range(stop_num, 0, -1):
        line2 = lines2[j]
        line2 = line2.split('\n')
        line2 = line2[0]
        data = line2.split(' ')
        t = data[0]
        t = t.split('=')
        t = t[1]
        type_int = data[1]

        if type_int == 'single-single':
            max = len(data)
            v_inf = data[max-1]
            v_inf = v_inf.split('=')
            v_inf = float(v_inf[1])
            b = data[max-2]
            b = b.split('=')
            b = float(b[1])
        else:
            v_inf = -1
            b = -1

        output = data[2].split('=')
        output_id = output[1]
        output_id = output_id.split('(')
        output_id = float(output_id[0])

        output_mass = output[2]
        output_mass = output_mass.split(')')
        output_mass = float(output_mass[0])

        output_kw = data[5].split('=')
        output_kw = float(output_kw[1])

        if int(output_id) == int(OUTPUT_ID):
        ### Now find all input ids for this collision and add them to the input_array
            input = data[3].split('=')
            
            if len(input) == 5:
                count_array.append(count)
                b_array.append(b)
                v_inf_array.append(v_inf)
                t_array.append(t)
                m_total_array.append(output_mass)
                collision_type.append(type_int)
                bh_ID.append(ID)
                input_id_1 = input[1].split('(')
                input_id_1 = float(input_id_1[0])
                input_mass_1 = input[2].split(')')
                input_mass_1 = float(input_mass_1[0])
                input_array.append(input_id_1)

                input_id_2 = input[3].split('(')
                input_id_2 = float(input_id_2[0])
                input_mass_2 = input[4].split(')')
                input_mass_2 = float(input_mass_2[0])
                input_array.append(input_id_2)

                input_kw_1 = data[6].split('=')
                input_kw_1 = float(input_kw_1[1])
                input_kw_2 = data[7].split('=')
                input_kw_2 = float(input_kw_2[1])

                
                if input_mass_1 < input_mass_2:
                    input_mass_array.append(input_mass_1)
                    input_mass_array2.append(input_mass_2)
                    OUTPUT_ID = input_id_2
                    input_ID_array.append(input_id_1)
                    input_ID_array2.append(input_id_2)
                    output_ID_array.append(output_id)

                    input_kw_array.append(input_kw_1)
                    input_kw_array2.append(input_kw_2)
                    output_kw_array.append(output_kw)
                else:
                    input_mass_array.append(input_mass_2)
                    input_mass_array2.append(input_mass_1)
                    OUTPUT_ID = input_id_1
                    input_ID_array.append(input_id_2)
                    input_ID_array2.append(input_id_1)
                    output_ID_array.append(output_id)

                    input_kw_array.append(input_kw_2)
                    input_kw_array2.append(input_kw_1)
                    output_kw_array.append(output_kw)

            if len(input) == 7:
                count_array.append(count)
                count_array.append(count)

                b_array.append(b)
                b_array.append(b)

                v_inf_array.append(v_inf)
                v_inf_array.append(v_inf)

                t_array.append(t)
                t_array.append(t)

                m_total_array.append(output_mass)
                m_total_array.append(output_mass)

                collision_type.append(type_int)
                collision_type.append(type_int)
                
                bh_ID.append(ID)
                bh_ID.append(ID)

                id_array = []
                mass_array = []
                input_id_1 = input[1].split('(')
                input_id_1 = float(input_id_1[0])
                input_mass_1 = input[2].split(')')
                input_mass_1 = float(input_mass_1[0])
                id_array.append(input_id_1)
                mass_array.append(input_mass_1)

                input_id_2 = input[3].split('(')
                input_id_2 = float(input_id_2[0])
                input_mass_2 = input[4].split(')')
                input_mass_2 = float(input_mass_2[0])
                id_array.append(input_id_2)
                mass_array.append(input_mass_2)

                input_id_3 = input[5].split('(')
                input_id_3 = float(input_id_3[0])
                input_mass_3 = input[6].split(')')
                input_mass_3 = float(input_mass_3[0])
                id_array.append(input_id_3)
                mass_array.append(input_mass_3)

                input_kw_1 = data[6].split('=')
                input_kw_1 = float(input_kw_1[1])
                input_kw_2 = data[7].split('=')
                input_kw_2 = float(input_kw_2[1])
                input_kw_3 = data[8].split('=')
                input_kw_3 = float(input_kw_3[1])

                if input_mass_1 < input_mass_2 and input_mass_1 < input_mass_3:
                    input_mass_array.append(input_mass_1)
                    input_ID_array.append(input_id_1)
                    input_kw_array.append(input_kw_1)
                    if input_mass_2 > input_mass_3:
                        input_mass_array.append(input_mass_3)
                        input_mass_array2.append(input_mass_2)
                        input_mass_array2.append(input_mass_2)
                        OUTPUT_ID = input_id_2
                        input_ID_array.append(input_id_3)
                        input_ID_array2.append(input_id_2)
                        input_ID_array2.append(input_id_2)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_3)
                        input_kw_array2.append(input_kw_2)
                        input_kw_array2.append(input_kw_2)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    else:
                        input_mass_array.append(input_mass_2)
                        input_mass_array2.append(input_mass_3)
                        input_mass_array2.append(input_mass_3)
                        OUTPUT_ID = input_id_3
                        input_ID_array.append(input_id_2)
                        input_ID_array2.append(input_id_3)
                        input_ID_array2.append(input_id_3)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_2)
                        input_kw_array2.append(input_kw_3)
                        input_kw_array2.append(input_kw_3)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                if input_mass_2 < input_mass_1 and input_mass_2 < input_mass_3:
                    input_mass_array.append(input_mass_2)
                    input_ID_array.append(input_id_2)
                    input_kw_array.append(input_kw_2)
                    if input_mass_3 > input_mass_1:
                        input_mass_array.append(input_mass_1)
                        input_mass_array2.append(input_mass_3)
                        input_mass_array2.append(input_mass_3)
                        OUTPUT_ID = input_id_3
                        input_ID_array.append(input_id_1)
                        input_ID_array2.append(input_id_3)
                        input_ID_array2.append(input_id_3)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)


                        input_kw_array.append(input_kw_1)
                        input_kw_array2.append(input_kw_3)
                        input_kw_array2.append(input_kw_3)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    else:
                        input_mass_array.append(input_mass_3)
                        input_mass_array2.append(input_mass_1)
                        input_mass_array2.append(input_mass_1)
                        OUTPUT_ID = input_id_1
                        input_ID_array.append(input_id_3)
                        input_ID_array2.append(input_id_1)
                        input_ID_array2.append(input_id_1)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_3)
                        input_kw_array2.append(input_kw_1)
                        input_kw_array2.append(input_kw_1)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                if input_mass_3 < input_mass_2 and input_mass_3 < input_mass_1:
                    input_mass_array.append(input_mass_3)
                    input_ID_array.append(input_id_3)
                    input_kw_array.append(input_kw_3)

                    if input_mass_2 > input_mass_1:
                        input_mass_array.append(input_mass_1)
                        input_mass_array2.append(input_mass_2)
                        input_mass_array2.append(input_mass_2)
                        OUTPUT_ID = input_id_2
                        input_ID_array.append(input_id_1)
                        input_ID_array2.append(input_id_2)
                        input_ID_array2.append(input_id_2)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_1)
                        input_kw_array2.append(input_kw_2)
                        input_kw_array2.append(input_kw_2)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    else:
                        input_mass_array.append(input_mass_2)
                        input_mass_array2.append(input_mass_1)
                        input_mass_array2.append(input_mass_1)
                        OUTPUT_ID = input_id_1
                        input_ID_array.append(input_id_2)
                        input_ID_array2.append(input_id_1)
                        input_ID_array2.append(input_id_1)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_2)
                        input_kw_array2.append(input_kw_1)
                        input_kw_array2.append(input_kw_1)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
            if len(input) == 9:
                count_array.append(count)
                count_array.append(count)
                count_array.append(count)

                b_array.append(b)
                b_array.append(b)
                b_array.append(b)


                v_inf_array.append(v_inf)
                v_inf_array.append(v_inf)
                v_inf_array.append(v_inf)


                t_array.append(t)
                t_array.append(t)
                t_array.append(t)

                m_total_array.append(output_mass)
                m_total_array.append(output_mass)
                m_total_array.append(output_mass)
                collision_type.append(type_int)
                collision_type.append(type_int)
                collision_type.append(type_int)
                
                bh_ID.append(ID)
                bh_ID.append(ID)
                bh_ID.append(ID)

                id_array = []
                mass_array = []

                input_id_1 = input[1].split('(')
                input_id_1 = float(input_id_1[0])
                input_mass_1 = input[2].split(')')
                input_mass_1 = float(input_mass_1[0])
                id_array.append(input_id_1)
                mass_array.append(input_mass_1)

                input_id_2 = input[3].split('(')
                input_id_2 = float(input_id_2[0])
                input_mass_2 = input[4].split(')')
                input_mass_2 = float(input_mass_2[0])
                id_array.append(input_id_2)
                mass_array.append(input_mass_2)

                input_id_3 = input[5].split('(')
                input_id_3 = float(input_id_3[0])
                input_mass_3 = input[6].split(')')
                input_mass_3 = float(input_mass_3[0])
                id_array.append(input_id_3)
                mass_array.append(input_mass_3)

                input_id_4 = input[7].split('(')
                input_id_4 = float(input_id_4[0])
                input_mass_4 = input[8].split(')')
                input_mass_4 = float(input_mass_4[0])
                id_array.append(input_id_4)
                mass_array.append(input_mass_4)

                input_kw_1 = data[6].split('=')
                input_kw_1 = float(input_kw_1[1])
                input_kw_2 = data[7].split('=')
                input_kw_2 = float(input_kw_2[1])
                input_kw_3 = data[8].split('=')
                input_kw_3 = float(input_kw_3[1])
                input_kw_4 = data[9].split('=')
                input_kw_4 = float(input_kw_4[1])
                if input_mass_1 < input_mass_2 and input_mass_1 < input_mass_3 and input_mass_1 < input_mass_4:
                    input_mass_array.append(input_mass_1)
                    input_ID_array.append(input_id_1)
                    input_kw_array.append(input_kw_1)
                    if input_mass_2 > input_mass_3 and input_mass_2 > input_mass_4:
                        input_mass_array.append(input_mass_3)
                        input_mass_array.append(input_mass_4)
                        input_mass_array2.append(input_mass_2)
                        input_mass_array2.append(input_mass_2)
                        input_mass_array2.append(input_mass_2)

                        OUTPUT_ID = input_id_2
                        input_ID_array.append(input_id_3)
                        input_ID_array.append(input_id_4)
                        input_ID_array2.append(input_id_2)
                        input_ID_array2.append(input_id_2)
                        input_ID_array2.append(input_id_2)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_3)
                        input_kw_array.append(input_kw_4)
                        input_kw_array2.append(input_kw_2)
                        input_kw_array2.append(input_kw_2)
                        input_kw_array2.append(input_kw_2)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    elif input_mass_4 > input_mass_3 and input_mass_4 > input_mass_2:
                        input_mass_array.append(input_mass_3)
                        input_mass_array.append(input_mass_2)
                        input_mass_array2.append(input_mass_4)
                        input_mass_array2.append(input_mass_4)
                        input_mass_array2.append(input_mass_4)
                        OUTPUT_ID = input_id_4
                        input_ID_array.append(input_id_3)
                        input_ID_array.append(input_id_2)
                        input_ID_array2.append(input_id_4)
                        input_ID_array2.append(input_id_4)
                        input_ID_array2.append(input_id_4)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_3)
                        input_kw_array.append(input_kw_2)
                        input_kw_array2.append(input_kw_4)
                        input_kw_array2.append(input_kw_4)
                        input_kw_array2.append(input_kw_4)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    else:
                        input_mass_array.append(input_mass_2)
                        input_mass_array.append(input_mass_4)
                        input_mass_array2.append(input_mass_3)
                        input_mass_array2.append(input_mass_3)
                        input_mass_array2.append(input_mass_3)
                        OUTPUT_ID = input_id_3
                        input_ID_array.append(input_id_2)
                        input_ID_array.append(input_id_4)
                        input_ID_array2.append(input_id_3)
                        input_ID_array2.append(input_id_3)
                        input_ID_array2.append(input_id_3)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_2)
                        input_kw_array.append(input_kw_4)
                        input_kw_array2.append(input_kw_3)
                        input_kw_array2.append(input_kw_3)
                        input_kw_array2.append(input_kw_3)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                if input_mass_2 < input_mass_1 and input_mass_2 < input_mass_3 and input_mass_2 < input_mass_4:
                    input_mass_array.append(input_mass_2)
                    input_ID_array.append(input_id_2)
                    input_kw_array.append(input_kw_2)
                    if input_mass_1 > input_mass_3 and input_mass_1 > input_mass_4:
                        input_mass_array.append(input_mass_3)
                        input_mass_array.append(input_mass_4)
                        input_mass_array2.append(input_mass_1)
                        input_mass_array2.append(input_mass_1)
                        input_mass_array2.append(input_mass_1)
                        OUTPUT_ID = input_id_1
                        input_ID_array.append(input_id_3)
                        input_ID_array.append(input_id_4)
                        input_ID_array2.append(input_id_1)
                        input_ID_array2.append(input_id_1)
                        input_ID_array2.append(input_id_1)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_3)
                        input_kw_array.append(input_kw_4)
                        input_kw_array2.append(input_kw_1)
                        input_kw_array2.append(input_kw_1)
                        input_kw_array2.append(input_kw_1)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)

                    elif input_mass_4 > input_mass_3 and input_mass_4 > input_mass_1:
                        input_mass_array.append(input_mass_3)
                        input_mass_array.append(input_mass_1)
                        input_mass_array2.append(input_mass_4)
                        input_mass_array2.append(input_mass_4)
                        input_mass_array2.append(input_mass_4)
                        OUTPUT_ID = input_id_4
                        input_ID_array.append(input_id_3)
                        input_ID_array.append(input_id_1)
                        input_ID_array2.append(input_id_4)
                        input_ID_array2.append(input_id_4)
                        input_ID_array2.append(input_id_4)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_3)
                        input_kw_array.append(input_kw_1)
                        input_kw_array2.append(input_kw_4)
                        input_kw_array2.append(input_kw_4)
                        input_kw_array2.append(input_kw_4)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    else:
                        input_mass_array.append(input_mass_1)
                        input_mass_array.append(input_mass_4)
                        input_mass_array2.append(input_mass_3)
                        input_mass_array2.append(input_mass_3)
                        input_mass_array2.append(input_mass_3)
                        OUTPUT_ID = input_id_3

                        input_ID_array.append(input_id_1)
                        input_ID_array.append(input_id_4)
                        input_ID_array2.append(input_id_3)
                        input_ID_array2.append(input_id_3)
                        input_ID_array2.append(input_id_3)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_1)
                        input_kw_array.append(input_kw_4)
                        input_kw_array2.append(input_kw_3)
                        input_kw_array2.append(input_kw_3)
                        input_kw_array2.append(input_kw_3)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                if input_mass_3 < input_mass_2 and input_mass_3 < input_mass_1 and input_mass_3 < input_mass_4:
                    input_mass_array.append(input_mass_3)
                    input_ID_array.append(input_id_3)
                    input_kw_array.append(input_kw_3)

                    if input_mass_1 > input_mass_2 and input_mass_1 > input_mass_4:
                        input_mass_array.append(input_mass_2)
                        input_mass_array.append(input_mass_4)
                        input_mass_array2.append(input_mass_1)
                        input_mass_array2.append(input_mass_1)
                        input_mass_array2.append(input_mass_1)
                        OUTPUT_ID = input_id_1
                        input_ID_array.append(input_id_2)
                        input_ID_array.append(input_id_4)
                        input_ID_array2.append(input_id_1)
                        input_ID_array2.append(input_id_1)
                        input_ID_array2.append(input_id_1)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_2)
                        input_kw_array.append(input_kw_4)
                        input_kw_array2.append(input_kw_1)
                        input_kw_array2.append(input_kw_1)
                        input_kw_array2.append(input_kw_1)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    elif input_mass_4 > input_mass_2 and input_mass_4 > input_mass_1:
                        input_mass_array.append(input_mass_2)
                        input_mass_array.append(input_mass_1)
                        input_mass_array2.append(input_mass_4)
                        input_mass_array2.append(input_mass_4)
                        input_mass_array2.append(input_mass_4)
                        OUTPUT_ID = input_id_4
                        input_ID_array.append(input_id_2)
                        input_ID_array.append(input_id_1)
                        input_ID_array2.append(input_id_4)
                        input_ID_array2.append(input_id_4)
                        input_ID_array2.append(input_id_4)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_2)
                        input_kw_array.append(input_kw_1)
                        input_kw_array2.append(input_kw_4)
                        input_kw_array2.append(input_kw_4)
                        input_kw_array2.append(input_kw_4)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    else:
                        input_mass_array.append(input_mass_1)
                        input_mass_array.append(input_mass_4)
                        input_mass_array2.append(input_mass_2)
                        input_mass_array2.append(input_mass_2)
                        input_mass_array2.append(input_mass_2)
                        OUTPUT_ID = input_id_2
                        input_ID_array.append(input_id_1)
                        input_ID_array.append(input_id_4)
                        input_ID_array2.append(input_id_2)
                        input_ID_array2.append(input_id_2)
                        input_ID_array2.append(input_id_2)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_1)
                        input_kw_array.append(input_kw_4)
                        input_kw_array2.append(input_kw_2)
                        input_kw_array2.append(input_kw_2)
                        input_kw_array2.append(input_kw_2)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)

                if input_mass_4 < input_mass_2 and input_mass_4 < input_mass_1 and input_mass_4 < input_mass_3:
                    input_mass_array.append(input_mass_4)
                    input_ID_array.append(input_id_4)
                    input_kw_array.append(input_kw_4)

                    if input_mass_1 > input_mass_2 and input_mass_1 > input_mass_3:
                        input_mass_array.append(input_mass_2)
                        input_mass_array.append(input_mass_3)
                        input_mass_array2.append(input_mass_1)
                        input_mass_array2.append(input_mass_1)
                        input_mass_array2.append(input_mass_1)
                        OUTPUT_ID = input_id_1
                        input_ID_array.append(input_id_2)
                        input_ID_array.append(input_id_3)
                        input_ID_array2.append(input_id_1)
                        input_ID_array2.append(input_id_1)
                        input_ID_array2.append(input_id_1)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_2)
                        input_kw_array.append(input_kw_3)
                        input_kw_array2.append(input_kw_1)
                        input_kw_array2.append(input_kw_1)
                        input_kw_array2.append(input_kw_1)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)

                    elif input_mass_3 > input_mass_2 and input_mass_3 > input_mass_1:
                        input_mass_array.append(input_mass_2)
                        input_mass_array.append(input_mass_1)
                        input_mass_array2.append(input_mass_3)
                        input_mass_array2.append(input_mass_3)
                        input_mass_array2.append(input_mass_3)
                        OUTPUT_ID = input_id_3
                        input_ID_array.append(input_id_2)
                        input_ID_array.append(input_id_1)
                        input_ID_array2.append(input_id_3)
                        input_ID_array2.append(input_id_3)
                        input_ID_array2.append(input_id_3)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_2)
                        input_kw_array.append(input_kw_1)
                        input_kw_array2.append(input_kw_3)
                        input_kw_array2.append(input_kw_3)
                        input_kw_array2.append(input_kw_3)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                    else:
                        input_mass_array.append(input_mass_1)
                        input_mass_array.append(input_mass_3)
                        input_mass_array2.append(input_mass_2)
                        input_mass_array2.append(input_mass_2)
                        input_mass_array2.append(input_mass_2)
                        OUTPUT_ID = input_id_2
                        input_ID_array.append(input_id_1)
                        input_ID_array.append(input_id_3)
                        input_ID_array2.append(input_id_2)
                        input_ID_array2.append(input_id_2)
                        input_ID_array2.append(input_id_2)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)
                        output_ID_array.append(output_id)

                        input_kw_array.append(input_kw_1)
                        input_kw_array.append(input_kw_3)
                        input_kw_array2.append(input_kw_2)
                        input_kw_array2.append(input_kw_2)
                        input_kw_array2.append(input_kw_2)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)
                        output_kw_array.append(output_kw)

            count += 1
            continue

    final_bh_mass = np.zeros(len(input_mass_array))
    if len(final_bh_mass) > 0:
        final_bh_mass[0] = mass_bh
    else:
        final_bh_mass = []
    
    
    
    df2 = pd.DataFrame({"Type": collision_type, "BH ID": bh_ID, "Time" : t_array, "Mass1" : input_mass_array, "Mass2" : input_mass_array2, "MergerMass" : m_total_array, "Mass1_IDs" : input_ID_array, "Mass2_IDs" : input_ID_array2, "Merger_IDs" : output_ID_array, "Mass1_type" :input_kw_array, "Mass2_type" : input_kw_array2 , "Merger_type" : output_kw_array, "Final_bh_mass" : final_bh_mass})

    df = df.append(df2,ignore_index = True)[df2.columns.tolist()]
    ids_to_check  = (input_ID_array2)
    
    return OUTPUT_ID, df, ids_to_check

print("Finished defining functions...")

#Model Parameters:
rvs = ['1']#str
fbhs = ['0.05']#str
iterations = [0]#int
path = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/'
DataDir = './data/'
# path = '/projects/b1095/egp8636/elena_cmc/rundir/triples/'

print("Finished model parameters...")

m = 0 #iteration index count
for rv in rvs:
    print(f'rv = {rv}')
    for hmbf in fbhs:
        print(f'hmbf = {hmbf}')
        for e in range(0,iterations[m] + 1):
            print(f'e = {e}')
            m += 1

            #change this based on your files
            model = f'{path}rv{rv}rg8z0.0017n8e5/w5fb0.05fbh{hmbf}alpha3-2.0/{e}_{e}/'
            print(f"model = {path}rv{rv}rg8z0.0017n8e5/w5fb0.05fbh{hmbf}alpha3-2.0/{e}_{e}/")

            #units 
            units = read_units(model+'initial')
            time_units = units[0]['t_myr']
            mass_units = units[0]['mstar_msun']


            # Determining the BHs we want to look at:
            bh_idss = np.genfromtxt(model + 'initial.bhformation.dat', usecols = (3))
            bh_time = np.genfromtxt(model+ 'initial.bhformation.dat', usecols = (0))
            bh_masses = np.genfromtxt(model + 'initial.bhformation.dat', usecols = (6))

            # change here the mass at which you want to filter 
            bh_mass = [81.8923]
            bh_ids  = [864151]
            bh_time = [0.0678917901117370715]

            # Ignore these, they were failed attempts
            # bh_mass = [42.9484]
            # bh_ids =  [864150]
            # bh_time = [0.00241723323394524032]
            # masses_in_range = [bh_masses[i] for i in range(len(bh_idss)) if (bh_masses[i] > 50 and bh_masses[i] < 90)]
            # print(f'masses_in_range = {masses_in_range}')
            # bh_time = [bh_time[i] for i in range(len(bh_idss)) if bh_masses[i] == bh_mass[0]]
            # bh_ids = [bh_idss[i] for i in range(len(bh_idss)) if bh_masses[i] == bh_mass[0]]
            
            # These were the original three lines
            # bh_time = [bh_time[i] for i in range(len(bh_idss)) if bh_masses[i] >= 100.0]
            # bh_ids = [bh_idss[i] for i in range(len(bh_idss)) if bh_masses[i] >= 100.0]
            # bh_mass = [bh_masses[i] for i in range(len(bh_masses)) if bh_masses[i] >= 100.0]

            print(f'bh_time = {bh_time}')
            print(f'bh_ids = {bh_ids}')
            print(f'bh_mass = {bh_mass}')

            # THIS RUNS THE FUNCTIONS TO GET BH FORMATION HISTORY

            for z in range(len(bh_time)):

                #data frame set up
                columns= ["BH ID", "Type" ,"Time", "Mass1" , "Mass2" , "MergerMass" , "Mass1_IDs" , "Mass2_IDs" , "Merger_IDs" , "Mass1_type" , "Mass2_type" , "Merger_type" , "Final_bh_mass"]
                df = pd.DataFrame(columns=columns)

                ID = bh_ids[z]
                stop_time = bh_time[z]
                mass_bh = bh_mass[z]

                print('IMBH ID = ', ID, ' with mass ', mass_bh, ' at time ', stop_time*time_units)

                #Access the collision file and find the line number at which you exceed BH formation time
                f = open(model+'initial.collision.log','r')
                lines = f.readlines()
                stop_num = 0
                for j in range(1,len(lines)):
                    line = lines[j]
                    line = line.split('\n')
                    line = line[0]
                    data = line.split(' ')
                    t = data[0]
                    t = t.split('=')
                    t = float(t[1])
                    if t >= stop_time:
                        if j < len(lines)-1:
                            stop_num = j+1#+1 is just in case the next line has the same time
                        else:
                            stom_num = j
                        break

                if stop_num == 0:
                    stop_num = j

                # Now find collision history of star with given ID
                OUTPUT_IDs = [ID] #start by putting the BH ID
                ids_tested = [] #will be used to track repeated IDs 
                ids_tested2 = []

                for OUTPUT_ID in (OUTPUT_IDs):

                    semerge_ids = []
                    semerge_ids.append(OUTPUT_ID)

                    if OUTPUT_ID not in ids_tested:

                        print('Looking at collisions...')

                        OUTPUT_ID2 , df2, ids_to_check  = collision(OUTPUT_ID, stop_num,lines, mass_bh, ID)
                        df = df.append(df2,ignore_index = True)[df2.columns.tolist()]
                        semerge_ids.extend(ids_to_check)
                        if OUTPUT_ID2 != OUTPUT_ID: #if a new ID is outputted
                            OUTPUT_IDs.append(int(OUTPUT_ID2))
                        ids_tested.append(OUTPUT_ID)

                    for id in semerge_ids:
                        if id not in ids_tested2:
                            print('Looking in semerge file...')
                            OUTPUT_ID3, df2 = check_semerge(model, id,stop_time, ID)
                            df = df.append(df2,ignore_index = True)[df2.columns.tolist()]
                            ids_tested2.append(id)
                            
                            if OUTPUT_ID3 != OUTPUT_ID: #if a new ID is outputted
                                OUTPUT_IDs.append(int(OUTPUT_ID3))


                rslt_df = df.sort_values(by = 'Time')
                rslt_df[columns].to_csv(f'{DataDir}IMBH_histories/id={ID}_mass={bh_mass[0]}.csv', index=False)











