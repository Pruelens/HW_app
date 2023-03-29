import random
import numpy as np
import streamlit as st
import pandas as pd
import plotly.express as px

def simulate_HW_twolocus(pop_size, generations, start_p_A, link_AB, W_genotypes_Aa, W_genotypes_Bb,mig_rate, mu_rate,rec_rate):
    # Set initial parameters
    pop_size = pop_size              # population size
    genotypes = ['AA','Aa','aa']     # possible genotypes in the population
    genotypes_2 = ['BB','Bb','bb']  
    rec_rate = rec_rate #recombination rate

    genotypes_linked = [['AA', 'BB'], ['Aa', 'Bb'], ['aa', 'bb']]

    W_AA, W_Aa, W_aa = W_genotypes_Aa  # fitness values for each genotype
    W_BB, W_Bb, W_bb = W_genotypes_Bb
    p_A = start_p_A                      # initial frequency of allele A
    p_B = start_p_A #start_p_B                      # initial frequency of allele B
    _weight_A = [p_A*p_A, 2*p_A*(1-p_A), (1-p_A)**2]  # calculate genotype frequencies based on allele frequency
    _weight_B = [p_B*p_B, 2*p_B*(1-p_B), (1-p_B)**2]  # calculate genotype frequencies based on allele frequency
    


    nr_replicates = 5

    migration_pop_structure = ([1, 1, 1],[1, 1, 1]) #Probability of 'AA','Aa','aa' in the source population
    Population_A_out = random.choices(genotypes,weights = migration_pop_structure[0], k = pop_size)  
    Population_B_out = random.choices(genotypes_2,weights = migration_pop_structure[1], k = pop_size)
    Population_outside = [list(a) for a in zip(Population_A_out, Population_B_out)]

    progress_text = "Simulation in progress. Please wait."
    my_bar = st.progress(0, text=progress_text)
    percent_complete = 1
    
    
    # Iterate over generations
    generations = generations
    output_population  = pd.DataFrame(columns = ['gen','Replicate','freq_AA', 'freq_aa', 'freq_Aa', 'freq_A', 
                                                 'freq_a', 'freq_BB', 'freq_bb', 'freq_Bb', 'freq_B', 'freq_b'])# store the population frequencies after each generation
    
    K = pop_size                     # carying capacity'

    for repl in range(nr_replicates):
        if link_AB == True:
            Population = random.choices(genotypes_linked,weights = _weight_A, k = pop_size)
        else:
            Population_A = random.choices(genotypes,weights = _weight_A, k = pop_size)  
            Population_B = random.choices(genotypes_2,weights = _weight_B, k = pop_size)
            Population = [list(a) for a in zip(Population_A, Population_B)]
        
        if (len(Population)%2) == 0:
            pass
        else:
            del Population[-1]

        for gen in range(generations):
            progress_per = percent_complete/(nr_replicates * generations)
            if progress_per == 1:
                my_bar.progress(progress_per, text='Simulation Finished')
            else:
                my_bar.progress(progress_per, text=progress_text)
            percent_complete += 1

            Pop_len = len(Population)    # get the current population size
            
            # Calculate genotype and allele frequencies in the population
            A_pop,B_pop = np.array(Population).T
            A_pop,B_pop = list(A_pop),list(B_pop)

            freq_AA = A_pop.count('AA') / Pop_len
            freq_aa = A_pop.count('aa')/ Pop_len
            freq_Aa = (A_pop.count('Aa') + A_pop.count('aA')) / Pop_len
            freq_A = ''.join(A_pop).count('A') / (Pop_len*2)
            freq_a = ''.join(A_pop).count('a') / (Pop_len*2)
            
            freq_BB = B_pop.count('BB') / Pop_len
            freq_bb = B_pop.count('bb')/ Pop_len
            freq_Bb = (B_pop.count('Bb') + B_pop.count('bB')) / Pop_len
            freq_B = ''.join(B_pop).count('B') / (Pop_len*2)
            freq_b = ''.join(B_pop).count('b') / (Pop_len*2)
            # Add the frequencies to the output list
            output_population.loc[len(output_population)] = list((gen,repl,freq_AA, freq_aa, freq_Aa, freq_A, 
                                                                freq_a, freq_BB, freq_bb, freq_Bb, freq_B, freq_b))
            
            # Generate new population by randomly mating individuals in the current population
            random.shuffle(Population)   # shuffle the population for random mating

            new_generation = []
            for i in range(0, Pop_len, 2):
                p1 = Population[i] #e.g. ['Aa', 'Bb']
                p2 = Population[i+1]

                #parent 1 p1
                p1_gametes = []
                for _i in range(4): #number of gametes
                    if random.random() < rec_rate:
                        _r = random.sample([[0,1],[1,0]],k=1)[0]
                        p1_gametes.append((p1[0][_r[0]], p1[1][_r[1]]))
                    else:
                        if random.random() < 0.5:
                            p1_gametes.append((p1[0][0], p1[1][0]))
                        else: 
                            p1_gametes.append((p1[0][1], p1[1][1]))

                

                #parent 2 p2
                p2_gametes = []
                for _i in range(4): #number of gametes
                    if random.random() < rec_rate:
                        _r = random.sample([[0,1],[1,0]],k=1)[0]
                        p2_gametes.append((p2[0][_r[0]], p2[1][_r[1]]))
                    else:
                        if random.random() < 0.5:
                            p2_gametes.append((p2[0][0], p2[1][0]))
                        else: 
                            p2_gametes.append((p2[0][1], p2[1][1]))

                for offspring in range(4):
                    new_individual = [p1_gametes[offspring][0]+ p2_gametes[offspring][0], p1_gametes[offspring][1]+ p2_gametes[offspring][1]]
                    new_generation.append(new_individual)


            Population = new_generation


            # Create a new population by selecting individuals based on their fitness values
            new_pop = []

            f_AABB = sum((W_AA,W_BB))/2* 1/((Pop_len/K)*2)
            f_AABb = sum((W_AA,W_Bb))/2 * 1/((Pop_len/K)*2)
            f_AAbb = sum((W_AA,W_bb))/2 * 1/((Pop_len/K)*2)

            f_AaBB = sum((W_Aa,W_BB))/2 * 1/((Pop_len/K)*2)
            f_AaBb = sum((W_Aa,W_Bb))/2 * 1/((Pop_len/K)*2)
            f_Aabb = sum((W_Aa,W_bb))/2 * 1/((Pop_len/K)*2)

            f_aaBB = sum((W_aa,W_BB))/2 * 1/((Pop_len/K)*2)
            f_aaBb = sum((W_aa,W_Bb))/2 * 1/((Pop_len/K)*2)
            f_aabb = sum((W_aa,W_bb))/2 * 1/((Pop_len/K)*2)
        


            for indiv in Population:
                # Introduce mutations
                mut_indiv = indiv
                # for i in range(2):
                #     if random.random() < mu_rate:
                #         if mut_indiv[i] == 'A':
                #             mut_indiv = mut_indiv[:i] + 'a' + mut_indiv[i+1:]
                #         else:
                #             mut_indiv = mut_indiv[:i] + 'A' + mut_indiv[i+1:]
                            
                r = random.random()
                if mut_indiv[0]+mut_indiv[1] == 'AABB':
                    if r < f_AABB:
                        new_pop.append(mut_indiv)
                elif mut_indiv[0]+mut_indiv[1] in ('AABb','AAbB') :
                    if r < f_AABb:
                        new_pop.append(mut_indiv)
                elif mut_indiv[0]+mut_indiv[1] == 'AAbb':
                    if r < f_AAbb:
                        new_pop.append(mut_indiv)  
                elif mut_indiv[0]+mut_indiv[1] in ('AaBB','aABB'):
                    if r < f_AaBB:
                        new_pop.append(mut_indiv)
                elif mut_indiv[0]+mut_indiv[1]  in ('AaBb','AabB','aABb','aAbB'):
                    if r < f_AaBb:
                        new_pop.append(mut_indiv) 
                elif mut_indiv[0]+mut_indiv[1] in ('Aabb','aAbb'):
                    if r < f_Aabb:
                        new_pop.append(mut_indiv)
                elif mut_indiv[0]+mut_indiv[1] == 'aaBB':
                    if r < f_aaBB:
                        new_pop.append(mut_indiv)    
                elif mut_indiv[0]+mut_indiv[1] in ('aaBb','aabB'):
                    if r < f_aaBb:
                        new_pop.append(mut_indiv)
                elif mut_indiv[0]+mut_indiv[1] == 'aabb':
                    if r < f_aabb:
                        new_pop.append(mut_indiv)    
                                
                else:
                    print(ut_indiv[0]+mut_indiv[1])

            #Introduce migrants
            Pop_len = len(new_pop)
            num_migrants = np.random.binomial(Pop_len, mig_rate)
            for i in range(num_migrants):
            #     # Create a migrant individual
                 migrant = random.sample(Population_outside,1)[0]
                 new_pop.append(migrant)

            if (len(new_pop)%2) == 0:
                pass
            else:
                del new_pop[-1]
            Population = new_pop
    

    
    return output_population

# Define the Streamlit app
def app():
    # Define the app header
    st.header("Explore Hardy-weinberg Equilibrium")

    st.markdown('This is a simulation app that models the evolution of two loci. '
            'It uses a Wright-Fisher model, where a population is randomly mated for a given number of generations. '
            'The app allows you to set various parameters such as population size, the starting allele frequencies of the A locus, '
            'fitness values for each genotype, migration rate, mutation rate, and recombination rate. The app then generates '
            'a plot of the frequencies of the different genotypes over the number of generations specified. '
            'This app can be used to understand how genetic drift, mutation, and selection affect the evolution of populations over time.')


    # Define the input fields
    start_p_A = st.slider("Allele frequency A (p)", min_value=0.0, max_value=1.0, step=0.01, value=0.5)

    N = st.slider("Population size (N)", min_value=10, max_value=1000, step=10, value=100)
    
    generations = st.slider("Generations", min_value=10, max_value=1000, step=10, value=100)

    st.markdown('Set to fitness values of each locus separately. Final fitness is the average fitness of A and B genotype.')
    link_AB = st.checkbox('Link loci A and B at the start?')
    col1, col2, col3 = st.columns(3)



    st.markdown('#')
    
    with col1:
        W_AA = st.number_input('Fitness AA',min_value=0.00, max_value=1.00, value=1.00, step=0.05)
        W_BB = st.number_input('Fitness BB',min_value=0.00, max_value=1.00, value=1.00, step=0.05)

    with col2:
        W_Aa = st.number_input('Fitness Aa',min_value=0.00, max_value=1.00, value=1.00, step=0.05)
        W_Bb = st.number_input('Fitness Bb',min_value=0.00, max_value=1.00, value=1.00, step=0.05)

    with col3:
        W_aa = st.number_input('Fitness aa',min_value=0.00, max_value=1.00, value=1.00, step=0.05)
        W_bb = st.number_input('Fitness bb',min_value=0.00, max_value=1.00, value=1.00, step=0.05)

    col1_2, col2_2,col2_3 = st.columns(3)
    with col1_2:
        mig_rate = st.number_input('Migration rate',min_value=0.00, max_value=1.00, value=0.00, step=0.05)
    with col2_2:
        mu_rate = st.number_input('Mutation rate',min_value=0.00, max_value=1.00, value=0.00, step=0.05)  
    with col2_3:
        rec_rate = st.number_input('Recombination rate',min_value=0.00, max_value=1.00, value=0.00, step=0.05)  
    W_genotypes_Aa = [W_AA, W_Aa, W_aa]
    W_genotypes_Bb = [W_BB, W_Bb, W_bb]

    # Run the simulation

    if st.button('Run Simulation'):
        data = simulate_HW_twolocus(N, generations, start_p_A, link_AB, W_genotypes_Aa, W_genotypes_Bb,mig_rate, mu_rate,rec_rate)

        # Show the results
        st.subheader("Simulation Results")
        #st.write(f"Genotype frequencies: ")

        fig_allele_A = px.line(data, x="gen", y="freq_A", color='Replicate', 
                    range_y=[0,1], range_x=[0,generations])
        fig_allele_A.update_layout(xaxis_title="Generations",yaxis_title='Allele Frequency A', showlegend=False)

        fig_allele_B = px.line(data, x="gen", y="freq_B", color='Replicate', 
                    range_y=[0,1], range_x=[0,generations])
        fig_allele_B.update_layout(xaxis_title="Generations",yaxis_title='Allele Frequency B', showlegend=False)


        tab1, tab2 = st.tabs(["Allele Frequencies", "Genotype Frequencies"])
        with tab1:

            tab1_1, tab1_2 = st.tabs(["Allele A", "Allele B"])
            with tab1_1:
                st.plotly_chart(fig_allele_A, use_container_width=False, sharing="streamlit", theme="streamlit")
            with tab1_2:
                st.plotly_chart(fig_allele_B, use_container_width=False, sharing="streamlit", theme="streamlit")

        with tab2:
            tab2_1, tab2_2, tab2_3,tab2_4, tab2_5, tab2_6 = st.tabs(["Frequency AA", "Frequency Aa","Frequency aa",
                                              "Frequency BB", "Frequency Bb","Frequency bb"])
            with tab2_1:
                fig_AA = px.line(data, x="gen", y="freq_AA", color='Replicate', 
                        range_y=[0,1], range_x=[0,generations])
                fig_AA.update_layout(xaxis_title="Generations",yaxis_title='Frequency AA', showlegend=False)
                st.plotly_chart(fig_AA, use_container_width=False, sharing="streamlit", theme="streamlit")
            with tab2_2:
                fig_Aa = px.line(data, x="gen", y="freq_Aa", color='Replicate', 
                        range_y=[0,1], range_x=[0,generations])
                fig_Aa.update_layout(xaxis_title="Generations",yaxis_title='Frequency Aa', showlegend=False)
                st.plotly_chart(fig_Aa, use_container_width=False, sharing="streamlit", theme="streamlit")
            with tab2_3:
                fig_aa = px.line(data, x="gen", y="freq_aa", color='Replicate',
                        range_y=[0,1], range_x=[0,generations])
                fig_aa.update_layout(xaxis_title="Generations",yaxis_title='Frequency aa', showlegend=False)
                st.plotly_chart(fig_aa, use_container_width=False, sharing="streamlit", theme="streamlit")

            with tab2_4:
                fig_BB = px.line(data, x="gen", y="freq_BB", color='Replicate', 
                        range_y=[0,1], range_x=[0,generations])
                fig_BB.update_layout(xaxis_title="Generations",yaxis_title='Frequency BB', showlegend=False)
                st.plotly_chart(fig_BB, use_container_width=False, sharing="streamlit", theme="streamlit")
            with tab2_5:
                fig_Bb = px.line(data, x="gen", y="freq_Bb", color='Replicate', 
                        range_y=[0,1], range_x=[0,generations])
                fig_Bb.update_layout(xaxis_title="Generations",yaxis_title='Frequency Bb', showlegend=False)
                st.plotly_chart(fig_Bb, use_container_width=False, sharing="streamlit", theme="streamlit")
            with tab2_6:
                fig_bb = px.line(data, x="gen", y="freq_bb", color='Replicate',
                        range_y=[0,1], range_x=[0,generations])
                fig_bb.update_layout(xaxis_title="Generations",yaxis_title='Frequency bb', showlegend=False)
                st.plotly_chart(fig_bb, use_container_width=False, sharing="streamlit", theme="streamlit")


app()