function[Record_Best_Chromasomes]=GA_Selection_of_Terms_for_Test_Simple_DOE11(X_Full_Train,Y_Train_Start,NofVars,NofTrials,...
TP,ML,Pop,Gen_1,Fitness_Flag,No_of_Iterations)

%XXXXXXXXXXXXXXXXXXXXXXX
%This version simply runs the GA for 15 trials (driven by R^2 or Q^2) and returns the
%Chromasome that represents the best model for each trial. 
clearvars -except X_Full_Train Y_Train_Start NofVars NofTrials...
 TP ML Pop Gen_1 Fitness_Flag No_of_Iterations


%First load the data and the full set of Z-scales (1224 variables)
%clear all

%This function calls two other functions:-
    %f(x)=MLR_fitting_and_Q2_by_L10out_Variable_Selection_Test
    %f(x)=Deletion_of_Excess_Variables_in_a_GA

    
Z_Scales=X_Full_Train;
[Nrow, Ncol]=size(Z_Scales);
Record_Best_Chromasomes=zeros(NofTrials,Ncol+1);


No_Children=Pop;
WP=1;
Transfer_Pop=TP;

for Trial=1:NofTrials 

%completeness = (Trial/15);
    
%h=waitbar(completeness);

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%1)SET UP INITIAL POPULATION
%A) Set up the Initial Population (random set of "Pop" Chromasomes)
%Set up a random initial population by randomly switching NofVar genes in 
%each chromasome "on". Fitness assessed by R^2 of the model represented by
%the chromasomes.
IP_Chromasome = zeros(Pop,Ncol);
for IP=1:Pop
    Gene_Numbers=randperm(Ncol,NofVars);
    for Gene_No=1:NofVars
        for gene=1:Ncol
            if gene==Gene_Numbers(Gene_No)
            IP_Chromasome(IP,gene)=1;
            else           
            end
        end
    end
end


%Build the Z_scales matrix represented by these chromasomes and determine
%the R^2 value for a MLR model based upon these variables.

Fitness_IP=zeros(Pop,1);
for IP_1=1:Pop
    Counter_2=1;
    Trial_Z_Scales=zeros(Nrow,NofVars);
    for Gene_1=1:Ncol
        if IP_Chromasome(IP_1,Gene_1)==1
            Trial_Z_Scales(:,Counter_2)=Z_Scales(:,Gene_1);
            Counter_2=Counter_2+1;
        else
        end
    end
    Units=ones(Nrow,1);
    
    %Fit a MLR model (using regress)on the Z_scales and return the R2 value
    %
    Regress_Z_Scales=[Units Trial_Z_Scales];
        
    [~,~,~,~,stats]=regress(Y_Train_Start,Regress_Z_Scales);
    
    
    %[Q_Squared]=MLR_fitting_and_Q2_by_L10out_Variable_Selection_Test(Trial_Z_Scales,Y_Train_Start,No_of_Iterations);
    R_Squared=stats(1,1);
    
    
    Fitness_IP(IP_1,1)=R_Squared;
end

%Asess the Initial Population and select the best individual to progress to
%the first generation.
Initial_Population=[IP_Chromasome, Fitness_IP];
Sorted_IP=sortrows(Initial_Population,-(Ncol+1));

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%2) START THE GA
Flag=zeros(Gen_1,1);
%Record_Q_Squared=zeros(Nrow_4,Gen_1);
%Record_R_Squared=zeros(Nrow_4,Gen_1);

for i=1:Gen_1
    clearvars -except i Trial WP NAAR NofVars NofTrials No_of_Iterations TP ML...
Pop Gen_1 Fitness_Flag Flag Sorted_IP Units Y_Train_Start Record_Best_Chromasomes...
Record_Best_Fitness Record_Flag X_Full_Train Nrow Ncol No_Children Z_Scales...
Transfer_Pop Sorted_Population_3 Record_Full_Population_1 Record_R_Squared Record_Q_Squared...
Run_1 Record_Consensus_Chromasomes Best_Chromasomes Best_Consensus_Chromasome Best_Q_Squared
    
    
    
    if i==1
        Fit_Population=Sorted_IP(1:WP,:);
    else
        Fit_Population=Sorted_Population_3;
    end
    Fittest_Individual = Fit_Population(1:1,:);
    
    Parent_Population=Fit_Population(:,1:Ncol);
    No_Parents=size(Parent_Population,1);
    Children=zeros(No_Children,Ncol);
    %Mating via single point crossover
        if No_Parents==1    
            for kj=1:No_Children
                Children(kj,:)=Parent_Population(1,:);
            end
        elseif No_Parents==2
            for kk=1:No_Children/2;
                Parents = Parent_Population(1:2,:);
                Cop1= randi((Ncol-1),1);
                Children(2*kk-1:2*kk,1:Ncol) = [Parents(1,1:Cop1), Parents(2,Cop1+1:Ncol); Parents(2,1:Cop1), Parents(1,Cop1+1:Ncol)];
            end
        elseif No_Parents>2
            if No_Parents==3
                Tournament_K=2;
            else
                Tournament_K=6;
            end
            %TOURNAMENT SELECTION-AND FORMATION OF CHILDREN
            %Use tournament selection (K=2 or K=6) to select the best
            %parents and then cross to form new children
                
            
            Nrow_1=No_Parents;
            for kk=1:No_Children/2;
                for parent=1:2
                    Random_assignment = rand([Tournament_K,1]);
                    Test_set=zeros(Tournament_K,Ncol);
                    for kl=1:Tournament_K
                         Row = ceil(Random_assignment(kl)*Nrow_1);
                         Test_set(kl,:) = Parent_Population(Row,:);
                    end
                    Sorted_TS = sortrows(Test_set, -(NAAR+1));
                    Parents(parent,:) = Sorted_TS(1,:);
                end
                Cop1= randi((Ncol-1),1);
                Children(2*kk-1:2*kk,1:Ncol) = [Parents(1,1:Cop1), Parents(2,Cop1+1:Ncol); Parents(2,1:Cop1), Parents(1,Cop1+1:Ncol)];
            end
                
        end
     
    
    
    Nrow_3=size(Children,1);        
    
    %Up to six mutations (determined by ML?) are randomly inserted into
    %each of the children.
    
    Input_Sequences=Children;
    for ab=1:ML
        for bc=1:Nrow_3
            Mutation_Individual=Input_Sequences(bc,:);
            mutation_test = rand(1,1);
            
            Mutated_Gene = ceil(mutation_test*Ncol);
            Gene_Value = Mutation_Individual(1,Mutated_Gene);
            I_to_be_Modified=Mutation_Individual;
            if Gene_Value==0
                I_to_be_Modified(1,Mutated_Gene)=1;
            elseif Gene_Value==1
                I_to_be_Modified(1,Mutated_Gene)=0;
            end
            Modified_Individual=I_to_be_Modified;
            Input_Sequences(bc,:)=Modified_Individual;
        end
    end
    New_Individuals=Input_Sequences;
         
    
    Individual_to_Beat=Fittest_Individual(1,1:Ncol);    
    Full_Population=[New_Individuals;Individual_to_Beat]; 
    
    %Having got "Full_Population" need to check that the maximum
    %number of varaibles is not exceeded and if it is the number of
    %variables is reduced.
    Population_to_be_Measured=Deletion_of_Excess_Variables_in_a_GA(Full_Population,NofVars);
    Nrow_4=size(Population_to_be_Measured,1);
    
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %Convert chroamsomes into Z_Scales and obtain a fitness value for each
    %member of the population.
    Cell_Array_Trial_Z_Scales=cell(Nrow_4);
    No_Columns=sum(Population_to_be_Measured,2);
    Fitness_WP_R=zeros(Nrow_4,1);
    for Working_Pop=1:Nrow_4
        tic;
        Counter_3=1;
        Trial_Z_Scales=zeros(Nrow,No_Columns(Working_Pop,1));
        for Gene_2=1:Ncol
            if Population_to_be_Measured(Working_Pop,Gene_2)==1
            Trial_Z_Scales(:,Counter_3)=Z_Scales(:,Gene_2);
            Counter_3=Counter_3+1;
            else
            end
        end
        oldtimer=toc;
        tic;
        Trial_Z_Scales=Z_Scales(:,find(Population_to_be_Measured(Working_Pop,:)==1));
        newtimer=toc;
        sprintf('%s%g','Speedup=',double(oldtimer)/double(newtimer))
        Cell_Array_Trial_Z_Scales{Working_Pop}=Trial_Z_Scales;
        
        %Fit a MLR model (using regress)on the Z_scales and return the Q2 value
        %as calculated by Leave 10% out Cross Validation. This is the mean of
        %100 cross-validation runs.
        
        %[Q_Squared]=MLR_fitting_and_Q2_by_L10out_Variable_Selection_Test(Trial_Z_Scales,Y_Train_Start,No_of_Iterations);
        
        Regress_Z_Scales=[Units Trial_Z_Scales];
        
        [~,~,~,~,stats]=regress(Y_Train_Start,Regress_Z_Scales);
        
        Fitness_WP_R(Working_Pop,1)=stats(1,1);
        
    end
    
    Fitness_WP_Q=zeros(Nrow_4,1);
    if Fitness_Flag==+1
    %if Max_Fitness>0.0 %DOE requires that the GA is driven by Q^2 or a combination of Q^2 and R^2.
         
         for W_Pop=1:Nrow_4
             Trial_Z_Scales_1=Cell_Array_Trial_Z_Scales{W_Pop};
             [Q_Squared]=MLR_fitting_and_Q2_by_L10out_Variable_Selection_Test_vs2(Trial_Z_Scales_1,Y_Train_Start,No_of_Iterations);
             Fitness_WP_Q_LPO(W_Pop,1)=Q_Squared;
             %Record_Q_Squared(W_Pop,i)=Q_Squared;
         end
    elseif Fitness_Flag==0
         
         for W_Pop=1:Nrow_4
             Trial_Z_Scales_1=Cell_Array_Trial_Z_Scales{W_Pop};
            [Q_Squared]=MLR_fitting_and_Q2_by_LOO_from_X_Array_vs2(Trial_Z_Scales_1,Y_Train_Start);
            Fitness_WP_Q_LOO(W_Pop,1)=Q_Squared;
         end
    end
           
    if Fitness_Flag==-1       
        Measured_Population=[Population_to_be_Measured, Fitness_WP_R];
    elseif Fitness_Flag==0
        Measured_Population=[Population_to_be_Measured, Fitness_WP_Q_LOO];
    elseif Fitness_Flag==1
        Measured_Population=[Population_to_be_Measured, Fitness_WP_Q_LPO];
    end
    
    
   
    Record_Full_Population_1{i}=Measured_Population;
    Sorted_Population=sortrows(Measured_Population,-(Ncol+1));    
    Sorted_Population_3=Sorted_Population(1:Transfer_Pop,:); 
    Best_Q_Squared(i,1)=Sorted_Population_3(1,Ncol+1);
    
end
Record_Best_Chromasomes(Trial,:)=Sorted_Population_3(1,:);
end
