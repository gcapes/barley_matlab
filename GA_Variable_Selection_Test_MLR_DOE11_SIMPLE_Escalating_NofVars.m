%This script does the following:-
%A) generate a random set of 16AA sequences
%B) split this set of random sequences into a Test set and a Training set
%(roughly 50/50.

%Sections A and B have been previously done and the Test/Train splits are
%stored in 'Twenty_Training_and_Test_sets_for_GA_Variable_Sel_200'

%C) generate extended Z scales (includes squared terms and interactions)
%for both training and test set.
%D) use a GA to select variables to provide a good fit to the training
%data. Fitness values provided by FL-12.
%E) assess fits using bootstrapping and leave 10% out CV.
%F test final fits using the test set. collect statistically significant
%data on the impact of different fitting strategies.

%This program calls many other functions


%f(x)=Multiple_Models_Fitness_Value_16AA_vs2_with_error
%f(x)=GA_Selection_of_Terms_for_Test_Simple_DOE11- which in-turn calls:-
%f(x)=MLR_fitting_and_Q2_by_L10out_Variable_Selection_Test_vs2
%f(x)=MLR_fitting_and_Q2_by_LOO_from_X_Array_vs2
%f(x)=Deletion_of_Excess_Variables_in_a_GA

% clear all
rng('shuffle') %Randomize the seed for the random number generator.



load('Twenty_Training_and_Test_sets_for_GA_Variable_Sel_200.mat');

model=13;

for Run_V=1:1
    Design=...
        [1	0.05	100	20	-1	0
        2	0.05	100	20	0	0
        3	0.05	100	20	1	5
        4	0.05	100	20	1	20
        5	0.05	200	40	-1	0
        6	0.05	200	40	0	0
        7	0.05	200	40	1	5
        8	0.2     200	40	1	20
        9	0.2     100	20	-1	0
        10	0.2     100	20	0	0
        11	0.2     100	20	1	5
        12	0.2     100	20	1	20
        13	0.2     200	40	-1	0
        14	0.2     200	40	0	0
        15	0.2     200	40	1	5
        16	0.2     200	40	1	20];
    
    
    Selected_Design=Design(Run_V,:);
    Fitness_Flag=Selected_Design(1,5); %GA Driven by R^2 (-1) or Q^2-LOO(0) or Leave 10% out(+1)
    No_of_Iterations=Selected_Design(1,6); %Iterations for Q^2
    Gen_1=Selected_Design(1,3); %Number of Generations
    Pop=Selected_Design(1,4); %Mutation Level
    Error=Selected_Design(1,2);
    
    for Example=1:15
        clear X_Predict X_Train
        
        %completeness = ((Run_V-1)*15+Example)/240;
        completeness=Example/15;
        h=waitbar(completeness);
        
        NAAR=16;
        Number_Required=200;
        NofVars=18;
        NofTrials=15;
        %No_of_Iterations=20;
        TP=2;
        ML=1;
        %Pop=20;
        %Gen_1=60;
        Y_Flag=-1; %Y Response with error
        
        
        
        
        Test_Set_Start=Test_Set_Cell_Array{Example};
        Training_Set_Start=Training_Set_Cell_Array{Example};
        
        %clear Test_Set Training_Set
        %6) The Z Scales for the Training set (Start) are expanded to include squared terms
        %and all interactions.
        clear Squared_Terms Cross_Terms
        for m=1:48
            Squared_Terms(:,m) = Training_Set_Start(:,m).^2;
        end
        Counter_3=1;
        for o=1:48
            for p=1:o-1
                Cross_Terms(:,Counter_3) = Training_Set_Start(:,o).*Training_Set_Start(:,p);
                Counter_3=Counter_3+1;
            end
        end
        
        X_Full_Train=[Training_Set_Start, Squared_Terms, Cross_Terms];
        [Nrow_1,Ncol]=size(X_Full_Train);
        
        %7) The same thing can be done for the Test set
        clear Squared_Terms Cross_Terms
        for r=1:48
            Squared_Terms(:,r) = Test_Set_Start(:,r).^2;
        end
        
        Counter_4=1;
        for s=1:48
            for t=1:s-1
                Cross_Terms(:,Counter_4) = Test_Set_Start(:,s).*Test_Set_Start(:,t);
                Counter_4=Counter_4+1;
            end
        end
        clear X_Full_Test Fitness_2 Fitness_2_raw Y_Train
        X_Full_Test=[Test_Set_Start, Squared_Terms, Cross_Terms];
        [Nrow_2,Ncol]=size(X_Full_Test);
        %8) The original set of Z_scales for the training set is used to obtain
        %Fitness values for the fitting. The Z_scales are found in "Training_Set_Start".
        
        
        STS=size(Training_Set_Start,1);
        
        if Y_Flag==1
            Fitness_2 = Multiple_Models_Fitness_Value_16AA_vs2(model,STS,NAAR,Training_Set_Start);
        elseif Y_Flag==-1
            Fitness_2_raw = Multiple_Models_Fitness_Value_16AA_vs2_with_error(model,STS,NAAR,Training_Set_Start,Error);
            Fitness_2=Fitness_2_raw(:,2); %Values with error are in second column.
        end
        Y_Train_Start=Fitness_2; %Activity_for_Training_Set
        
        %9) Similarly for the Test Set. The original set of Z_scales is used to obtain
        %Fitness values for testing the correlation. The Z_scales are found in "Test_Set_Start".
        clear Fitness_3 Y_Test_Start Fitness_3_raw
        
        MR=size(Test_Set_Start,1);
        
        if Y_Flag==1
            Fitness_3 = Multiple_Models_Fitness_Value_16AA_vs2(model,MR,NAAR,Test_Set_Start);
        elseif Y_Flag==-1
            Fitness_3_raw = Multiple_Models_Fitness_Value_16AA_vs2_with_error(model,MR,NAAR,Test_Set_Start,Error);
            Fitness_3=Fitness_3_raw(:,2); %Values with error are in second column.
        end
        Y_Test_Start=Fitness_3; %Activity_for_Test_Set
        
        %10)X_Full_Train and Y_Train are used as inputs to the Variable Selection GA which returns a
        %set of chromasomes (selected variables) with accompanying R^2 or Q^2 values
        clear Sorted_GA_Chromasomes Best_GA_Chromasome_NofVars
        for i=1:6
            NofVars=i*3;
            [Record_GA_Chromasomes]=GA_Selection_of_Terms_for_Test_Simple_DOE11(X_Full_Train,Y_Train_Start,NofVars,NofTrials,...
                TP,ML,Pop,Gen_1,Fitness_Flag,No_of_Iterations);
            
            
            %Whether the GA has been driven by R^2 or Q^2 the first step is to
            %determine the best chromasome based upon a simple sorting of the values in
            %column 1225.
            Sorted_GA_Chromasomes=sortrows(Record_GA_Chromasomes,-(Ncol+1));
            Best_GA_Chromasome_NofVars(i,:)=Sorted_GA_Chromasomes(1,1:Ncol+1);
        end
        Record_Best_Chromasomes{Example}=Best_GA_Chromasome_NofVars;
        %Then select the best chromasome of the set with NofVars=3,6,9,12,15,and
        %18.
        Sorted_Chromasomes_Escalating_NofVars=sortrows(Best_GA_Chromasome_NofVars,-(Ncol+1));
        Best_GA_Chromasome=Sorted_Chromasomes_Escalating_NofVars(1,:);
        %11) The Best Chromasome is used to select the best variables from which
        %a correlation is derived. Likely to be less than 21.
        
        clear Units Trial_Z_Scales
        Counter_5=1;
        for Gene_1=1:Ncol
            if Best_GA_Chromasome(1,Gene_1)==1
                Trial_Z_Scales(:,Counter_5)=X_Full_Train(:,Gene_1);
                Counter_5=Counter_5+1;
            else
            end
        end
        
        Units=ones(Nrow_1,1); %Nrow_1=size of training set
        X_Train=[Units, Trial_Z_Scales];
        %Fit a MLR model (using regress)on the
        %Z_scales and return the R2 and Q2 value.
        %R Squared
        [b,bint,r,rint,stats]=regress(Y_Train_Start,X_Train);
        
        Coefficients=b';
        %12) Having got the coefficients that represent the best model use them to predict values for the test set.
        %First generate the Z_Scales array for the test set.
        clear Y_Test_Predict Trial_Z_Scales Units Predict_Array
        Counter_6=1;
        for Gene_1=1:Ncol
            if Best_GA_Chromasome(1,Gene_1)==1
                Trial_Z_Scales(:,Counter_6)=X_Full_Test(:,Gene_1);
                Counter_6=Counter_6+1;
            else
            end
        end
        %Then calculate predicted values for the testset.
        
        Units=ones(Nrow_2,1);
        X_Predict=[Units, Trial_Z_Scales];
        for ij=1:Nrow_2
            for ik=1:Counter_6
                Predict_Array(ij,ik)=X_Predict(ij,ik).*Coefficients(1,ik);
            end
        end
        Y_Test_Predict=sum(Predict_Array,2);
        
        %13) Finally calculate a R^2(Predict) value for the test set.
        clear Sum_of_Squares_Vector Residuals
        Mean_Test_Set=mean(Y_Test_Start);
        Sum_of_Squares_Vector=(Y_Test_Start-Mean_Test_Set).^2;
        SOS_Total=sum(Sum_of_Squares_Vector,1); %Total Sum of Squares for Test Set
        
        Residues=Y_Test_Predict-Y_Test_Start;
        SOS_Predict=sum(Residues.^2,1);
        
        R_Squared_Testset=1-(SOS_Predict/SOS_Total);
        
        
        Summary(Example,:)=[Best_GA_Chromasome(1,1:Ncol+1),R_Squared_Testset];
        
    end
    Record_Summary_Cell_Array{Run_V}=Summary;
end
