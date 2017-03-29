function[Corrected_Chromasome]=Deletion_of_Excess_Variables_in_a_GA(Chromasome,NofVar)
%This script is to be used with "Reetz_Data_PLS_Model_by_GA_selection_of_Terms" 
% to put an upper limit on the number of variables selected. It takes the 
%chromasomes for a population and counts how many "1s" are present and if
%more than a preset number (NofVar) randomly converts 1s to 0s until the
%required number of variables is reached.
clearvars -except Chromasome NofVar
%Chromasome=Full_Population;


[Nrow,Ncol]=size(Chromasome);
Single_Chromasome=zeros(1,Ncol);
Total_No_Vars_coded=sum(Chromasome,2);
for ff=1:Nrow
    Single_Chromasome=Chromasome(ff,:);
    
    if Total_No_Vars_coded(ff,1)>NofVar
        Cycles=Total_No_Vars_coded(ff,1)-NofVar;
        for ee=1:Cycles 
            Counter_12=1;
            Total_No_Vars=sum(Single_Chromasome,2);
            Variable_Position=zeros(Total_No_Vars,1);
            %Determine where all the selected variables are.
            tic
            for dd=1:Ncol
                if Single_Chromasome(1,dd)==1 
                    Variable_Position(Counter_12,1)=dd;
                    Counter_12=Counter_12+1;
                else
                end
            end
            oldtimer=toc
            tic
            Variable_Position_vec=find(Single_Chromasome==1).';
            newtimer=toc;
            sprintf('%s%g','Speedup=',double(oldtimer)/double(newtimer))
            assert(isequal(Variable_Position_vec,Variable_Position))
            %Hence "Variable_Position" is a list of the selected genes (=1) in
            %Chromasome ff
        
            
            %Now randomly select one of the selected genes (=1) to deselect (change
            %to zero)
            mutation_test=[];
            mutation_test = rand(1,1);
            Nrow_10=size(Variable_Position,1);
            Deleted_Variable_Number = ceil(mutation_test*Nrow_10);
            Deleted_Variable_Position=Variable_Position(Deleted_Variable_Number,1);
            
            Single_Chromasome(1,Deleted_Variable_Position)=0;
            
        end
    else
    end
    Chromasome(ff,:)=Single_Chromasome;
end
%SUM=sum(Chromasome,2);
Corrected_Chromasome=Chromasome;

            
            
            
            
            
            
