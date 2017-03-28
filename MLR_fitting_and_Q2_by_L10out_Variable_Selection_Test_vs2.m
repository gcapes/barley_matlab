function[Q_Squared_Mean]=MLR_fitting_and_Q2_by_L10out_Variable_Selection_Test_vs2(X_Archive,Y_Archive,No_of_Iterations)
%This script takes an array of X values and a corresponding set of Y values
%and calculates a Q^2 type statistic for the MLR model mapping Y onto X
%using 10 fold (ie:- leave 10% out) Cross validation.

%20/10/14:- Corrected script so that Row_Numbers was defined outside the
%sample loop- otherwise summation of the testsets within one Q^2
%calculation doesn't cover the whole dataset
%of the Test set (not the full set).

%Using FL-12 (Model 13)
%X_Archive=Trial_Z_Scales_1;
%Y_Archive=Y_Train_Start;

%If appropriate autoscale X and Y variables;%%%%%%%%%%%
%X_Scaled=Autoscale_with_Reported_Mean_and_STD(X);
%The autoscale function returns the scaled X array with X_Mean and X_STD
%tacked onto the bottom.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%calculate the total sum-of-squares for the Y variable
%
%Calculate Q^2 by Leave 10% out Cross Validation
%Y_Archive=Y_Train_Start;
%X_Archive=Trial_Z_Scales;

Y_mean=mean(Y_Archive);
Total_SofS=sum((Y_Archive-Y_mean).^2);

[Nrow_1,Ncol_1]=size(X_Archive);

Basic_Length=floor(Nrow_1/10);
Remainder=Nrow_1-10*Basic_Length;
%Accuracy of the mean of the Q^2 values can be improved by increasing the number of
%iterations

Q_Squared=zeros(No_of_Iterations,1);
for Calculation=1:No_of_Iterations
    clear X Y Row_Numbers
    
        X=X_Archive;
        Y=Y_Archive;
        %Row_Numbers definition moved outside sample loop 20/10/14
        Row_Numbers=randperm(Nrow_1,Nrow_1)';
        %randperm(n,k) returns a row vector containing k unique integers selected
        %randomly from 1 to n inclusive.
        Part_RSS=zeros(10,1);
        Sample_Row_Numbers=cell(10,1);
        
        %Remainder is a number between 0 and 9. The first "Remainder" number 
        %of sets have "Basic_Length"+1 numbers in them, the remainder have Basic_Length numbers.    
        for Sample_1=1:Remainder
            Sample_Row_Numbers{Sample_1}=Row_Numbers((Sample_1-1)*(Basic_Length+1)+1:Sample_1*(Basic_Length+1),1);
        end
        for Sample_2=1:10-Remainder
            Sample_Row_Numbers{Sample_2+Remainder}=Row_Numbers(Remainder*(Basic_Length+1)+(Sample_2-1)*Basic_Length+1:Remainder*(Basic_Length+1)+(Sample_2*Basic_Length));
        end
        
            %The "Sample_Row_Numbers" cell array contains a random selection of
            %row numbers (between 1 and Nrow_1) in 10 sets- but importantly
            %none of the row numbers are duplicated across all 10 sets
            
        for Sample=1:10
            clear Y_Test_Set X_Test_Set Y_Train X_Train X_Train_Raw Y_Predict b Y_Test_Set...
            Y_Predict_Raw Coefficients Control_Variable X_ones_Train X_ones_Test...
            Sample_Row_Numbers_1 Residuals_TS
            
            Sample_Row_Numbers_1=Sample_Row_Numbers{Sample};
            Sample_Size=size(Sample_Row_Numbers_1,1);
            Counter=1;
            Control_Variable=zeros(Nrow_1,1);
            X_Test_Set=zeros(Sample_Size,Ncol_1);
            Y_Test_Set=zeros(Sample_Size,1);
            X_Train_Raw=zeros(Nrow_1-Sample_Size,Ncol_1);
            Y_Train=zeros(Nrow_1-Sample_Size,1);
            Y_Predict_Raw=zeros(Sample_Size,Ncol_1+1);
            
            for n=1:Sample_Size
                X_Test_Set(Counter,:)=X(Sample_Row_Numbers_1(n,1),:);
                Y_Test_Set(Counter,:)=Y(Sample_Row_Numbers_1(n,1),:);
                Control_Variable(Sample_Row_Numbers_1(n,1),1)=1;
                Counter=Counter+1;
            end
            Counter_1=1;
            for Case_1=1:Nrow_1 
                if Control_Variable(Case_1,1)==0
                    X_Train_Raw(Counter_1,:)=X(Case_1,:);
                    Y_Train(Counter_1,:)=Y(Case_1,:);    
                    Counter_1=Counter_1+1;
                else
                end
            end
            Nrow_3=size(X_Train_Raw,1);
            X_ones_Train=ones(Nrow_3,1);
            X_Train=[X_ones_Train,X_Train_Raw];
            
            [b,~,~,~,~]=regress(Y_Train,X_Train);
            %[~,~,~,~,BETA,~,~,~] = plsregress(X_Train,Y_Train,No_LV);
            %Coefficients=BETA';
            Coefficients=b';
            
            
            X_ones_Test=ones(Sample_Size,1);
            X_Full_Test_Set=[X_ones_Test, X_Test_Set];
            
            for nn=1:Sample_Size
                Y_Predict_Raw(nn,:)=Coefficients(1,:).*X_Full_Test_Set(nn,:);
            end
            Y_Predict=sum(Y_Predict_Raw,2);
            Residuals_TS=Y_Test_Set-Y_Predict;
            
            Part_RSS(Sample,1)=sum(Residuals_TS.^2);
        end
        RSS=sum(Part_RSS);
        Q_Squared(Calculation,1)=1-(RSS/Total_SofS); 
end
Q_Squared_Mean=mean(Q_Squared,1);



    

