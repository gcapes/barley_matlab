function[Fitness_Result] = Multiple_Models_Fitness_Value_16AA_vs2_with_error(model, N_mutants, NAAR, Z_scales, Error)
%This program provides multiple models for the GA to be tested against.


%Fitness for the GA determined by the model for the full set of Reetz Data
%Modified 22/08/2014 so that all 8 AA sites contribute to fitness
%Model 1 to 6 have 6 interactions
%Model 7 is the "Base Model" which also has 6 interactions
%Model 8 is the "Simple Model" with no interactions
%Models 9 to 14 have 4 interactions.
%Models 15 to 20 have 2 interactions

%NOTE Models 7 and 8 have two squared correction terms. In the other models
%only the first term with a coefficient of 2.9411 is given the squared
%correction term.

clearvars -except model N_mutants NAAR Z_scales Record_of_Result Error
%MODELS 1 to 6
if model<7
Parameters=...
[1       1           1           1          1           1        %Constant
2.4      -0.20574	-0.20574	-0.721  	-0.20574	-0.20574 %I(2)
0.811485 -5.0284	2.9411      -0.721      2.9411      0.74884  %J(3)
-0.721	 2.9411     0.74884     2.4         -0.721      2.4      %K(2)
0.74884  0.74884	0.66608     0.66608     2.4 	    -0.20574 %L(3)
-0.20574 0.74884	0.66608     0.74884     0.811485	-0.721	 %M(1)
0.66608	 2.4	    2.4	        0.66608 	0.74884 	0.66608	 %N(3)
0.811485 0.66608	0.74884	    2.4         0.811485	2.9411	 %O(1)
2.9411	 -5.0284	0.811485	0.811485	2.4         -5.0284	 %P(1)
-0.721   2.4        2.4     	-0.20574	-5.0284     0.74884	 %Q(2)
0.74884	 0.66608	-5.0284     -5.0284     0.74884 	-0.721	 %R(3)
-5.0284	 0.811485	-5.0284     0.74884     0.66608     2.9411	 %S(2)
0.66608	 -0.721     2.9411      -0.20574	2.9411      0.66608	 %T(3)
-0.20574 0.811485	-0.20574	0.811485	-0.20574	0.811485 %U(1)
2.9411	 2.9411     0.811485	2.9411      -0.721      0.811485 %V(3)
2.4      -0.721     -0.721      2.9411      -5.0284     -5.0284	 %W(1)
-5.0284	 -0.20574	-0.721      -5.0284     0.66608     2.4];	 %X(1)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fitness Values are calculated in a series of steps:-
%1) Main effects
%2) Interactions
%3) Squared Term

%Part 1
X_Constant = ones(N_mutants,1);
X_Input =[X_Constant, Z_scales(:,2), Z_scales(:,6), Z_scales(:,8),...
    Z_scales(:,12), Z_scales(:,13), Z_scales(:,18), Z_scales(:,19), Z_scales(:,22),...
    Z_scales(:,26), Z_scales(:,30), Z_scales(:,32), Z_scales(:,36), Z_scales(:,37),...
    Z_scales(:,42), Z_scales(:,43), Z_scales(:,46)];

for J=1:N_mutants
    Fitness_Part_1(J,1)=sum(X_Input(J,1:17).*Parameters(:,model)');
end

%Part 2
Interactions=...
[11	5	8	14	7	9
8	3	13	4	16	15
15	8	10	4	6	16
3	9	13	12	3	9
5	7	14	6	10	12
4	8	2	11	13	3
12	6	1	8	10	4
2	10	10	11	9	11
2	6	4	5	4	10
7	5	13	14	1	2
12	10	4	5	16	13
13	12	1	1	2	15];

Interaction_Parameters=...
[0.72066	0.72066     -0.72066	-0.72066	-0.72066	-0.72066
 -1.139     1.139       -1.139      1.139       1.139       1.139
-0.6696 	-0.6696     -0.6696     0.6696      0.6696      -0.6696
0.72066     -0.72066	-0.72066	-0.72066	-0.72066	-0.72066
-1.139      1.139       1.139       1.139       1.139       -1.139
-0.6696     -0.6696     -0.6696     -0.6696     0.6696      -0.6696];


Interactions_model=Interactions(:,model);
X_Interactions=...
    [X_Input(:,Interactions_model(1,1)+1).*X_Input(:,Interactions_model(2,1)+1),...
    X_Input(:,Interactions_model(3,1)+1).*X_Input(:,Interactions_model(4,1)+1),...
    X_Input(:,Interactions_model(5,1)+1).*X_Input(:,Interactions_model(6,1)+1),...
    X_Input(:,Interactions_model(7,1)+1).*X_Input(:,Interactions_model(8,1)+1),...
    X_Input(:,Interactions_model(9,1)+1).*X_Input(:,Interactions_model(10,1)+1),...
    X_Input(:,Interactions_model(11,1)+1).*X_Input(:,Interactions_model(12,1)+1)];
%X_Interactions=X_Interactions_dash';

Model_Interaction_Parameters_dash=Interaction_Parameters(:,model);
Model_Interaction_Parameters=Model_Interaction_Parameters_dash';
for K=1:N_mutants
    Fitness_Vector=Model_Interaction_Parameters.*X_Interactions(K,:); 
    Fitness_Part_2(K,1)= sum(Fitness_Vector);
end

%Part 3
for ij=2:17
    if Parameters(ij,model)==2.9411;
        Squared_Term=ij;
        break
    end
end
Squared_Z_Scales=X_Input(:,Squared_Term).*X_Input(:,Squared_Term);
for L=1:N_mutants
    Fitness_Part_3(L,1)=-0.24995*Squared_Z_Scales(L,1);
end

Fitness = Fitness_Part_1 + Fitness_Part_2 + Fitness_Part_3;

%Fitness of Optimal Sequence		
%Model	Fitness     Sequence
%1      275.7971	
%2      147.1037    
%3      161.415
%4      150.7107
%5      180.8312
%6      103.9629

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%Model 7 is the "Base Model" and Model 8 is a version with no interactions
elseif model<9
    
Parameters_2=...
[1       1          %Constant
2.4      0.66608    %I(2)
2.9411	 -0.721     %J(3)
-5.0284	 2.9411     %K(2)
0.74884	 -5.0284    %L(3)
-0.20574 -0.20574	%M(1)
0.811485 0.66608	%N(3)
0.66608	 0.811485	%O(1)
-0.721	 2.4        %P(1)
2.4      2.9411	    %Q(2)
2.9411	 0.74884	%R(3)
-5.0284	 0.74884	%S(2)
0.74884	 -0.20574	%T(3)
-0.20574 -5.0284	%U(1)
0.811485 -0.721     %V(3)
0.66608	 2.4        %W(1)
-0.721	 0.811485];	%X(1)

%Part 1
X_Constant = ones(N_mutants,1);
X_Input =[X_Constant, Z_scales(:,2), Z_scales(:,6), Z_scales(:,8),...
    Z_scales(:,12), Z_scales(:,13), Z_scales(:,18), Z_scales(:,19), Z_scales(:,22),...
    Z_scales(:,26), Z_scales(:,30), Z_scales(:,32), Z_scales(:,36), Z_scales(:,37),...
    Z_scales(:,42), Z_scales(:,43), Z_scales(:,46)];

for J=1:N_mutants
    Fitness_Part_1(J,1)=sum(X_Input(J,1:17).*Parameters_2(:,model-6)');
end

%Part 2
if model==7
Interactions=...
    [1
    2
    1
    3
    4
    5
    9
    10
    9
    11
    12
    13];

    Interaction_Parameters_2=...
    [-0.72066
    1.139
    -0.6696
    -0.72066
    1.139
    -0.6696];



    X_Interactions=...
    [X_Input(:,Interactions(1,1)+1).*X_Input(:,Interactions(2,1)+1),...
    X_Input(:,Interactions(3,1)+1).*X_Input(:,Interactions(4,1)+1),...
    X_Input(:,Interactions(5,1)+1).*X_Input(:,Interactions(6,1)+1),...
    X_Input(:,Interactions(7,1)+1).*X_Input(:,Interactions(8,1)+1),...
    X_Input(:,Interactions(9,1)+1).*X_Input(:,Interactions(10,1)+1),...
    X_Input(:,Interactions(11,1)+1).*X_Input(:,Interactions(12,1)+1)];
    %X_Interactions=X_Interactions_dash';


    Model_Interaction_Parameters=Interaction_Parameters_2';
    for K=1:N_mutants
        Fitness_Vector=Model_Interaction_Parameters.*X_Interactions(K,:); 
        Fitness_Part_2(K,1)= sum(Fitness_Vector);
    end
    %Part 3- For the "Base Model".
    Squared_Z_Scales_1=X_Input(:,3).*X_Input(:,3);
    Squared_Z_Scales_2=X_Input(:,11).*X_Input(:,11);
    for L=1:N_mutants
        Fitness_Part_3(L,1)=-0.24995*Squared_Z_Scales_1(L,1)-0.24995*Squared_Z_Scales_2(L,1);
    end
    
    
elseif model==8
    Fitness_Part_2=zeros(N_mutants,1);
    
    %Part 3
    Squared_Z_Scales_1=X_Input(:,4).*X_Input(:,4);
    Squared_Z_Scales_2=X_Input(:,10).*X_Input(:,10);
    for L=1:N_mutants
        Fitness_Part_3(L,1)=-0.24995*Squared_Z_Scales_1(L,1)-0.24995*Squared_Z_Scales_2(L,1);
    end
end


Fitness = Fitness_Part_1 + Fitness_Part_2 + Fitness_Part_3;

%Model	Fitness	    Sequence
%7  	253.5575	GCGCFCDFGCGCFCDF
%8	    103.6071	RRRRFCDDRCRRFRDD

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%Models 9 to 14 inclusive:- randomized coefficients with 4 interactions
elseif model<15

Parameters_3=...
[1       1        1        1        1        1          %Constant
0.811485 0.74884  2.9411   -0.20574 0.74884	 0.811485	%I(2)
-0.20574 2.9411	  2.4	   2.4	    0.66608	 -0.20574	%J(3)
0.66608	 -0.721	  -5.0284  0.66608  -0.20574 2.9411	    %K(2)
0.66608	 2.9411	  -0.20574 -5.0284  0.66608	 2.4	    %L(3)
2.4	     -5.0284  0.66608  -0.721   0.811485 0.811485	%M(1)
-0.721	 -0.721	  0.66608  0.74884  2.4	     0.66608	%N(3)
0.74884	 0.66608  2.9411   -0.20574 2.4	     -5.0284	%O(1)
2.9411	 0.66608  -0.721   0.811485 0.811485 -0.721	    %P(1)
-0.20574 0.74884  2.4	   2.9411   2.9411	 0.66608	%Q(2)
-0.721	 -0.20574 0.811485 2.9411	-0.20574 2.9411	    %R(3)
-5.0284	 0.811485 -0.721   0.66608	-0.721	 -0.721     %S(2)
-5.0284	 -0.20574 0.811485 -0.721	0.74884	 -0.20574	%T(3)
0.811485 0.811485 -0.20574 0.74884	-5.0284	 -5.0284	%U(1)
0.74884	 2.4	  -5.0284  2.4	    -5.0284	 2.4        %V(3)
2.9411	 -5.0284  0.74884  0.811485	-0.721	 0.74884	%W(1)
2.4	     2.4	  0.74884  -5.0284	2.9411	 0.74884];	%X(1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Part 1
X_Constant = ones(N_mutants,1);
X_Input =[X_Constant, Z_scales(:,2), Z_scales(:,6), Z_scales(:,8),...
    Z_scales(:,12), Z_scales(:,13), Z_scales(:,18), Z_scales(:,19), Z_scales(:,22),...
    Z_scales(:,26), Z_scales(:,30), Z_scales(:,32), Z_scales(:,36), Z_scales(:,37),...
    Z_scales(:,42), Z_scales(:,43), Z_scales(:,46)];

for J=1:N_mutants
    Fitness_Part_1(J,1)=sum(X_Input(J,1:17).*Parameters_3(:,model-8)');
end

%Part 2
Interactions_3=...
[5	11	7	11	5	3
4	1	15	16	15	10
9	10	14	6	7	5
13	13	2	9	3	2
1	4	16	13	3	14
12	7	6	6	15	11
6	11	11	8	16	2
7	3	15	7	10	13];

Interaction_Parameters_3=...
[0.72066	0.72066	-0.72066	-0.72066	-0.72066	-0.72066
-1.139      1.139	-1.139      1.139       1.139       1.139
-0.6696     -0.6696	-0.6696     0.6696      0.6696      -0.6696
-0.6696     0.6696	0.72066     0.6696      -1.139      -0.72066];


Interactions_model=Interactions_3(:,model-8);
X_Interactions=...
    [X_Input(:,Interactions_model(1,1)+1).*X_Input(:,Interactions_model(2,1)+1),...
    X_Input(:,Interactions_model(3,1)+1).*X_Input(:,Interactions_model(4,1)+1),...
    X_Input(:,Interactions_model(5,1)+1).*X_Input(:,Interactions_model(6,1)+1),...
    X_Input(:,Interactions_model(7,1)+1).*X_Input(:,Interactions_model(8,1)+1)];
%X_Interactions=X_Interactions_dash';

Model_Interaction_Parameters_dash=Interaction_Parameters_3(:,model-8);
Model_Interaction_Parameters=Model_Interaction_Parameters_dash';
for K=1:N_mutants
    Fitness_Vector=Model_Interaction_Parameters.*X_Interactions(K,:); 
    Fitness_Part_2(K,1)= sum(Fitness_Vector);
end

%Part 3
for ij=2:17
    if Parameters_3(ij,model-8)==2.9411;
        Squared_Term=ij;
        break
    end
end
Squared_Z_Scales=X_Input(:,Squared_Term).*X_Input(:,Squared_Term);
for L=1:N_mutants
    Fitness_Part_3(L,1)=-0.24995*Squared_Z_Scales(L,1);
end

Fitness = Fitness_Part_1 + Fitness_Part_2 + Fitness_Part_3;

%Model      Fitness     Sequence
%9          151.5645	RRRCDRDDGRGRDCDD
%10         156.9155	GHGCFRFDRCGRDCFD
%11         169.4052	RCGRDCDFRCGCFRFF
%12         177.2177	GCRRFCDDRCRRDCDF
%13         150.3016	RCRCFCDDRRGCFRDD Note GA gets trapped with 13 
%14         194.7765	RCVCDCFFRCGRFCDD


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%Models 15 to 20 inclusive:- randomized coefficients with 2 interactions
elseif model<21

Parameters_4=...
[1       1          1        1          1           1           %Constant
0.74884	 0.66608	-0.721	 2.4        0.811485	-0.20574	%I(2)
0.66608	 -0.721     2.4      -0.721     0.74884     0.74884     %J(3)
-0.20574 2.4        0.811485 -0.721     0.74884     0.66608     %K(2)
0.66608	 -0.20574	-0.20574 0.66608	2.9411      2.9411      %L(3)
-5.0284	 0.74884	0.811485 -5.0284	-0.721      -5.0284     %M(1)
-0.721	 0.66608	0.66608	 0.74884	0.811485	2.9411      %N(3)
0.74884	 -0.721     -5.0284	 -0.20574	2.4         2.4         %O(1)
0.811485 0.74884	-0.721	 -0.20574	-5.0284     0.811485	%P(1)
2.9411	 -5.0284	0.66608	 2.4        -5.0284     -0.721      %Q(2)
2.9411	 2.4        2.4	     0.74884	-0.721      0.74884     %R(3)
-0.20574 -0.20574	0.74884	 0.811485	0.66608     -5.0284     %S(2)
-5.0284	 2.9411     -0.20574 0.66608	-0.20574	-0.721      %T(3)
2.4	     -5.0284	0.74884	 2.9411     2.9411      2.4         %U(1)
0.811485 0.811485	-5.0284	 2.9411     2.4         0.811485	%V(3)
-0.721	 2.9411     2.9411	 -5.0284	0.66608     -0.20574	%W(1)
2.4	     0.811485	2.9411	 0.811485	-0.20574	0.66608];	%X(1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Part 1
X_Constant = ones(N_mutants,1);
X_Input =[X_Constant, Z_scales(:,2), Z_scales(:,6), Z_scales(:,8),...
    Z_scales(:,12), Z_scales(:,13), Z_scales(:,18), Z_scales(:,19), Z_scales(:,22),...
    Z_scales(:,26), Z_scales(:,30), Z_scales(:,32), Z_scales(:,36), Z_scales(:,37),...
    Z_scales(:,42), Z_scales(:,43), Z_scales(:,46)];

for J=1:N_mutants
    Fitness_Part_1(J,1)=sum(X_Input(J,1:17).*Parameters_4(:,model-14)');
end

%Part 2
Interactions_4=...
[9	6	4	7	6	4
3	9	2	1	2	7
14	7	3	15	13	2
10	2	4	16	7	3];

Interaction_Parameters_4=...
[0.72066	0.72066	 1.139	0.72066     -0.72066	0.6696
-1.139      -0.6696	 0.6696	-0.72066	0.72066     1.139];


Interactions_model=Interactions_4(:,model-14);
X_Interactions=...
    [X_Input(:,Interactions_model(1,1)+1).*X_Input(:,Interactions_model(2,1)+1),...
    X_Input(:,Interactions_model(3,1)+1).*X_Input(:,Interactions_model(4,1)+1)];
%X_Interactions=X_Interactions_dash';

Model_Interaction_Parameters_dash=Interaction_Parameters_4(:,model-14);
Model_Interaction_Parameters=Model_Interaction_Parameters_dash';
for K=1:N_mutants
    Fitness_Vector=Model_Interaction_Parameters.*X_Interactions(K,:); 
    Fitness_Part_2(K,1)= sum(Fitness_Vector);
end

%Part 3
for ij=2:17
    if Parameters_4(ij,model-14)==2.9411;
        Squared_Term=ij;
        break
    end
end
Squared_Z_Scales=X_Input(:,Squared_Term).*X_Input(:,Squared_Term);
for L=1:N_mutants
    Fitness_Part_3(L,1)=-0.24995*Squared_Z_Scales(L,1);
end

Fitness = Fitness_Part_1 + Fitness_Part_2 + Fitness_Part_3;
end

%Model      Fitness         Sequence
%15         127.3676	RCRCFRDDRCGRDRFD
%16         125.6362	RRRRDRFDGCGCFCDD
%17         135.5818	GCRCDCFFRCRRDRSD
%18         119.3987	RRGCFCDFRCRCDCFD
%19         145.8370	RCRCFRDFGRRRDCDF
%20         178.1557	GCRCFCDDGCGRDCFD

for J=1:N_mutants
    %Error is of the form 0.01 (=1%), 0.02 (=2%) etc. randn returns random
    %numbers from the normal distribution with mean=0 and variance=1.
    Fitness_with_error = Fitness(J,1)*(1+Error*randn);
    Fitness_Result(J,:)=[Fitness(J,1), Fitness_with_error];
end