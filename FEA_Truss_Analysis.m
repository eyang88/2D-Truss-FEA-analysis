%%-------------------------------------------------------------------------
%%FEA analysis of a 2D truss system. Calculates and plots important values 
%%such as deflection, stress, strain, and safety factor. Pre-processing of 
%%loads must be done by hand prior to analysis. Supported by external .txt 
%%files that can be found in the repository.
%%-------------------------------------------------------------------------
clear; clc;

%% Main Code
[alpha,Length,I,A,E,NOD1,NOD2,NUMEL,No_Nodes,qy,qx,Py, h]= read_inputs;
[BC_Node,BC_dof,BC_magnitude,Total_imposed_BCs]= read_BCs;
[S,R] = Assemble_Global_Variables(NUMEL,A,E,I,Length,alpha,NOD1,NOD2,No_Nodes,qy,qx,Py);
[S,R] = Impose_BCs(Total_imposed_BCs,BC_magnitude,BC_dof,BC_Node,S,R);
D = Gauss_elimination(No_Nodes, R, S);
[strain, membrane_stress, shear_stress, max_bending_stress, von_mises, SF] = post_process(NUMEL,E, D, I, A, Length, NOD1, NOD2, h);

%% Reads the inputs of the supporting input geometry file
function [alpha,Length,I,A,E,NOD1,NOD2,NUMEL,No_Nodes,qy,qx,Py, h]= read_inputs
input_values = load('input_geometry.txt');
NUMEL = input_values(1,1);
No_Nodes = input_values(1,2);
for i=2:NUMEL+1
    NOD1(i-1) = input_values(i,2);
    NOD2(i-1) = input_values(i,3);
    E(i-1) = input_values(i,5);
    Length(i-1) = input_values(i,4);
    alpha(i-1) = input_values(i,6);
    A(i-1) = input_values(i,7);
    I(i-1) = input_values(i,8);
    qx(i-1) = input_values(i,9);
    qy(i-1) = input_values(i,10);
    Py(i-1) = input_values(i,11);
    h(i-1) = input_values(i,12);
end
end

%% Defines the boundary conditions through supporting .txt file
function [BC_Node,BC_dof,BC_magnitude,Total_imposed_BCs]= read_BCs
Input_BCs = load ('boundary_conditions.txt') ;
Total_imposed_BCs = Input_BCs(1,1);

for i=2:Total_imposed_BCs+1
    BC_Node(i-1) = Input_BCs(i,1);
    BC_dof(i-1) = Input_BCs(i,2);
    BC_magnitude(i-1) = Input_BCs(i,3);
end
end

%% Finds the local stiffness matrix for each element
function [SE] = get_element_stiffness(A, E, L, I, alpha)
km = A*E/L;
kb = E*I/L^3;
SE_Local = [km   0        0           -km  0         0;
            0    12*kb    6*kb*L      0    -12*kb    6*kb*L;
            0    6*kb*L   4*kb*L^2    0    -6*kb*L   2*kb*L^2;
            -km  0        0           km   0         0;
            0    -12*kb   -6*kb*L     0    12*kb     -6*kb*L;
            0    6*kb*L   2*kb*L^2    0    -6*kb*L   4*kb*L^2];
T = [cos(alpha)  sin(alpha) 0   0           0          0;
     -sin(alpha) cos(alpha) 0   0           0          0;
     0           0          1   0           0          0;
     0           0          0   cos(alpha)  sin(alpha) 0;
     0           0          0   -sin(alpha) cos(alpha) 0;
     0           0          0   0           0          1];

SE = transpose(T)*SE_Local*T;
end

%% Finds the local forces applied on each element
function [RE] = get_element_nodal_force(qy, qx, Py, L, alpha)
rm1 = qx*L/2+Py*2;
rb1 = qy*L/2;
Mb1 = qy*L^2/12;
rm2 = qx*L/2+Py*2;
rb2 = qy*L/2;
Mb2 = -qy*L^2/12;
RE_Local = [rm1; rb1; Mb1; rm2; rb2; Mb2];

T = [cos(alpha)  sin(alpha) 0   0           0          0;
     -sin(alpha) cos(alpha) 0   0           0          0;
     0           0          1   0           0          0;
     0           0          0   cos(alpha)  sin(alpha) 0;
     0           0          0   -sin(alpha) cos(alpha) 0;
     0           0          0   0           0          1];
RE = transpose(T)*RE_Local;
end

%% Compiles the global stiffness matrix and force vector using each element
function [S,R] = Assemble_Global_Variables(NUMEL,A,E,I,Length,Alpha,NOD1,NOD2,No_nodes,qy,qx,Py)
R = zeros(No_nodes*3,1);
S = zeros(No_nodes*3,No_nodes*3);
for i = 1:NUMEL
    SE = get_element_stiffness(A(i),E(i),Length(i),I(i),Alpha(i));
    RE = get_element_nodal_force(qy(i),qx(i),Length(i),Alpha(i),Py(i));
    KK(3) = 3*NOD1(i); 
    KK(2) = KK(3)-1; 
    KK(1) = KK(2)-1;
    KK(6) = 3*NOD2(i); 
    KK(5) = KK(6) - 1; 
    KK(4) = KK(5) - 1;

    for M=1:6
        R(KK(M)) = R(KK(M))+RE(M);
        
        for J=1:6
            S(KK(M),KK(J)) = S(KK(M),KK(J)) + SE(M,J); 
        end
    end    
end
end

%% Utilizes the penalty method to simulate imposing the boundary conditions
function [S,R] = Impose_BCs(Total_imposed_BCs,BC_magnitude,BC_dof,BC_Node,S,R)
penalty_number = 1E5*max(S,[],'all');
for i = 1:Total_imposed_BCs
    deflection = penalty_number*BC_magnitude(i); 
    Boundary_dof = 3*(BC_Node(i)-1)+BC_dof(i);
    R(Boundary_dof) = R(Boundary_dof) + deflection;
    S(Boundary_dof,Boundary_dof) = S(Boundary_dof,Boundary_dof)+penalty_number;
end
end
 
%% Uses the Gaussian elimination method to solve the system of equations to get the nodal deflections
function [D] = Gauss_elimination(No_Nodes, R, S)
for k=1:3*No_Nodes-1
   for i=k+1: 3*No_Nodes
       c = S(i,k)/S(k,k);
       for j=k+1:3*No_Nodes
           S(i,j) = S(i,j) - c*S(k,j);
       end
       R(i)= R(i)-c*R(k);
   end
end 
R(3*No_Nodes) = R(3*No_Nodes)/S(3*No_Nodes,3*No_Nodes);
for ii = 1:3*No_Nodes-1
    i = 3*No_Nodes-ii;
    sum = 0;
    for j = i+1:3*No_Nodes
        sum = sum+S(i,j)*R(j);
    end
    R(i) = (R(i)-sum)/S(i,i);
end
D = R;
end

%% Conducts the post-processing of the analysis, including stress, strain, and factor of safety
function [strain, membrane_stress, shear_stress, max_bending_stress, von_mises, SF] = post_process(NUMEL,E, D, I, A,L, NOD1, NOD2, h)
for i = 1:NUMEL
    strain(i) = (D(3*NOD2(i)-2) - D(3*NOD1(i)-2))/L(i);
    membrane_stress(i) = strain(i)*E(i);
    shear_stress(i) = -E(i)*I(i)/A(i)*(12/L(i)^3*D(3*NOD1(i)-1)+6/L(i)^2*D(3*NOD1(i))-12/L(i)^3*D(3*NOD2(i)-1)+6/L(i)^2*D(3*NOD2(i)));
   
    k=0;
    bending_stress = zeros(101,101);
    for x = 0:L(i)/100:L(i)
        k = k+1;
        j = 0;
        for y = -h(i)/2:h(i)/100:h(i)/2
            j = j+1;
            B1 = -6/L(i)^2+12*x/L(i)^3;
            B2 = -4/L(i)+6*x/L(i)^2;
            B3 = 6/L(i)^2-12*x/L(i)^3;
            B4 = -2/L(i)+6*x/L(i)^2;
            bending_stress(k,j) = -E(i)*y/L(i)^2*(B1*D(3*NOD1(i)-1)+B2*D(3*NOD1(i))+B3*D(3*NOD2(i)-1)+B4*D(3*NOD2(i)));
        end
    end
    max_bending_stress(i) = max(bending_stress,[],'all');
    von_mises(i) = (max_bending_stress(i)^2-2*max_bending_stress(i)*membrane_stress(i)+membrane_stress(i)^2+3*shear_stress(i)^2)^(1/2);
    SF(i) = 300E6/von_mises(i);
end
end
