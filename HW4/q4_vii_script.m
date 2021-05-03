%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT: HW4: Q2.vii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Branch information & constants
% fbus tbus r x b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_BUS = 1;
T_BUS = 2;
BR_R = 3;
BR_X = 4;
BR_B = 5;
branch = [1 2 0.01 0.05 0.02*2];
% branch = [1 2 0.004 0.0533 0;
%     2 3 0.02 0.25 0.22;
%     3 4 0.02 0.25 0.22;
%     2 4 0.01 0.15 0.11;
%     4 5 0.006 0.08 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bus information & constants
% bus_no type pd qd vm va
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BUS_I = 1;
BUS_TYPE = 2;
PQ = 1;
REF = 3;
PD = 3;
QD = 4;
VM = 5;
VA = 6;
bus = [2 3 0 0 1 0.524;
       1 1 2 0.75 1 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_branch = size(branch, 1);
n_bus = max(max(branch(:, F_BUS)), max(branch(:, T_BUS)));
Y=zeros(n_bus,n_bus);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y-bus formation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_s = 1./(branch(:, BR_R) + 1j * branch(:, BR_X));
B_c = branch(:, BR_B);
Y_tt = Y_s + 1j*B_c/2;
for br_idx=1:1:n_branch
    fb = branch(br_idx, F_BUS);
    tb = branch(br_idx, T_BUS);
    Y(fb,tb)=Y(fb,tb)-Y_s(br_idx);
    Y(tb,fb)=Y(fb,tb);
end
for bus_idx=1:n_bus
    for br_idx=1:n_branch
        if branch(br_idx, F_BUS)== bus_idx 
            Y(bus_idx,bus_idx)=Y(bus_idx,bus_idx)+Y_tt(br_idx);
        elseif branch(br_idx, T_BUS)== bus_idx
            Y(bus_idx,bus_idx)=Y(bus_idx,bus_idx)+Y_tt(br_idx);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NR method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = bus(:,VM);
Theta = bus(:,VA);
P = bus(:, PD);
Q = bus(:, QD);
G = real(Y);
B = imag(Y);
PQ_bus_index = find(bus(:,BUS_TYPE)==PQ);
n_pq_bus = length(PQ_bus_index);
n_iteration = 5;
iter_idx = 0;
for itr_idx=1:1:n_iteration
    % Calculate power based on current values of V, Theta
    [P_val,Q_val]=get_power_val(G,B,n_bus,Theta,V);
    % Find how off the power values are from the given power values
    [mismatch_PQ]=get_power_mismatch(P,Q,P_val,Q_val,...
        n_pq_bus,PQ_bus_index);
    % Form Jacobian Matrix based on current V, Theta
    [J]=get_jacobian(V,Theta,G,B,n_bus,PQ_bus_index,n_pq_bus);
    % Find Delta V and Delta Theta
    d_theta_v=inv(J)*mismatch_PQ;
    d_theta=d_theta_v(1:n_bus-1);
    d_v=d_theta_v(n_bus:end);
    % Update Theta
    Theta(2:n_bus)=Theta(2:n_bus)+d_theta;
    % Update V
    for i=1:n_pq_bus
        j=PQ_bus_index(i);
        V(j)=V(j)+d_v(i);
    end
end
disp('After 5 NR iterations...');
fprintf('Bus 1 voltage = %f\n', V(2));
fprintf('Bus 1 angle = %f\n', Theta(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q4 approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_bar = Y(2,1);
Y_mod = Y(2,2);
V_nom = -inv(Y_mod)*Y_bar*V(2)*(cos(Theta(2)) + 1j*sin(Theta(2)));
V_nom_mag = abs(V_nom);
Theta_nom = angle(V_nom);
K = diag(V_nom)*conj(Y_mod);
J = [real(K) imag(K);
    imag(K) real(K)];
delta_V = inv(J)*[P(1);Q(1)];
Theta_approx = Theta_nom + inv(diag(V_nom_mag))...
    *[-diag(sin(Theta_nom)) diag(cos(Theta_nom))]*delta_V;
V_approx = V_nom_mag...
    + [diag(cos(Theta_nom)) diag(sin(Theta_nom))]*delta_V;
fprintf('\n');
disp('Approximation result...');
fprintf('Bus 1 voltage = %f\n', V_approx);
fprintf('Bus 1 angle = %f\n', Theta_approx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definition: Power calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P_val, Q_val]=get_power_val(G,B,n_bus,Theta,V)
P_val=zeros(n_bus,1);
Q_val=zeros(n_bus,1);
for i = 1:n_bus
    for j = 1:n_bus
        P_val(i) = P_val(i) + V(i)*V(j)*( G(i,j)*cos(Theta(i)-Theta(j))...
            + B(i,j)*sin(Theta(i)-Theta(j)));
        Q_val(i) = Q_val(i) + V(i)*V(j)*( G(i,j)*sin(Theta(i)-Theta(j))...
            - B(i,j)*cos(Theta(i)-Theta(j)));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definition: Get power diff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mismatch_PQ]=get_power_mismatch(P_act,Q_act,P_val,...
    Q_val,n_pq_bus,pq_bus_idx)
    dP=P_act-P_val;
    dQ=Q_act-Q_val;
    mismatch_P=zeros(n_pq_bus,1);
    mismatch_Q=zeros(n_pq_bus,1);
    for k=1:1:n_pq_bus
        n=pq_bus_idx(k);
        mismatch_Q(k, 1)=dQ(n);
        mismatch_P(k, 1)=dP(n);
    end
    mismatch_PQ=[mismatch_P;mismatch_Q];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definition: Get Jacobian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J]=get_jacobian(V,Theta, G, B, n_bus, PQ_bus_idx, n_pq_bus)
% Build J11, J12, J21, J22
J11 = zeros(n_bus-1,n_bus-1);
J12 = zeros(n_bus-1,n_pq_bus);
J21 = zeros(n_pq_bus,n_bus-1);
J22 = zeros(n_pq_bus,n_pq_bus);
%%
for k=2:n_bus
    for j=2:n_bus
        if (j==k) % Diagonal elements
        for m=1:n_bus
            J11(k-1,j-1)=J11(k-1,j-1)+(V(k)*V(m))...
                *(-G(k,m)*sin(Theta(k)-Theta(m))...
                +B(k,m)*cos(Theta(k)-Theta(m)));
        end
            J11(k-1,j-1)=J11(k-1,j-1)...
                -(V(k)^2)*B(k,k);
        else % Non digonal element
            J11(k-1,j-1)=V(k)*V(j)*(G(k,j)*sin(Theta(k)-Theta(j))...
                -B(k,j)*cos(Theta(k)-Theta(j)));
        end
    end
end
%%
for k=2:n_bus
    for j=1:n_pq_bus
        m=PQ_bus_idx(j);
        if (m==k) % Diagonal element
        for m=1:n_bus
            J12(k-1,j)=J12(k-1,j)+V(m)*(G(k,m)*cos(Theta(k)-Theta(m))...
                +B(k,m)*sin(Theta(k)-Theta(m)));
        end
            J12(k-1,j)=J12(k-1,j)+(V(k))*G(k,k);
        else % Non diagonal element
            J12(k-1,j)=V(k)*(G(k,m)*cos(Theta(k)-Theta(m))...
                +B(k,m)*sin(Theta(k)-Theta(m)));
        end
    end
end
%%
for k=1:n_pq_bus
    n=PQ_bus_idx(k);
    for j=2:n_bus
        if (j==n) % Diagonal element
        for m=1:n_bus
            J21(k,j-1)=J21(k,j-1)+(V(n)*V(m))...
                *(G(n,m)*cos(Theta(n)-Theta(m))...
                +B(n,m)*sin(Theta(n)-Theta(m)));
        end
            J21(k,j-1)=J21(k,j-1)-(V(n)^2)*G(n,n);
        else % Non diagonal element
            J21(k,j-1)=V(n)*V(j)*(-G(n,j)*cos(Theta(n)-Theta(j))...
                -B(n,j)*sin(Theta(n)-Theta(j)));
        end
    end
end
%%
for k=1:n_pq_bus
    n=PQ_bus_idx(k);
    for j=1:n_pq_bus
        m=PQ_bus_idx(j);
        if (m==n) % Diagonal element
        for m=1:n_bus
            J22(k,j) = J22(k,j) + V(m)*(G(n,m)*sin(Theta(n)-Theta(m))...
                - B(n,m)*cos(Theta(n)-Theta(m)));
        end
            J22(k,j)=J22(k,j)-V(n)*B(n,n);
        else % Non diagonal element
            J22(k,j)=V(n)*(G(n,m)*sin(Theta(n)-Theta(m))...
                - B(n,m)*cos(Theta(n)-Theta(m)));
        end
    end
end
%% Form Jacobian
J=[J11 J12;J21 J22];
end

