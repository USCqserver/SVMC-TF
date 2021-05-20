function svmc_p2_TF_ML(step, num_samples)
%step=sweeps

%60 qubits instance
abc = load('Hamiltonian_sanity_check3_pert2_lambda1_025_train_2.mat');

abc.Hamiltonians{1}

abcH = abc.Hamiltonians{1}

abcj = abcH.Hamiltonian_J
abch = abcH.Hamiltonian_h

natom = 60;
j1 = zeros(natom, natom);
j2 = zeros(natom, natom);
j3 = zeros(natom, natom);

h1 = zeros(1, natom);
h2 = zeros(1, natom);
h3 = zeros(1, natom);

j1 = abcj{1};
j2 = abcj{2};
j3 = abcj{3};

h1 = abch{1};
h2 = abch{2};
h3 = abch{3};


% j1(0)

h = heatmap(h1);
J = heatmap(j1);
% % Define Pauli matrices
rng('shuffle')
sX = [0 1; 1 0];
sZ = [1 0; 0 -1];
unit = speye(2);

natom = 60;
updates = [1:1:natom];



beta = (1/1.57)*1e-9;

dlm = dlmread('DW_2000Q_6_annealing_schedule.txt');
slist = dlm(:,1).';
A_s = dlm(:,2).';
B_s = dlm(:,3).';
A_sp1 = @(s)interp1(slist,A_s,s);
B_sp1 = @(s)interp1(slist,B_s,s);

tau = 1e-6; 
% tf = 2*tau*(1-sstar);
tf = tau;

inslist = linspace(0,1,step+1);
length(inslist)
dt_me = tf/step;
tstep_me = 0:dt_me:tf;

tic
% Quantum trajectory by using waiting time distribution
dt_svmc = tf/step;
tstep_svmc = 0:dt_svmc:tf;

counter = 0;

num_cycles = numel(tstep_me);
num_moves = 1;
num_vars = natom;

p=2;


%new_states = 2*pi*rand(num_samples,num_cycles,num_vars);
new_states = zeros(num_samples,num_cycles,num_vars);
states = zeros(num_samples,num_cycles,num_vars);
pop_stats = zeros(num_samples,num_vars);
num_cycles 

ablist = zeros(1, numel(tstep_me));

for sample = 1:num_samples
    for index = 1:num_cycles
        for n = 1:num_vars
            A_sp1(inslist(index))/B_sp1(inslist(index));
            a = -min(1,A_sp1(inslist(index))/B_sp1(inslist(index)))*pi;
            b = min(1,A_sp1(inslist(index))/B_sp1(inslist(index)))*pi;
            r = a + (b-a) .* rand(1,1);
            new_states(sample,index,n) = r;
        end
    end
end


for sample = 1:num_samples
    state = pi/2.*ones(1,num_vars);
    for index = 1:num_cycles
        updates = updates(randperm(length(updates)));
        sum_of_sin = sum(sin(state));

        Hp = zeros(num_vars, num_vars);
        for k = 1:num_vars
            Hp(k,k) = -h1(k)*cos(state(k));
        end

        for m = 1:num_vars
            for n = m:num_vars
                Hp(m,n) = j1(m,n)*cos(state(m))*cos(state(n));
            end
        end

        Hp = sum(Hp,'all');

        energy = -1e9.*A_sp1(inslist(index)).*(2*pi/2).*(sum_of_sin) + 1e9.*B_sp1(inslist(index)).*(2*pi/2).*Hp;
        ablist(1,index) = A_sp1(inslist(index))/B_sp1(inslist(index));
        
        for spin = updates
            new_sum_of_sin = sum_of_sin-sin(state(spin))+sin(state(spin)+new_states(sample,index,spin));
            old_j = 0;
            for n = spin:num_vars
                old_j = old_j + j1(spin,n)*cos(state(spin))*cos(state(n));
            end
            new_j = 0;
            for n = spin:num_vars
                new_j = new_j + j1(spin,n)*cos(state(spin)+new_states(sample,index,spin))*cos(state(n));
            end            

            old_hj  = -h1(spin)*cos(state(spin)) + old_j;
            new_hj  = -h1(spin)*cos(state(spin)+new_states(sample,index,spin)) + new_j;
            
            new_Hp = Hp -old_hj +new_hj;

            new_energy = -1e9.*A_sp1(inslist(index)).*(2*pi/2).*(new_sum_of_sin) + 1e9.*B_sp1(inslist(index)).*(2*pi/2).*new_Hp;
            ch_energy = (new_energy - energy);
            delta = (ch_energy)*beta;
            exp(-delta);

            if ch_energy <= 0
                state(spin)= state(spin) + new_states(sample,index,spin);
                energy = new_energy;
            elseif rand <= exp(-delta)
                state(spin)= state(spin) + new_states(sample,index,spin);
                energy = new_energy;
            end

        end
        
        states(sample,index,:) = state;
    end
    
    cos(state);
    comp_state = floor(cos(state));

    for index = 1:num_vars
        if comp_state(index) == -1
            comp_state(index) = 1;
        elseif comp_state(index) == 0 | comp_state(index) == 1
            comp_state(index) = 0;
        end
    end

    comp_state;    
    pop_stats(sample,:) = comp_state;
end

eptime = toc

floor(cos(states(end,end,:)));
pop_stats

%dlmwrite('dataexport.txt',pop_stats,'delimiter','\t','newline','pc')
dlmwrite('dataexport.txt',pop_stats,'newline','pc','delimiter','')




