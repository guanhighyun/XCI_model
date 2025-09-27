directory = 'Smoldyn_simulation'; 
mkdir(directory);

% End point of the simulation. Units in seconds.
tstop = 8000;
random_seeds = 1;
dt = [0.0002];

a_act = 1.1; % activator synthesis rate from single X chromosome
d_act = 0.72; % degradation rate of free activator
K_n = 2.2; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 9.4; % Hill coefficient for SPEN supressing activator synthesis rate. 
a_x = 3.2; % xist synthesis rate from single X chromosome. 
d_x = 0.0076; % degradation rate of Xist
m = 3.09; % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 2.06; % quantity of SPEN at which Xist dissociation is half max
K_a = 16.1; % quantity of activator at which Xist transcription is half max
P1 = [0.05]; % Probability for Xist binding to DNA
P3 = [1]; % Binding probability to chromosome for SPEN 
k2 = 8.7; % maximum dissociation rate for Xist
k4 = 10.2; % dissoication rate for bound SPEN
n_SPEN = 1000; % Total SPEN in the nucleus

rho_Xist_synthesis = [0.1]; % Reaction radius for bimolecular reactions

fid=fopen(sprintf('%s/run.sh',directory),'w');
fprintf(fid,'#!/bin/bash\n\n');
fprintf(fid,'module load python/3.8.8\n\n');
for j = random_seeds
    for h = 1:numel(P1)
        for p = 1:numel(k4)
            for t = 1:numel(a_x)
                for u = 1:numel(d_x)
                    for v = 1:numel(d_act)
                        for i = 1:numel(n_SPEN)
                            for k = 1:numel(rho_Xist_synthesis)
                                for mm = 1:numel(K_S)
                                    for nn = 1:numel(a_act)
                                        for q = 1:numel(K_n)
                                            for hh = 1:numel(P5)
                                                % How many time steps should the system record molecular coordinates
                                                samplingrate = round(0.1/dt);
                                                curr_fileprefix = sprintf('%s-P1_%g-k4_%g-P5_%g-a_x_%g-d_x_%g-d_act_%g-n_SPEN_%g-rho_Xist_%g-K_S_%g-a_act_%g-K_n_%g-%d',...
                                                filebase,P1(h),k4(p),P5(hh),a_x(t),d_x(u),d_act (v),n_SPEN(i),rho_Xist_synthesis(k),K_S(mm),a_act(nn),K_n(q),j);
                                                
                                                % Create the Smoldyn config file and bash file
                                                smoldyn_cfg(curr_fileprefix,directory,P1(h),k4(p),P5(hh),a_x(t),d_x(u),d_act(v),n_SPEN(i),rho_Xist_synthesis(k),tstop,dt,samplingrate,j);
                                                fprintf(fid,'sbatch -p general -N 1 -J Python_smoldyn -t 264:00:00 --mem=4g --wrap="python3 Execute_smoldyn.py %g %s %g %g %g %g %g %g %g %g"\n',dt,curr_fileprefix,K_S(mm),a_act(nn),K_n(q),n,k2,m,a_x,K_a);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end                               
fclose(fid);