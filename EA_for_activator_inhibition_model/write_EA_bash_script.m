% Create bash script to execute parallel EA on HPC clusters.
directory = 'EA_for_activator_inhibition_model'; mkdir(directory);
fid=fopen(sprintf('%s/run.sh',directory),'w');
fprintf(fid,'#!/bin/bash\n\n');
fprintf(fid,'module load python/3.8.8\n\n');
random_seeds = 1:500; % Number of total rounds of EA.
for j = random_seeds
    fprintf(fid,'sbatch -p general -N 1 -J Python -t 24:00:00 --mem=1g --wrap="python3 EA.py %d"\n',j);
end                               
fclose(fid);