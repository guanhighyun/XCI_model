function smoldyn_cfg(fileprefix,directory,P1,k4,P3,Pa_x,d_x,d_act,n_SPEN,rho_Xist_synthesis,tstop,dt,samplingrate,random_seed)

% Configuration and output file names
cfg_name = [directory '/' fileprefix '.cfg']; 
xyz_name = [fileprefix '.txt']; 

% Start to write the configuration file
fid=fopen(cfg_name,'w');

% Model parameters
fprintf(fid, 'random_seed %d\n',random_seed);
fprintf(fid, 'variable rho_react = %g\n',0.1);
fprintf(fid, 'variable rho_Xist_synthesis = %g\n',rho_Xist_synthesis);
fprintf(fid, 'variable rho_eps = 0.000010\n');
fprintf(fid, 'variable P1 = %g\n',P1);
fprintf(fid, 'variable k4 = %g\n',k4);
fprintf(fid, 'variable P3 = %g\n',P3);
fprintf(fid, 'variable Pa_x = %g\n',Pa_x);
fprintf(fid, 'variable d_x = %g\n',d_x);
fprintf(fid, 'variable d_act = %g\n',d_act);

% Boundaries
fprintf(fid, 'dim 3\n');
fprintf(fid, 'boundaries x -8.1 8.1 r\n');
fprintf(fid, 'boundaries y -8.1 8.1 r\n');
fprintf(fid, 'boundaries z -8.1 8.1 r\n');

% Model species.
% Activator: Xist activator
% TC: Xist transcription site
% SPEN
% Xist
% Xistb: Xist bound to the chromosome
% BS: Xist binding sites on the chromosome
% SPENb: SPEN bound to the chromosome
% Act_gene: activator gene that produces activators
% XS: Xist bound with one SPEN
% XS2: Xist bound with two SPEN
% ...
% XS10: Xist bound with ten SPEN
fprintf(fid, 'species Activator TC SPEN Xist Xistb BS SPENb Act_gene XS XS1 XS2 XS3 XS4 XS5 XS6 XS7 XS8 XS9 XS10\n');

% Diffusion coefficient. Unit: um^2/min
fprintf(fid, 'difc Xist(soln) %g\n',0.4);
fprintf(fid, 'difc TC %g\n',0);
fprintf(fid, 'difc Act_gene %g\n',0);
fprintf(fid, 'difc Xistb %g\n',0);
fprintf(fid, 'difc BS %g\n',0);
fprintf(fid, 'difc SPENb %g\n',0);
fprintf(fid, 'difc Activator(soln) %g\n',15);
fprintf(fid, 'difc SPEN(soln) %g\n',15);

% Build periodic boundaries
fprintf(fid, 'start_surface outer_domain\n');
fprintf(fid, 'polygon both edge\n');
fprintf(fid, 'action both all(all) jump\n');
fprintf(fid, 'panel rect +x 8 8 -8 -16 16 r1\n');
fprintf(fid, 'panel rect +y -8 8 -8 16 16 r2\n');
fprintf(fid, 'panel rect +z -8 8 8 16 -16 r3\n');
fprintf(fid, 'panel rect -x -8 -8 -8 16 16 r4\n');
fprintf(fid, 'panel rect -y -8 -8 -8 16 16 r5\n');
fprintf(fid, 'panel rect -z -8 8 -8 16 -16 r6\n');
fprintf(fid,'jump r1 back <-> r4 back\n');
fprintf(fid,'jump r2 back <-> r5 back\n');
fprintf(fid,'jump r3 back <-> r6 back\n');
fprintf(fid, 'end_surface\n');

% Build transparent compartment of three chromosomes. All molecules can transmit
% through the boundaries. Just for definition of reactions within the
% compartment.
fprintf(fid, 'start_surface X1_boundary\n');
fprintf(fid, 'polygon both edge\n');
fprintf(fid, 'action both all(all) transmit\n');
fprintf(fid, 'panel sphere -5 0 0 2.5-0.01 10 10 panel_sphere1\n');
fprintf(fid, 'end_surface\n');

fprintf(fid, 'start_surface X2_boundary\n');
fprintf(fid, 'polygon both edge\n');
fprintf(fid, 'action both all(all) transmit\n');
fprintf(fid, 'panel sphere 0 0 0 2.5-0.01 10 10 panel_sphere2\n');
fprintf(fid, 'end_surface\n');

fprintf(fid, 'start_surface X3_boundary\n');
fprintf(fid, 'polygon both edge\n');
fprintf(fid, 'action both all(all) transmit\n');
fprintf(fid, 'panel sphere 5 0 0 2.5-0.01 10 10 panel_sphere3\n');
fprintf(fid, 'end_surface\n');

fprintf(fid, 'start_compartment chromosome_X1\n');
fprintf(fid, 'surface X1_boundary\n');
fprintf(fid, 'point -5 0 0\n');
fprintf(fid, 'end_compartment\n');

fprintf(fid, 'start_compartment chromosome_X2\n');
fprintf(fid, 'surface X2_boundary\n');
fprintf(fid, 'point 0 0 0\n');
fprintf(fid, 'end_compartment\n');

fprintf(fid, 'start_compartment chromosome_X3\n');
fprintf(fid, 'surface X3_boundary\n');
fprintf(fid, 'point 5 0 0\n');
fprintf(fid, 'end_compartment\n');

% Define the compartment of the full simulation domain
fprintf(fid, 'start_compartment full_domain\n');
fprintf(fid, 'surface outer_domain\n');
fprintf(fid, 'point rho_eps 0 0\n');
fprintf(fid, 'end_compartment\n');

% Create molecule lists
fprintf(fid, 'molecule_lists list1 list2 list3 list4 list5 list6 list7\n');
fprintf(fid, 'mol_list Activator list1\n');
fprintf(fid, 'mol_list TC list2\n');
fprintf(fid, 'mol_list SPEN list3\n');
fprintf(fid, 'mol_list Xist list4\n');
fprintf(fid, 'mol_list Xistb list5\n');
fprintf(fid, 'mol_list BS list6\n');
fprintf(fid, 'mol_list SPENb list7\n');

% Write reactions and rates.

% Activator synthesis within the whole simulation domain
fprintf(fid, 'reaction compartment=full_domain Activator_production 0 -> Activator(soln)\n');

% Activator degradation
fprintf(fid, 'reaction compartment=full_domain Activator_degredation Activator(soln) -> 0 d_act\n');

% Xist transcription from the transcription site
fprintf(fid, 'reaction XIST_transcription TC(soln) -> TC(soln) + Xist(soln)\n');
fprintf(fid, 'product_placement XIST_transcription offset TC 0 0 0\n');
fprintf(fid, 'product_placement XIST_transcription offset Xist 0.01 0 0\n'); % Place Xist 0.01 um beyond the transcription site

% Xist degradation
fprintf(fid, 'reaction XIST_degredation Xist(soln) -> 0 d_x\n');

% Xist bind to binding sites on Chromosome X1
fprintf(fid, 'reaction compartment=chromosome_X1 XIST_bindTo_chromosome1 BS(soln) + Xist(soln) -> Xistb(soln)\n');
fprintf(fid, 'reaction_probability XIST_bindTo_chromosome1 P1\n'); % probability of binding
fprintf(fid, 'binding_radius XIST_bindTo_chromosome1 rho_react\n'); % binding radius = 0.1 um

% Xist dissociate from the binding sites on Chromosome X1
fprintf(fid, 'reaction compartment=chromosome_X1 Xistb_dissociation_1  Xistb(soln) -> BS(soln) + Xist(soln)\n');
fprintf(fid, 'product_placement Xistb_dissociation_1 offset product_1 0 0 0\n');
% Place dissociated Xist 0.1 + 0.000010 um from the binding site to avoid
% repeated binding events to the same binding site
fprintf(fid, 'product_placement Xistb_dissociation_1 offset product_2 rho_react+rho_eps 0 0\n'); 

% Xist bind to Chromosome X2
fprintf(fid, 'reaction compartment=chromosome_X2 XIST_bindTo_chromosome2 BS(soln) + Xist(soln) -> Xistb(soln)\n');
fprintf(fid, 'reaction_probability XIST_bindTo_chromosome2 P1\n');
fprintf(fid, 'binding_radius XIST_bindTo_chromosome2 rho_react\n');

% Xist dissociate from Chromosome X2
fprintf(fid, 'reaction compartment=chromosome_X2 Xistb_dissociation_2  Xistb(soln) -> BS(soln) + Xist(soln)\n');
fprintf(fid, 'product_placement Xistb_dissociation_2 offset product_1 0 0 0\n');
fprintf(fid, 'product_placement Xistb_dissociation_2 offset product_2 rho_react+rho_eps 0 0\n');

% Xist bind to Chromosome X3
fprintf(fid, 'reaction compartment=chromosome_X3 XIST_bindTo_chromosome3 BS(soln) + Xist(soln) -> Xistb(soln)\n');
fprintf(fid, 'reaction_probability XIST_bindTo_chromosome3 P1\n');
fprintf(fid, 'binding_radius XIST_bindTo_chromosome3 rho_react\n');

% Xist dissociate from Chromosome X3
fprintf(fid, 'reaction compartment=chromosome_X3 Xistb_dissociation_3  Xistb(soln) -> BS(soln) + Xist(soln)\n');
fprintf(fid, 'product_placement Xistb_dissociation_3 offset product_1 0 0 0\n');
fprintf(fid, 'product_placement Xistb_dissociation_3 offset product_2 rho_react+rho_eps 0 0\n');

% One SPEN protein bind to an Xist molecule
fprintf(fid, 'reaction 1SPEN_bindTo_Xist Xistb(soln) + SPEN(soln) <-> XS(soln)\n');
fprintf(fid, 'reaction_probability 1SPEN_bindTo_Xistfwd P3\n'); % forward reaction rate
fprintf(fid, 'binding_radius 1SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 1SPEN_bindTo_Xistrev k4\n'); % reverse reaction rate
fprintf(fid, 'product_placement 1SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 1SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% A second SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 2SPEN_bindTo_Xist XS(soln) + SPEN(soln) <-> XS2(soln)\n');
fprintf(fid, 'reaction_probability 2SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 2SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 2SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 2SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 2SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% A third SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 3SPEN_bindTo_Xist XS2(soln) + SPEN(soln) <-> XS3(soln)\n');
fprintf(fid, 'reaction_probability 3SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 3SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 3SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 3SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 3SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% A fourth SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 4SPEN_bindTo_Xist XS3(soln) + SPEN(soln) <-> XS4(soln)\n');
fprintf(fid, 'reaction_probability 4SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 4SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 4SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 4SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 4SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% A fifth SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 5SPEN_bindTo_Xist XS4(soln) + SPEN(soln) <-> XS5(soln)\n');
fprintf(fid, 'reaction_probability 5SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 5SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 5SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 5SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 5SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% A sixth SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 6SPEN_bindTo_Xist XS5(soln) + SPEN(soln) <-> XS6(soln)\n');
fprintf(fid, 'reaction_probability 6SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 6SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 6SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 6SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 6SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% A seventh SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 7SPEN_bindTo_Xist XS6(soln) + SPEN(soln) <-> XS7(soln)\n');
fprintf(fid, 'reaction_probability 7SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 7SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 7SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 7SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 7SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% An eighth SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 8SPEN_bindTo_Xist XS7(soln) + SPEN(soln) <-> XS8(soln)\n');
fprintf(fid, 'reaction_probability 8SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 8SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 8SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 8SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 8SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% A ninth SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 9SPEN_bindTo_Xist XS8(soln) + SPEN(soln) <-> XS9(soln)\n');
fprintf(fid, 'reaction_probability 9SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 9SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 9SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 9SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 9SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% A tenth SPEN protein bind to the Xist molecule
fprintf(fid, 'reaction 10SPEN_bindTo_Xist XS9(soln) + SPEN(soln) <-> XS10(soln)\n');
fprintf(fid, 'reaction_probability 10SPEN_bindTo_Xistfwd P3\n');
fprintf(fid, 'binding_radius 10SPEN_bindTo_Xistfwd rho_react\n');
fprintf(fid, 'reaction_rate 10SPEN_bindTo_Xistrev k4\n');
fprintf(fid, 'product_placement 10SPEN_bindTo_Xistrev offset product_1 0 0 0\n');
fprintf(fid, 'product_placement 10SPEN_bindTo_Xistrev offset product_2 rho_react+rho_eps 0 0\n');

% Put 1000 SPEN molecules into the whold simulation domain at the initial time point
fprintf(fid, 'compartment_mol %g SPEN(soln) full_domain\n',n_SPEN);

% Chromosome centers (Xist transcription sites)
centers = [-5, 0, 0; 0, 0, 0; 5, 0, 0 ];
for i = 1:3
    center = centers(i,:);
    mat_name = sprintf('%s_chr%d.mat',cfg_name,i);
    [x,y,z] = random_walk_gene(mat_name,center); % Place Xist binding sites according to a random walk algorithm
    for j = 1:numel(x)
        fprintf(fid, 'mol 1 BS(soln) %g %g %g\n',x(j),y(j),z(j));
    end
end

fprintf(fid, 'time_start 0\n'); % Start time of the simulation
fprintf(fid, 'time_stop %g\n',tstop); % Finish time of the simulation
fprintf(fid, 'time_step %g\n',dt); % Time step

% Put Xist transcription sites onto each chromosome
fprintf(fid, 'mol 1 TC(soln) %g %g %g\n',centers(1,1),centers(1,2),centers(1,3));
fprintf(fid, 'mol 1 TC(soln) %g %g %g\n',centers(2,1),centers(2,2),centers(2,3));
fprintf(fid, 'mol 1 TC(soln) %g %g %g\n',centers(3,1),centers(3,2),centers(3,3));

% Put activator transcription sites onto each chromosome
fprintf(fid, 'mol 1 Act_gene(soln) %g %g %g\n',centers(1,1),centers(1,2),centers(1,3));
fprintf(fid, 'mol 1 Act_gene(soln) %g %g %g\n',centers(2,1),centers(2,2),centers(2,3));
fprintf(fid, 'mol 1 Act_gene(soln) %g %g %g\n',centers(3,1),centers(3,2),centers(3,3));

% Output files of molecular coordinates
fprintf(fid, 'output_files %s\n', xyz_name);
fprintf(fid, 'output_files Activator_%s\n', xyz_name);
fprintf(fid, 'cmd N %d molpos Activator(all) Activator_%s\n', samplingrate, xyz_name);

fprintf(fid, 'output_files XS_%s\n', xyz_name);
fprintf(fid, 'cmd N %d molpos XS(all) XS_%s\n', samplingrate, xyz_name);
for i = 1:10
    fprintf(fid, 'output_files XS%d_%s\n', i, xyz_name);
    fprintf(fid, 'cmd N %d molpos XS%d(all) XS%d_%s\n', samplingrate, i, i, xyz_name);
end

fprintf(fid, 'cmd N %d molpos Xist(all) %s\n', samplingrate, xyz_name);
fprintf(fid, 'cmd N %d molpos Xistb(all) %s\n', samplingrate, xyz_name);
end
