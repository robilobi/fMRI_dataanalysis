%__________________________________________________________________________
% Batch script:
%   - second level statistics ANOVA 2 x 2 within subject design
% Modified by Roberta Bianco 2017
%__________________________________________________________________________

%===========================================================================
% user specified settings
%===========================================================================
clc
clear
clear matlabbatch
addpath('D:\MAToolBox\spm12\')


%% run up spm
spm fmri;
spm_jobman('initcfg')

%% define your variables
rootdir='D:\Roberta\SCRIPTS';
rsltdir='D:\Roberta\DATA\MRI\RFX\ANOVA\';
firstlev = 'D:\Roberta\DATA\MRI\';

  % if the intended stats directory does not exits, make it and go into it
        if exist(rsltdir) == 0
            mkdir(rsltdir);
            %cd(rsltdir);
        elseif exist(rsltdir) == 7
            % if it does exist, go into it and delete its contents.
            cd(rsltdir);
            delete('*.*');
        end

subject =  {'01','02','03'};

%% run batch
matlabbatch{1}.spm.stats.factorial_design.dir = {rsltdir};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'factor1';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'factor2';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

for s = 1:length(subject)
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans{s}= char([firstlev subject{s} '\analysis\con_0002.img,1']);
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1
                                                                    2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans{s}= char([firstlev subject{s} '\analysis\con_0003.img,1']);
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [2
                                                                    1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans{s}= char([firstlev subject{s} '\analysis\con_0006.img,1']);
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2
                                                                    2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans{s}= char([firstlev subject{s} '\analysis\con_0007.img,1']);
end 


matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'l>s';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 1 -1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 's>l';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 -1 1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 't>n';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [1 -1 1 -1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'n>t';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.convec = [-1 1 -1 1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'inter l>s';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.convec = [-1 1 1 -1];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'inter s>l';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.convec = [1 -1 -1 1];
matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
print results
matlabbatch{4}.spm.stats.results.spmmat ={[rsltdir filesep 'SPM.mat']};

for k=1:length(matlabbatch{3}.spm.stats.con.consess)
    matlabbatch{4}.spm.stats.results.conspec(k).titlestr=matlabbatch{3}.spm.stats.con.consess{k}.tcon.name;
    matlabbatch{4}.spm.stats.results.conspec(k).contrasts=k;
    matlabbatch{4}.spm.stats.results.conspec(k).threshdesc='none';
    matlabbatch{4}.spm.stats.results.conspec(k).thresh=0.00100;
    matlabbatch{4}.spm.stats.results.conspec(k).extent=5;
    matlabbatch{4}.spm.stats.results.conspec(k).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
end
matlabbatch{1}.spm.stats.results.print=0;

spm_jobman('run',matlabbatch)
