function afxLnmMean(input,output)
    
    if exist(output,'file')
        fprintf('afxLnmMean: %s already exists.\n', output)
    else
        [pth,name,ext] = fileparts(output);
        if ~exist(pth,'dir'), mkdir(pth); end
        matlabbatch{1}.spm.util.imcalc.input = input;
        matlabbatch{1}.spm.util.imcalc.output = [name ext];
        matlabbatch{1}.spm.util.imcalc.outdir =  {pth};
        matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm('Defaults','fmri');
        spm_jobman('run',matlabbatch);
    end
end