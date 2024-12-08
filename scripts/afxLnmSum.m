function afxLnmSum(input,output,thresh)
    
     if exist(output,'file')
         fprintf('afxLnmSum: %s already exists.\n', output)
     else
        [pth,name,ext] = fileparts(output);
        if ~exist(pth,'dir'), mkdir(pth); end
        matlabbatch{1}.spm.util.imcalc.input = input;
        matlabbatch{1}.spm.util.imcalc.output = [name ext];
        matlabbatch{1}.spm.util.imcalc.outdir =  {pth};
        matlabbatch{1}.spm.util.imcalc.expression = strcat('sum(',num2str(sign(thresh)),'.*X>',num2str(thresh,10),')');
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm('Defaults','fmri');
        spm_jobman('run',matlabbatch);
     end
end