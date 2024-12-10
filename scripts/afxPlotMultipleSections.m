function afxPlotMultipleSections(fnamesNii,threshold,dimension,sections)
    % max.wawrzyniak@medizin.uniklinik-leipzig.de
       
    if ~exist('fnamesNii','var'), fnamesNii = spm_select([1 Inf],'image'); end
    if ~exist('threshold','var'), [threshold,~] = spm_input('Threshold ',1,'r','3.1764',1); end% df = 98; p < 0.001
    if ~exist('dimension','var'), [dimension,~] = spm_input('Dimension ',2,'i','3',1); end
    if ~exist('sections','var'),  [sections,~] = spm_input('Define planes ',3,'i','40 20 0 -30',[1 Inf]); end
    
    % output dir
    [pth,~,~] = fileparts(fnamesNii(1,:));
    outDir = strcat(pth,'_sections');
    mkdir(outDir);
    
    for i = 1:size(fnamesNii,1)
        curNii = fnamesNii(i,:);
        [~, name, ~] = fileparts(curNii);
        curPNG = fullfile(outDir,strcat(name,'.png'));
        f = afxPlotSections(curNii,sections,{},[threshold Inf],[threshold Inf],1,dimension);
        set(f,'Position',[2000 100 1500/2 360/2]);
        print(f,curPNG,'-dpng','-r150');
        close(f);
    end
 end