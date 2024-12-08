function data = afxCsvRead(fname,delimiter)
    % Import CSV-Files
    % array = mycsvread(filename,delimiter)
    %
    % max.wawrzyniak@medizin.uniklinikum-leipzig.de
    
    if nargin < 2, delimiter = char(9); end % Tab
    data = {};
    % check if file exists
    if ~exist(fname,'file')
        disp(['File >' fname '< does not exist']);
        return;
    end
    % open csv-file
    fid = fopen(fname,'r');
    % read line by line
    tic;
    row = 1;
    while 1
        line = fgetl(fid);
        if ~ischar(line) , break, end
        if ~isempty(line)
            tmp = textscan(line,'%s','Delimiter',delimiter);
            for i = 1:length(tmp{1})
                num = str2double(tmp{1}{i});
                if ~isnan(num) && isreal(num), tmp{1}{i} = num; end
                data{row,i} = tmp{1}{i};
            end
        else
            data{row,1} = [];
        end
        row = row + 1;
    end
    t = toc;
    %disp(['Read CSV-File >' fname '< with ' num2str(size(data,1)) ' rows and ' num2str(size(data,2)) ' columns in ' num2str(round(t,3)) ' seconds.']); 
    % close csv-file
	fclose(fid);
    return;
end