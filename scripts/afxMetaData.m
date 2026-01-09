function meta = afxMetaData()
    meta.date = datestr(now);
    if isfile("VERSION")
        meta.version = strtrim(fileread("VERSION"));
    else
        meta.version = "unknown";
    end
    meta.matlab = version();
    meta.toolboxes = ver();
end