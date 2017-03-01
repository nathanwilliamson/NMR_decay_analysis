function experiment = readT1data(filename_dx)

experiment                  = struct();

%% Read measurements.

data                        = fileread(filename_dx);

startStr                    = '##XYPOINTS= (XY..XY)';
start                       = strfind(data,startStr);

stopStr                     = '##END= $$ End of List for Peak';
stop                        = strfind(data,stopStr);

titleStr                    = '##TITLE= Peak'; 
title                       = strfind(data,titleStr);

nData                       = numel(start);

t                           = cell(nData,1);
I                           = cell(nData,1);
name                        = cell(nData,1);

for currentData = 1:nData
    dataStr                     = data(start(currentData)+numel(startStr):stop(currentData)-1);

    t1t2                        = str2num(dataStr);
    
    t{currentData}              = t1t2(:,1);
    
    I{currentData}              = t1t2(:,2);

    name{currentData}           = data(title(currentData)+9:title(currentData)+31);
end

experiment.I                = I;
experiment.t                = t;
experiment.name             = name;


end
