function experiment = readdiffdatabruker(filename_dx,filename_xml)

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

g                           = cell(nData,1);
I                           = cell(nData,1);
name                        = cell(nData,1);

for currentData = 1:nData
    dataStr                     = data(start(currentData)+numel(startStr):stop(currentData)-1);

    t1t2                        = str2num(dataStr);
    
    g{currentData}              = t1t2(:,1);
    g{currentData}              = 1e-2 * g{currentData}; % Now in T/m.
    
    I{currentData}              = t1t2(:,2);

    name{currentData}           = data(title(currentData)+9:title(currentData)+31);
end

experiment.I                = I;
experiment.g                = g;
experiment.name             = name;

%% Read experimental parameters.

xDoc                        = xmlread(filename_xml);

% Extract gradient pulse shape type.
listItem_gradShapeType      = xDoc.getElementsByTagName('gradShapeType');
thisElement_gradShapeType   = listItem_gradShapeType.item(0);
gradShapeType               = str2double(thisElement_gradShapeType.getFirstChild.getData);
switch gradShapeType
    case 0
        experiment.gradientPulseShape    = 'sine'; % Pure sine, not sine-plateau-sine.
    case 1
        experiment.gradientPulseShape    = 'trap'; % Ramp-plateau-ramp.
    case 2
        experiment.gradientPulseShape    = 'opt'; % Sine-plateau-ramp.
end

% Extract experimental parameters.
% Note that what is internally represented as delta is the true duration of
% the entire pulse, NOT the effective duration (an equivalent rectangular
% pulse). 

listItem_gradPulseTime      = xDoc.getElementsByTagName('gradPulseTime');
thisElement_gradPulseTime   = listItem_gradPulseTime.item(0);
gradPulseTime               = str2double(thisElement_gradPulseTime.getFirstChild.getData);
experiment.delta            = gradPulseTime; 
experiment.delta            = experiment.delta / 1000; % unit conversion to seconds.

listItem_DELTA              = xDoc.getElementsByTagName('DELTA');
thisElement_DELTA           = listItem_DELTA.item(0);
DELTA                       = str2double(thisElement_DELTA.getFirstChild.getData);
experiment.DELTA            = DELTA;
experiment.DELTA            = experiment.DELTA / 1000; % unit conversion to seconds.

listItem_rampUpTime         = xDoc.getElementsByTagName('rampUpTime');
thisElement_rampUpTime      = listItem_rampUpTime.item(0);
rampUpTime                  = str2double(thisElement_rampUpTime.getFirstChild.getData);
experiment.eps1             = rampUpTime;
experiment.eps1             = experiment.eps1 / 1000; % unit conversion to seconds.

listItem_rampDownTime       = xDoc.getElementsByTagName('rampDownTime');
thisElement_rampDownTime    = listItem_rampDownTime.item(0);
rampDownTime                = str2double(thisElement_rampDownTime.getFirstChild.getData);
experiment.eps2             = rampDownTime;
experiment.eps2             = experiment.eps2 / 1000; % unit conversion to seconds.

%% Compute value of independent variable b as a function of e.g. gradient.

% Note that what is internally represented as delta is the true duration of
% the entire pulse, not the effective duration (an equivalent rectangular
% pulse). Thus, Bruker TopSpin's delta variable is not the same as what is
% in here defined as delta.

GYROMAGNETICRATIO       = 2.6752e8; % rad/T/s.

b = cell(nData,1);

for currentData = 1:nData
    switch experiment.gradientPulseShape
        case 'sine'
            b{currentData}           = GYROMAGNETICRATIO^2*experiment.g{currentData}.^2*experiment.delta^2 * 4/pi^2*(experiment.DELTA-experiment.delta/4);
            
        case 'trap'            
            term        = (1/60).*((-20).*experiment.delta.^3+30.*experiment.delta.^2.*(2.*experiment.DELTA+experiment.eps1+experiment.eps2)+(-10) ...
                .*experiment.delta.*(2.*experiment.eps1.^2+3.*experiment.eps1.*experiment.eps2+2.*experiment.eps2.^2+6.*experiment.DELTA.*(experiment.eps1+ ...
                experiment.eps2))+(experiment.eps1+experiment.eps2).*(6.*experiment.eps1.^2+4.*experiment.eps1.*experiment.eps2+6.*experiment.eps2.^2+15.* ...
                experiment.DELTA.*(experiment.eps1+experiment.eps2)));
            
            b{currentData}           = GYROMAGNETICRATIO^2*experiment.g{currentData}.^2 * term;
            
        case 'opt'
            term       = (1/60).*pi.^(-3).*(20.*experiment.eps1.*(3.*experiment.delta.^2+(-3).*experiment.delta.*(2.*experiment.DELTA+ ...
                experiment.eps2)+experiment.eps2.*(3.*experiment.DELTA+experiment.eps2)).*((-2)+pi).*pi.^2+((-20).*experiment.delta.^3+ ...
                30.*experiment.delta.^2.*(2.*experiment.DELTA+experiment.eps2)+(-20).*experiment.delta.*experiment.eps2.*(3.*experiment.DELTA+experiment.eps2)+ ...
                3.*experiment.eps2.^2.*(5.*experiment.DELTA+2.*experiment.eps2)).*pi.^3+20.*experiment.eps1.^3.*((-48)+(-12).* ...
                pi+6.*pi.^2+pi.^3)+30.*experiment.eps1.^2.*pi.*(2.*experiment.DELTA.*((-2)+pi).^2+(-2).* ...
                experiment.delta.*((-8)+pi.^2)+experiment.eps2.*((-8)+pi.^2)));
            
            b{currentData}           = GYROMAGNETICRATIO^2*experiment.g{currentData}.^2 * term;
            
    end
end

experiment.b        = b;

end
