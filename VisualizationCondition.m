%%%%%%%%%%%%%% for visualization

switch Dir
    case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
 case 'smallSync'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
    case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
  case 'largeSync'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
    case 'vs->vta'
       switch binwidthType
           case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
           case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
           case 'all'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.6]');
       end
    case 'vta->vs'
        switch binwidthType
           case 'small'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.25]');
           case 'large'
BinDirText=strcat(Dir,'--BinSize Range: [0.35, 0.6]');
            case 'all'
BinDirText=strcat(Dir,'--BinSize Range: [0.01, 0.6]');
        end
    otherwise
        BinDirText=strcat('Region--',Region,'--All Bin');
end


switch TrialsT
    case 'Hit'
        TrialsT = strcat('Hit');
    case 'Rej'
        TrialsT = strcat('Rej');
    case 'False Alarm'
        TrialsT =strcat('False Alarm');
    case 'Miss'
        TrialsT =strcat('Miss');
end


if (ismember(start,startOr) & ismember(stop,stopOr))
    PhaseStr = strcat('Original');
elseif (ismember(start,startRev) & ismember(stop,stopRev))
    PhaseStr = strcat('Reversal');
elseif (ismember(start,startOr) & ismember(stop,stopRev))
    PhaseStr = strcat('Or and Rev');   
end