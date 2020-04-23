%%not getting of the variables I'm supposed to

function S = replay_LoadSpikes(varargin)

fn = []; %filename

process_varargin(varargin); %allows for flexible inputs to functions

if isempty(fn);
    fn = FindFiles('*.ntt');
end

S = [];

for iF = 1:length(fn);
    [Timestamps, ScNumbers, CellNumbers, Features, Samples] =Nlx2MatSpike( fn{iF}, [1 1 1 1 1],0, 1, 1);

    nCells = max(CellNumbers);
    for iC = 1:nCells;
        temp = Timestamps(CellNumbers==iC);
        S = cat(1,S,{ts(temp./1e6)});
    end
end