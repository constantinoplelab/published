function colors = getColorScheme(datatype)
% get color scheme for visualization
% INPUT
% datatype (string): 'volume', 'side', 'delay'
% 'volume': light purple (smallest reward offer) to dark purple (largest 
% reward offer)
% 'side': orange (reward port contralateral to recording hemisphere) and
% light blue (reward port ipsilateral to recording hemisphere)
% 'delay': light grey (shortest reward delay quartile) to dark grey
% (longest reward delay quartile)
% 'block': low=blue, mixed=purple, high=red
%
% OUTPUT
% 'volume': 1x5 cell array ordered from light to dark purple
% 'side': struct 
% 'delay': 1x4 cell array ordered from pink to brown
% 'block': struct

if strcmpi(datatype, 'volume')
    colors = {'#BB7AF7', '#9A5BD0', '#793CAB', '#591E86', '#390061'};
elseif strcmpi(datatype, 'side')
    colors.contra = '#FF9C52';
    colors.ipsi = '#73ABD1';
elseif strcmpi(datatype, 'delay')
    colors = {'#E3B1C9', '#C89BA6', '#9E736F', '#603A20'};
elseif strcmpi(datatype, 'block')
    colors.low = 'b';
    colors.mixed = '#9933ff';
    colors.high = 'r';
end

