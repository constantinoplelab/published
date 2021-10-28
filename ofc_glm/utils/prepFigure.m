function f = prepFigure(ops,fignum,varargin)
    %helper function to make creation of figures easier
    %handles visibility, figure number, and lists of figures. 
    %this is just a lot of code to keep writing
    
    if isfield(ops,'fignum')
        assert( ~isempty(fignum) | fignum~=0 ,...
            'figure numbering should be supplied if ops.fignum exists');
    end
    
        
                
    if isfield(ops,'fignum')
        if isfield(ops,'vishandle')
            if ~ops.vishandle
                f = figure('Visible','off');
            else
                f = figure(fignum);
            end
        else
            f = figure(fignum);
        end
    else
        if isfield(ops,'vishandle')
            if ~ops.vishandle
                f = figure('Visible','off');
            else
                f = figure('Visible','on');
            end
        else
            f = figure('Visible','on');
        end
    end