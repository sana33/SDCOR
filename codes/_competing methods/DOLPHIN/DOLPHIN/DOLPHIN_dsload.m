function [ds,rows,cols] = DOLPHIN_dsload(fname,varargin)

    if (~exist(fname,'file'))
        error('File %s not found.', fname);
    end
    
    file=fopen(fname,'r');
    dim=fread(file,[1 2],'integer*4');
    rows = dim(2);
    cols = dim(1);
    a = 1;
    b = dim(2);

    if (nargin>1)
        a = varargin{1};
        fseek(file,(a-1)*dim(1)*4,'cof');
        if (nargin>2)
            b = varargin{2};
        end
    end

    ds=fread(file,[dim(1) b-a+1],'float');
    fclose(file);

    ds=ds';
    
end

