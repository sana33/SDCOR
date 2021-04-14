% dssave(fname,ds)

%

% fname = binary data set file name

% ds = matrix containing the data set

function DOLPHIN_dssave(fname,ds)

% cv=max(max(ds)-min(ds));

% mm=min(ds);

% %mm=min(min(ds));

% for i=1:d

%    ds(:,i)=(ds(:,i)-mm(i))./(1.0001*cv);

% end

ds=ds';

dim=size(ds);

file=fopen(fname,'w');

fwrite(file,[dim(1) dim(2)],'integer*4');

fwrite(file,ds,'float');

fclose(file);

