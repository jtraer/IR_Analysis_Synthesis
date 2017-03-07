function M=GtMtDt(Pth,Mt)
%* == GtMtDt.m i.e. Get Meta Data ==
%* Checks for Meta Data (in file specified by Pth) and if not all data is present the user is queried and the Metadata file is written

M.Path=Pth;
%* check for the existence of a text-file (Meta.txt) and read it if possible
if exist(Pth)==2;
    [tmp1,tmp2]=textread(Pth,'%s\t%s');
    %** Scroll through variables and add them to structure M if they are lsited in Mt
    for jt=1:min([length(tmp1) length(tmp2)]);
        if ~isempty(intersect(Mt,tmp1{jt}));
            eval(sprintf('M.%s=''%s'';',tmp1{jt},tmp2{jt}));
        end
    end
else
    % If no file exists preallocate a structure with a dummy variable and write the file
    fid=fopen(Pth,'w'); 
    fclose(fid);
end

%* Check the data we have (M) against the requested details (Mt) and if they don't exist query the user for them
%** Preallocate a flag of missing variables
Flg=zeros(length(Mt),1);
%** Scroll through requested details
for jp=1:length(Mt)
    %*** => specify structure and field names as strings (we will need to rename them later for nested structures)
    Mstr='M';
    Fstr=Mt{jp};
    %*** => find any .'s in the variable name => this means a nested structure
    ptndx=FndChr(Mt{jp},'.');
    %*** => if structure is nested we need to ensure the whole tree of variables exists
    if ~isempty(ptndx);
        cnt=0;
        ptndx=[0 ptndx];
        %*** => scroll through levels
        for jpt=1:(length(ptndx)-1);
            %**** ==> extract the name at this level
            Fstr_this_level=Fstr((ptndx(jpt)+1):(ptndx(jpt+1)-1));
            %**** ==> if this level doesn't exist make it.
            eval(sprintf('if ~isfield(%s,''%s''); %s.%s=[]; end',Mstr,Fstr_this_level,Mstr,Fstr_this_level))
            %**** ==> relabel the structure
            Mstr=sprintf('%s.%s',Mstr,Fstr_this_level);
            %**** ==> if this is the last level save the end variable name
            if jpt==length(ptndx)-1;
                Fstr=Fstr(ptndx(end)+1:end);
            end
        end %jpt
    end % if isempty(ptndx) 
    %Mtmp
    %tmp
    %** => See if the last field exists in the structure, and if not flag it. 
    eval(sprintf('Mvr=%s;',Mstr)); 
    if ~isfield(Mvr,Fstr); Flg(jp)=1;
        %*** => Query the missing value
        Prm=input(sprintf('What %s was %s?;\n\n',Mt{jp},Pth),'s');
        eval(sprintf('%s.%s=Prm;',Mstr,Fstr));
    else
        %*** => if value already exists print to screen and move on
        %eval(sprintf('Pvr=M.%s;',Mt{jp}));
        %fprintf('%s was %s?;\n\n',Pth,Pvr);
    end
end
%* update the text file as needed
%** See if any variables were flagged as missing
if sum(Flg)>0;
    %** Scroll through them
    for jp=find(Flg).';
        %*** Open metadata file and append the missing values to it
        fid=fopen(Pth,'a');
        fprintf(fid,'%s\t%s\n',Mt{jp},eval(sprintf('M.%s',Mt{jp})));
    end
    fprintf('Wrote %s to %s\n',Mt{jp},Pth);
    fclose(fid);
end
%* if necessary remove the dummy variable
if isfield(M,'blah')
    M=rmfield(M,'blah');
end
