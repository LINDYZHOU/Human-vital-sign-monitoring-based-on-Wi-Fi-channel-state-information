function output = CSIdata_in (oridata,format,parameter)
%this function is used to get CSI Mag or Phase data from PicoScenes
%document.
%writer Wang31 1142937225@qq.com
%
%oridata:PicoScenes sruct
%format:a string,including:'nonht','vht','ht','hesu','hemu'
%parameter:a string, including:'Mag','Phase'

%% collect data from thousands of struct consumes a lot of time.
%%In order to save time, the function will save the new output evertime in
%%the folder.

%safe
if nargin ~= 3
    error('输入变量数量错误!')
end
foldername=inputname(1);
filename=[format,'-',parameter];
%find whether have the folder
n=exist(foldername);
switch n
    case 0  %not have
        mkdir(foldername);%build and enter the folder
        cd(foldername);
    case 7  %already have
        cd(foldername);
        if exist([filename,'.mat'])==2  %whether have saved the data
            output=cell2mat(struct2cell(load([filename,'.mat'])));
            %return the existing data
            cd('..');
            return
        end
end
format=lower(format);%transfer into lowercase character
output_CSI=[];%initialize CSI and sequence container
output_seq=[];
%find the format and get the parameter required data
switch format
    case 'nonht'
        for i=1:length(oridata)
            data=oridata{i};
            pocketformat=getfield(getfield(data,'RxSBasic'),'PacketFormat');
            switch pocketformat
                case 0
                    CSIadd=getfield(getfield(data,'CSI'),parameter);
                    output=cat(1,output,CSIadd);
            end
        end
    case 'nt'
        for i=1:length(oridata)
            data=oridata{i};
            pocketformat=getfield(getfield(data,'RxSBasic'),'PacketFormat');
            switch pocketformat
                case 1
                    CSIadd=getfield(getfield(data,'CSI'),parameter);
                    output=cat(1,output,CSIadd);
            end
        end
    case 'vht'
        for i=1:length(oridata)
            data=oridata{i};
            pocketformat=getfield(getfield(data,'RxSBasic'),'PacketFormat');
            switch pocketformat
                case 2
                    CSIadd=getfield(getfield(data,'CSI'),parameter);
                    pocketnum = getfield(getfield(data,'StandardHeader'),'Sequence');
                    output_CSI=cat(1,output_CSI,CSIadd);
                    output_seq=cat(1,output_seq,pocketnum);
            end
        end
    case 'hesu'
        for i=1:length(oridata)
            data=oridata{i};
            pocketformat=getfield(getfield(data,'RxSBasic'),'PacketFormat');
            switch pocketformat
                case 3
                    CSIadd=getfield(getfield(data,'CSI'),parameter);
                    pocketnum = getfield(getfield(data,'StandardHeader'),'Sequence');
                    output_CSI=cat(1,output_CSI,CSIadd);
                    output_seq=cat(1,output_seq,pocketnum);
            end
        end
    case 'hemu'
        for i=1:length(oridata)
            data=oridata{i};
            pocketformat=getfield(getfield(data,'RxSBasic'),'PacketFormat');
            switch pocketformat
                case 4
                    CSIadd=getfield(getfield(data,'CSI'),parameter);
                    output=cat(1,output,CSIadd);
            end
        end
end

% output={output_CSI,output_seq};
%save
if isempty(output_CSI)
    cd('..');
    error('无该格式的包!')
else

    %interpolation
    [row,col] = size(output_CSI);
    % In PicoScense platform thesequence of package is form 0 to 4095 and
    % circulating.
    output=[];%initialize output
    output_seq=single(output_seq);%convert from uint to single
    [~,locs] = findpeaks(output_seq);% find the start and end of the circulation
    locs = [0;locs;length(output_seq)];%the first circulation is beginning from 0
    num_matrix = length(locs)-1;% the number of circulation.
    %for every circulation
    for num = 1:num_matrix
        matrix_csi=zeros(4096,col);%initialize csi container
        for i = (locs(num)+1):locs(num+1)
            matrix_csi(output_seq(i)+1,:) = output_CSI(i,:);
            %divide csi data by circulation
        end
        matrix_csi(all(matrix_csi==0,2),:)=[];%remove 0
        matrix_inter=[];%initialize interpolation container
        %for every column
        for n=1:col
            %interpolate in the position of 0 points
            intered=interp1(1+output_seq((locs(num)+1):locs(num+1)),matrix_csi(:,n),1:4096)';
            matrix_inter=cat(2,matrix_inter,intered);%conjunction
        end
        matrix_inter(any(isnan(matrix_inter),2),:)=[];%remove NaN
        output=cat(1,output,matrix_inter);
    end
    save(filename,'output');%save the data
end
cd('..');
end
