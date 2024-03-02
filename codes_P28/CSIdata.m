function output = CSIdata (oridata,format,parameter)
%this function is used to get CSI Mag or Phase data from PicoScenes
%document.
%writer Wang31 1142937225@qq.com
%
%oridata: the name of PicoScenes sruct 
%format: a string, including: 'nonht', 'vht', 'ht', 'hesu', 'hemu'
%parameter: a string, including: 'Mag', 'Phase'
%output: each row is at the same time, 
%             each col is from the same subcarrier.

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
output=[];%initialize output
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
    case 'ht'
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
                    output=cat(1,output,CSIadd);
            end
        end
    case 'hesu'
        for i=1:length(oridata)
            data=oridata{i};
            pocketformat=getfield(getfield(data,'RxSBasic'),'PacketFormat');
            switch pocketformat
                case 3
                    CSIadd=getfield(getfield(data,'CSI'),parameter);
                    output=cat(1,output,CSIadd);
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
%safe 
if length(output) == 0
    cd('..');
    error('无该格式的包!')
else
    save(filename,'output');%save the data
end
cd('..');
end
