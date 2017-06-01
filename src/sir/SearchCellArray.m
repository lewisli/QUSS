function [ Index ] = SearchCellArray( String, StringArray )

for i = 1:length(StringArray)
    if(strfind(String,StringArray{i}) == 1)
        Index = i;
        return;
    end
end

Index = 0;

end

