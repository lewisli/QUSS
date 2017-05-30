function [ ParameterStruct ] = ModelParameterToStruct( ParameterNames, ...
    ParameterMatrix )
%ModelParameterToStruct Summary of this function goes here
%   Detailed explanation goes here

ParameterStruct = struct();
if (length(ParameterNames) ~= size(ParameterMatrix,2))
   display('Mismatch in parameter number'); 
   return;
end


for i = 1:length(ParameterNames)
    ParameterStruct.(ParameterNames{i}) = ParameterMatrix(:,i);
end

end

