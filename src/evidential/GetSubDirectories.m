function [ SubDirectories ] = GetSubDirectories( ParentDirectory )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


FileList = dir(ParentDirectory);

% Remove '.' and '..' subdirectories
FileList(1:2) = [];
dirFlags = [FileList.isdir];
SubDirectories = FileList(dirFlags);

end

