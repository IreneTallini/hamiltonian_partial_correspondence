function Ms = load_sequence_mesh(paths, dataset_name)
% LOAD_SEQUENCE_MESH loads meshes whose path is specified in the list of string paths 
% from a single dataset specified in dataset_name 
%
%----------------------------------------------------------
% AUTHOR: Irene Tallini
%----------------------------------------------------------
%
% PLOT_FUNCTIONAL_CORRESPONDANCE(paths, dataset_name)
% INPUTS
%   - paths: list of paths to meshes
%   - dataset_name: FOR NOW ONLY TOSCA
% OUTPUT
%   - meshes with fields TRIV, VERT, m, n.

s = size(paths); s = s(2);
for i = 1:s
    tmp = load(strcat(paths(i), ".mat"));
    Ms(i) = rename_fields(tmp.surface, dataset_name);
end
end