classdef TriangleMesh < handle
    %TRIANGLEMESH Class for Triangle Meshes. 
    %   Loads triangle mesh from a given file, can perform basic
    %   operations and automatically stores computed data for later 
    %   usage (if wished). 
    % 
    %   The fields VERT, TRIV, n and m are always filled, other fields
    %   are empty by default. If fields are accessed via the corresponding
    %   function, it will check if the data already exists and make no 
    %   further computation or compute it if necessary. 
    %   Precomputed data can be safely removed by setting the field to [].
    % 
    %   written by Zorah Lähner (laehner@in.tum.de), 2016  
    
    properties
        VERT;
        TRIV;
        n;
        m;
        EVALS;
        EVECS;
        n_eigen = 0;
        LAPLACE;
        MASS;
        BARYCENTRIC;
        VORONOI;
        mass_mode = 'barycentric'; % options: barycentric, voronoi or full
        boundary_mode = 'neumann'; % options: neumann or dirichlet
        BOUNDARY;
        STIFFNESS;
        GRAD;
        DIV;
        TRI_AREAS;
        Dx; % differential, 3x2xm
        g; % first fundamental form, 2x2xm
        TRI_CENTER;
        
        X;
        Y;
        Z;
        GradOp;
        DivOp;
        grad_lbo;
        grad_lbo_RQt;
        edge_len;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONSTRUCTOR 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reads a triangle mesh out of a file. Possible types are
        % .off and .ply and tosca .mat format (reading into a surface
        % struct with TRIV, X, Y and Z.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = TriangleMesh(filename, options)
      
            if isstruct(filename) % construct TriangleMesh from struct
                obj.VERT = filename.VERT;
                obj.TRIV = filename.TRIV;
            else % read mesh from file
                % load file
                [pathstr, name, ext] = fileparts(filename);
                if strcmp('.off', ext)
                    [vert, tri] = read_off(filename);
                    obj.VERT = vert';
                    obj.TRIV = tri';
                elseif strcmp('.ply', ext)
                    [vert, tri] = read_ply(filename);
                    obj.VERT = vert;
                    obj.TRIV = tri;
                elseif strcmp('.mat', ext)
                    load(filename);
                    obj.TRIV = surface.TRIV;
                    obj.VERT = [surface.X, surface.Y, surface.Z];
                else
                    error('TriangleMesh: The file extension %s is not known.', ext);
                end
                
            end
            
            obj.n = size(obj.VERT, 1);
            obj.m = size(obj.TRIV, 1);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % READ VERTICES, X, Y, Z
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns 1st, 2nd or 3rd columns of obj.VERT, respectively.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = x(obj)
           res = obj.VERT(:,1);
        end
        
        function res = y(obj)
           res = obj.VERT(:,2);
        end
        
        function res = z(obj)
           res = obj.VERT(:,3);
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SETMASSMODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Changes the mass mode for future calculations.
        % Possible are lumped and full.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = setMassMode(obj, mode)
           if ~strcmp(mode, obj.mass_mode) 
              obj.mass(mode);
              fprintf('Mass mode is now %s. Warning: This will only affect future calculations. \n', mode)
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MASS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the mass matrix. 
        % Either the one belonging to input mode or (if no input mode
        % is given) to the obj.mass_mode setting.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = mass(obj, mode)
           if ~exist('mode', 'var')
              mode = obj.mass_mode; 
           end
            
           if (strcmp(mode, 'barycentric'))
               if isempty(obj.BARYCENTRIC)
                  obj.BARYCENTRIC = sparse(1:obj.n, 1:obj.n, vertex_voronoi(obj), obj.n, obj.n);
               end
               res = obj.BARYCENTRIC;
           elseif strcmp(mode, 'voronoi')
               if isempty(obj.VORONOI)
                   obj.VORONOI = sparse(1:obj.n, 1:obj.n, zeros(1,obj.n), obj.n, obj.n);
                   for i=1:obj.m
                      ids = obj.TRIV(i,:);
                      [~, areas] = triangle_hybridvoronoi( obj.VERT(ids, :)' );
                      for j=1:3
                          obj.VORONOI(ids(j), ids(j)) = obj.VORONOI(ids(j), ids(j)) + areas(j);
                      end
                   end
               end
               res = obj.VORONOI;
           elseif strcmp(mode, 'full')
               if isempty(obj.MASS)
                   obj.MASS = massmatrix(obj);
               end
               res = obj.MASS;
           else
               error('The mass_mode is set to incompatible type %s. Possible are barycentric, voronoi and full. \n', obj.mass_mode);
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STIFFNESS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the stiffness/cotan matrix. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = stiffness(obj, boundary_mode)
           if isempty(obj.STIFFNESS)
               obj.STIFFNESS = stiffnessmatrix(obj);
           end
           
           if exist('boundary_mode', 'var') && strcmp(boundary_mode, 'dirichlet')
               nbids = obj.nonboundary_vertices();
               obj.STIFFNESS = obj.STIFFNESS(nbids, :);
           end
           if exist('boundary_mode', 'var') && strcmp(boundary_mode, 'dirichlet_full')
               nbids = obj.nonboundary_vertices();
               bids = obj.boundary_vertices();
               S = sparse(obj.n, obj.n);
               S(bids, bids) = eye(length(bids));
               % Qui dovrebbe essere S(nbids, nbids) = obj.STIFFNESS(nbids, nbids);
               S(nbids,:) = obj.STIFFNESS(nbids, :);
               obj.STIFFNESS = S;
           end
           
           res = obj.STIFFNESS;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BOUNDARY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns the boundary part of the stiffness matrix.
        % NOTE: Not up to date.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res = boundary(obj)
           res = boundarymatrix(obj.VERT, obj.TRIV);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EIGEN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Returns n eigenvectors and corresponding eigen-
        % values of the laplacian of this mesh. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [evecs, evals] = eigen(obj, n, mass_mode, boundary_mode)
           if ~exist('boundary_mode', 'var')
             boundary_mode = obj.boundary_mode; 
           end
           
           if (n > obj.n_eigen) || (boundary_mode ~= obj.boundary_mode)
               obj.boundary_mode = boundary_mode;
               options.disp = 0;
               if (strcmp(boundary_mode, 'dirichlet'))
                   nbids = obj.nonboundary_vertices();
               else
                   nbids = 1:obj.n;
               end
               M = obj.mass(mass_mode);
               S = obj.stiffness();
               [evecs, obj.EVALS] = eigs(S(nbids, nbids), M(nbids, nbids), n, 1e-10, options);
               obj.EVALS = -diag(obj.EVALS);
               [obj.EVALS, sortIDs] = sort(obj.EVALS);
               evecs = evecs(:, sortIDs);
               
               obj.EVECS = zeros(obj.n, n);
               obj.EVECS(nbids, :) = evecs;
           end
           
           evecs = obj.EVECS;
           evals = obj.EVALS;
        end
        
        function res = evals(obj, n)
           if n > obj.n_eigen
              eigen(obj, n); 
           end
           
           res = obj.EVALS;
        end
        
        function res = evecs(obj, n)
           if n > obj.n_eigen
               eigen(obj, n);
           end
           
           res = obj.EVECS(:,1:n);
        end
        
        function res = laplacian(obj, mass_mode, mode_boundary)
           if ~exist('mode_boundary', 'var')
              mode_boundary = obj.boundary_mode; 
           end
           S = obj.stiffness(mode_boundary);
           M = obj.mass(mass_mode);
           if strcmp(mode_boundary, 'dirichlet')
               error('Dirichlet on Laplacian is deprecated.')
           end
           res = inv(M) * S;
        end
        
        %%%%%%%% GRADIENT AND DIVERGENCE %%%%%%%%
        
        function res = gradient(obj, f)
            if isempty(obj.Dx)
                obj.Dx = trimesh_differential(obj);
            end
            if isempty(obj.g)
                obj.g = trimesh_fff(obj.Dx);
            end
            res = trimesh_gradient(obj, f);
        end
        
        function res = divergence(obj, V)
            res = trimesh_divergence(obj, V);
        end
        
        function g = fff(obj)
           if isempty(obj.Dx)
                obj.Dx = trimesh_differential(obj);
           end
           if isempty(obj.g)
                obj.g = trimesh_fff(obj.Dx);
           end
           g = obj.g;
        end
        
        function Dx = differential(obj)
           if isempty(obj.Dx)
                obj.Dx = trimesh_differential(obj);
           end
           Dx = obj.Dx;
        end
        
        %%%%%%%% PROPERTIES %%%%%%%%
        function res = tri_areas(obj)
            if isempty(obj.TRI_AREAS)
                obj.TRI_AREAS = calc_tri_areas(obj);
            end
            
            res = obj.TRI_AREAS;
        end
        
        function res = tri_center(obj)
           if isempty(obj.TRI_CENTER)
              obj.TRI_CENTER = 1/3 * (obj.VERT(obj.TRIV(:,1),:) + obj.VERT(obj.TRIV(:,2),:) + obj.VERT(obj.TRIV(:,3),:));
           end
           
           res = obj.TRI_CENTER;
        end
        
        function res = boundary_vertices(obj)
            if isempty(obj.BOUNDARY)
                [~, ~, ~, obj.BOUNDARY] = get_boundary(obj.TRIV); 
                obj.BOUNDARY = sort(obj.BOUNDARY);
            end
           res = obj.BOUNDARY;
        end
        
        function res = nonboundary_vertices(obj)
            if isempty(obj.BOUNDARY)
                [~, ~, ~, obj.BOUNDARY] = get_boundary(obj.TRIV); 
            end
           res = ones(obj.n, 1);
           res(obj.BOUNDARY) = 0;
           res = logical(res);
        end
    end
    
end

