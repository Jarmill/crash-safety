classdef loc_sos_options
    %LOC_SOS_OPTIONS options data structure of pure SOS implementation of
    %system analysis. Meant only for polytopic input sets W
    %   this is a very bare-bones implementation
    %   only includes polytopic time-varying uncertainty in W: Aw >= b
    %   no switching, no box (box is a special type of polytopic)
    
    properties
        
        %% properties of run
        %terminal time   
        Tmax(1,1) double{mustBePositive}  = 5;   
        TIME_INDEP = 0; %formulation independent of time
        
        recover = 1; %recover the solution polynomials
        
        %% Variables and descriptors
        %variables defining sets (array of symbolic variables)
        
        vars = struct('t', [], 'x', [], 'w', [], 'tau', []);
        verbose = 0; %solver output: https://yalmip.github.io/faq/runinsilent/
        
        %% support sets
        X = [];         %valid set
        X_init = [];    %initial set
        X_term = [];    %terminal set                
        
        box = [];       %bounding box
        
        solver = 'mosek';
        
        t = [];
        x = [];
        
        %% Dynamics        
        f0 = []; %dynamics without uncertainty
        fw = {}; %elements are dynamics with each uncertainty
        
        %uncertainty set
        W = struct('A', [], 'b', [], 'G', []);
        
        %w | exists tau: Aw + Gtau <= b
        w = [];
        %scale time from [0, Tmax] to [0, 1]?
        scale = 1;
        
    end
    
    methods
        
        %% support set getters
        function Tsupp = get_t_supp(obj)
            %get time support
            Tsupp = struct('ineq', [], 'eq', []);
            t = obj.t;
            if obj.scale
                Tsupp.ineq = t*(1-t);
            else
                Tsupp.ineq = t*(obj.Tmax - t);
            end
        end
        
        function Xsupp = get_X(obj)
            Xsupp = obj.X;
        end                
        
        function Xsupp = get_X_init(obj)
            Xsupp = obj.X_init;
            if isempty(Xsupp)
                Xsupp = obj.get_X();
            end
        end
        
        function Xsupp = get_X_term(obj)
            Xsupp = obj.X_term;
            if isempty(Xsupp)
                Xsupp = obj.get_X();
            end
        end
        
        function allsupp = get_all_supp(obj)
            
            if obj.TIME_INDEP
                allsupp = obj.get_X();
            else
                Tsupp = obj.get_t_supp();
                Xsupp = obj.get_X();
            
                allsupp = struct('ineq', [Tsupp.ineq; Xsupp.ineq], 'eq', Xsupp.eq);
            end
        end
        
        function obj = set_box(obj, bounding_box)
            %set bounding box for X
            nx = length(obj.x);
            [box, box_center, box_half] = box_process(nx, bounding_box);
            
            X_box = struct('ineq', box_half.^2 - (obj.x - box_center).^2, 'eq', []);
            obj.X = X_box;                        
            obj.box = box;
        end
        
%         function W = bounded_noise(t, x, epsilon)
%             %produce polytopic constraints describing system with bounded
%             %noise epsilon in the H infinity sense
%             
%         end
    end
    
end
