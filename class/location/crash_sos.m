classdef crash_sos < location_crash_interface
    %CRASH_SOS worst-case crash-safety when a trajectory starts in an
    %initial set. Returns the peak-minimum control input needed to steer
    %some trajectory starting at X_init into X_term.
    

    
    methods
        function obj = crash_sos(opts)
            %CRASH_SOS Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj@location_crash_interface(opts);
        end
        
        %% set up the program
        
        function [poly_var, coeff_var] = make_poly(obj, d)
            
            %add a scalar gamma which is a lower bound on v0
            [poly_var, coeff_var] = make_poly@location_crash_interface(obj, d);
            
            gamma = sdpvar(1, 1);
            
            poly_var.gamma = gamma;
            coeff_var = [coeff_var; gamma];
            
        end
        
        %% constraints
        
        function [coeff, con, nonneg] = make_init_con(obj, d, poly)
            %Constraint on initial measure
            %auxiliary function v(t,x) lower bounded by gamma
            
            X0 = prep_space_cell(obj.opts.X_init); %not get_X_init()
 
            coeff = [];
            con = [];
            nonneg = [];

            Z = obj.opts.get_Z();
            for i = 1:length(X0)
                
                if isnumeric(X0{i})
                    v0_curr = replace(poly.v0, poly.x, X0{i});
                    XZ_curr = struct('ineq', Z, 'eq', []);
                    var_curr = poly.z;                    
                else
                    v0_curr = poly.v0;
                    XZ_curr = struct('ineq', [X0{i}.ineq; Z], 'eq', X0{i}.eq);
                    var_curr = [poly.x; poly.z];
                end
                
                [con_curr, coeff_curr] = obj.make_psatz(d, XZ_curr, v0_curr - poly.gamma, var_curr);
%                 X0_curr = [X0{i}; obj.opts.get_Z];
%                 
                
                if length(X0) == 1
                    init_tag = 'initial';
                else
                    init_tag = ['initial ', num2str(i)];
                end
                
                coeff = [coeff; coeff_curr];
                con   = [con; con_curr:init_tag];
%                 nonneg = [nonneg; init_pos];

                %the nonnegative function should hold for all time
                nonneg = [nonneg; poly.v-poly.gamma];
            end
            
        end
        
        function objective = get_objective(obj, poly_var, d)
            %fetch the SOS objective to minimize
            objective = -poly_var.gamma;
            
            if obj.opts.TIME_INDEP && (obj.opts.Tmax < Inf)
                objective = objective - obj.opts.Tmax * poly_var.alpha;
            end
        end
        
    end
end

