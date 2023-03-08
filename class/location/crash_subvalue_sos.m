classdef crash_subvalue_sos < location_crash_interface
    %CRASH_SUBVALUE_SOS create a subvalue map, returning the a lower bound
    %on the minimum control effort needed to crash into the unsafe set at
    %any given initial condition 
    

    
    methods
        function obj = crash_subvalue_sos(opts)
            %CRASH_SUBVALUE_SOS Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj@location_crash_interface(opts);
        end
        
        %% set up the program
        
        function [poly_var, coeff_var] = make_poly(obj, d)
            
            %add a polynomial q(x) which is a lower bound on v(0, x, z)
            [poly_var, coeff_var] = make_poly@location_crash_interface(obj, d);
            
%             gamma = sdpvar(1, 1);
            [q, cq] = polynomial(poly_var.x, d);
            
            poly_var.q= q;
            poly_var.cq = cq;
            coeff_var = [coeff_var; cq];
            
        end
        
        %% constraints
        
        function [coeff, con, nonneg] = make_init_con(obj, d, poly)
            %Constraint on initial measure
            %auxiliary function v(0,x,z) lower bounded by q(x)
            
            X0 = obj.prep_space_cell(obj.opts.X_init); %not get_X_init()
 
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
                
                [con_curr, coeff_curr] = obj.make_psatz(d, XZ_curr, v0_curr - poly.q, var_curr);
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
                nonneg = [nonneg; poly.v-poly.q];
            end
            
        end
        
        function objective = get_objective(obj, poly_var, d)
            %fetch the SOS objective to minimize
%             objective = -poly_var.gamma;
            mom_leb = obj.opts.mom_handle(d);            
            objective = -poly_var.cq'*mom_leb;
            
            if obj.opts.TIME_INDEP && (obj.opts.Tmax < Inf)
                objective = objective - obj.opts.Tmax * poly_var.alpha;
            end
        end
        
        function [poly_eval, func_eval] = recover_poly(obj, poly_var)
            
            [poly_eval, func_eval] = recover_poly@location_crash_interface(obj, poly_var);
            
            x = poly_var.x;
            [cq,mq] = coefficients(poly_var.q,poly_var.x);
            q_eval = value(cq)'*mq;  
            
            poly_eval.q = q_eval;
            func_eval.q = polyval_func(q_eval, x);
            
        end
        
    end
end

