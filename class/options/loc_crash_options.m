classdef loc_crash_options < loc_sos_options
    %LOC_CRASH_OPTIONS SOS options for the crash-safety problem
    
    properties
        %dominant change: add the 'z' variable to perform the peak
        %minimizing control.
        Zmax = 1;
        z = [];
        w = [];
    end
    
    methods
        function obj = loc_crash_options()
            %LOC_CRASH_OPTIONS 
            obj@loc_sos_options();
            
            obj.z = sdpvar(1,1);
        end
        
        function Xsupp = get_X(obj)
            Xsupp = get_X@loc_sos_options(obj);
            Xsupp.ineq = [Xsupp.ineq; obj.get_Z()];
        end
                
        function Xsupp = get_X_init(obj)
            Xsupp = get_X_init@loc_sos_options(obj);
            
            if isnumeric(Xsupp)
                %a kluge to deal with numeric data
                X0 = Xsupp;
                Xsupp = {X0; struct('ineq', obj.get_Z(), 'eq', [])};
            else
                Xsupp.ineq = [Xsupp.ineq; obj.get_Z()];
            end
        end
        
        function Xsupp = get_X_term(obj)
            %terminal set with x
            Xsupp = obj.X_term;
            if isempty(Xsupp)
                Xsupp = obj.X;
            end
            Xsupp.ineq = [Xsupp.ineq; obj.get_Z()];
        end
        
        function Z = get_Z(obj)
            z = obj.z;
            Z = z*(obj.Zmax - z);
        end
    end
end

